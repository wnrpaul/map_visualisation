# data_loader.py
import gpxpy
import geopandas as gpd
import numpy as np
import requests
import math
import overpy
import os
from shapely.geometry import LineString, Point, Polygon, MultiLineString, MultiPolygon, box
from utils import create_bounding_box

DEFAULT_BUFFER = 0.01
MIN_DISTANCE = 100
MAX_DISTANCE = 10000


class DataLoader:
  def __init__(self):
    self.coordinates = []
    self.elevations = []
    self.map_center = None
    self.map_objects = {}

  def load_gps_file(self, file_path):
    """Load and parse GPS track file (GPX, KML, GeoJSON)."""
    if not os.path.exists(file_path):
      raise FileNotFoundError(f"File not found: {file_path}")

    try:
      if file_path.endswith('.gpx'):
        with open(file_path, 'r') as gpx_file:
          gpx = gpxpy.parse(gpx_file)
          for track in gpx.tracks:
            for segment in track.segments:
              for point in segment.points:
                self.coordinates.append(
                    (point.latitude, point.longitude))
      elif file_path.endswith('.geojson'):
        gdf = gpd.read_file(file_path)
        self.coordinates = list(zip(gdf.geometry.y, gdf.geometry.x))
      elif file_path.endswith('.kml'):
        gdf = gpd.read_file(file_path, driver='KML')
        self.coordinates = list(zip(gdf.geometry.y, gdf.geometry.x))
      else:
        raise ValueError(
            "Unsupported file format. Please use GPX, KML, or GeoJSON.")

      self.map_center = (np.mean([coord[0] for coord in self.coordinates]),
                         np.mean([coord[1] for coord in self.coordinates]))
      print("GPS data loaded successfully.")
    except Exception as e:
      print(f"Error loading GPS file: {e}")

  def fetch_elevation_data(self, coordinates):
    """Fetch elevation data using the Open-Elevation API."""
    url = "https://api.open-elevation.com/api/v1/lookup"
    locations = [{"latitude": lat, "longitude": lon}
                 for lat, lon in coordinates]
    payload = {"locations": locations}

    try:
      response = requests.post(url, json=payload)
      response.raise_for_status()
      data = response.json()
      self.elevations = [result["elevation"] for result in data["results"]]
    except Exception as e:
      print(f"Error fetching elevation data: {e}")
      self.elevations = None

  def load_osm_features(self, distance=None):
    """
    Load all OpenStreetMap features within the specified radius using Overpass API.
    
    Parameters:
        distance (int): Radius around the map center to fetch features (in meters).
    """
    if not self.map_center:
        print("No GPS data available. Please upload a GPS track first.")
        return
    
    # Calculate default distance if not provided
    if distance is None:
        min_lon = min(lon for _, lon in self.coordinates)
        max_lon = max(lon for _, lon in self.coordinates)
        min_lat = min(lat for lat, _ in self.coordinates)
        max_lat = max(lat for lat, _ in self.coordinates)
        # Approximate distance as half the diagonal of the bounding box
        distance = int(
            math.sqrt((max_lon - min_lon)**2 + (max_lat - min_lat)**2) * 111000
        )
        print(f"Automatically calculated distance: {distance} meters.")
    
    if distance < MIN_DISTANCE or distance > MAX_DISTANCE:
        raise ValueError(f"Distance must be between {MIN_DISTANCE} and {MAX_DISTANCE} meters.")
    
    try:
        api = overpy.Overpass()
        lat, lon = self.map_center
        
        # Build the Overpass query to fetch all objects within the radius
        query = f"""
        [out:json];
        (
            node(around:{distance},{lat},{lon});
            way(around:{distance},{lat},{lon});
            relation(around:{distance},{lat},{lon});
        );
        out body;
        >;
        out skel qt;
        """
        # Debugging: Print the query
        print(f"Executing Overpass query:\n{query}")
        
        # Fetch data
        result = api.query(query)
        
        # Debugging: Print the number of nodes, ways, and relations returned
        print(f"Nodes: {len(result.nodes)}, Ways: {len(result.ways)}, Relations: {len(result.relations)}")
        
        # Process nodes, ways, and relations
        geometries_by_tag = {}
        for node in result.nodes:
            point = Point(float(node.lon), float(node.lat))
            for key, value in node.tags.items():
                tag_key = f"{key}={value}"
                if tag_key not in geometries_by_tag:
                    geometries_by_tag[tag_key] = []
                geometries_by_tag[tag_key].append(point)
        
        for way in result.ways:
            coords = [(float(node.lon), float(node.lat)) for node in way.nodes]
            if len(coords) > 1:  # Ensure there are at least two points to form a line
                line = LineString(coords)
                for key, value in way.tags.items():
                    tag_key = f"{key}={value}"
                    if tag_key not in geometries_by_tag:
                        geometries_by_tag[tag_key] = []
                    geometries_by_tag[tag_key].append(line)
        
        for relation in result.relations:
            outer_coords = []
            for member in relation.members:
                if member.role == "outer" and hasattr(member, "way"):
                    coords = [(float(node.lon), float(node.lat)) for node in member.way.nodes]
                    outer_coords.append(coords)
            
            if outer_coords:
                polygons = [Polygon(coords) for coords in outer_coords]
                multi_polygon = MultiPolygon(polygons) if len(polygons) > 1 else polygons[0]
                for key, value in relation.tags.items():
                    tag_key = f"{key}={value}"
                    if tag_key not in geometries_by_tag:
                        geometries_by_tag[tag_key] = []
                    geometries_by_tag[tag_key].append(multi_polygon)
        
        # Create GeoDataFrames for each tag
        self.map_objects = {}
        for tag_key, geometries in geometries_by_tag.items():
            gdf = gpd.GeoDataFrame(geometry=geometries, crs="EPSG:4326")
            self.map_objects[tag_key] = gdf
            print(f"Loaded {len(gdf)} features with tag '{tag_key}' from OSM.")
    except Exception as e:
        print(f"Error loading OSM features: {e}")

  def export_map_objects_to_geojson(self, output_dir="output"):
    """
    Export all map objects to GeoJSON files.
    
    Parameters:
        output_dir (str): Directory where the GeoJSON files will be saved.
    """
    import os
    os.makedirs(output_dir, exist_ok=True)  # Create the output directory if it doesn't exist
    
    for tag_key, gdf in self.map_objects.items():
        # Sanitize the tag key to create a valid filename
        sanitized_tag = tag_key.replace("=", "_").replace(":", "_")
        file_path = os.path.join(output_dir, f"{sanitized_tag}.geojson")
        
        # Save the GeoDataFrame as a GeoJSON file
        try:
            gdf.to_file(file_path, driver="GeoJSON")
            print(f"Exported {len(gdf)} features with tag '{tag_key}' to {file_path}")
        except Exception as e:
            print(f"Error exporting tag '{tag_key}': {e}")
