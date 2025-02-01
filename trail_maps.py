import os
import geopandas as gpd
import matplotlib.pyplot as plt
import folium
from folium.plugins import HeatMap, Fullscreen
import gpxpy
import requests
import numpy as np
from shapely.geometry import LineString, Point, Polygon, MultiLineString, box
import osmnx as ox
import overpy
import math
from geopy.distance import geodesic

class HikeMapper:
    def __init__(self):
        self.gps_data = None
        self.map_center = None
        self.track_coordinates = []
        self.map_objects = {}  # Store map objects (e.g., trees, buildings)

    def calculate_cumulative_distance(self, coordinates):
      """Calculate cumulative distance along a GPS track."""
      cumulative_distance = [0]  # Start with 0 distance
      for i in range(1, len(coordinates)):
          prev_coord = coordinates[i - 1]
          curr_coord = coordinates[i]
          segment_distance = geodesic(prev_coord, curr_coord).meters  # Distance in meters
          cumulative_distance.append(cumulative_distance[-1] + segment_distance)
      return cumulative_distance
    
    def plot_elevation_profile(self, output_file="elevation_profile.png"):
      """Plot elevation profile alongside the static map."""
      if not self.track_coordinates:
          print("No GPS data available. Please upload a GPS track first.")
          return
      
      # Fetch elevation data
      elevations = self.fetch_elevation_data(self.track_coordinates)
      if elevations is None:
          print("Failed to fetch elevation data. Skipping elevation profile.")
          return
      
      # Calculate cumulative distance
      cumulative_distance = self.calculate_cumulative_distance(self.track_coordinates)
      
      # Plot elevation profile
      fig, ax = plt.subplots(figsize=(10, 4))
      ax.plot(cumulative_distance, elevations, color="blue", linewidth=2)
      ax.set_title("Elevation Profile")
      ax.set_xlabel("Distance (meters)")
      ax.set_ylabel("Elevation (meters)")
      ax.grid(True)
      
      # Save the elevation profile
      plt.tight_layout()
      plt.savefig(output_file, dpi=300)
      print(f"Elevation profile saved as {output_file}")



    def load_osm_features_overpass(self, feature_type="trees", tags=None, distance=None):
      """
      Load OpenStreetMap features using Overpass API.
      
      Parameters:
          feature_type (str): A name for the feature type (e.g., "trees", "buildings").
          tags (dict): OSM tags to filter features (e.g., {"natural": "tree"}).
          distance (int): Radius around the map center to fetch features (in meters).
      """
      if not self.map_center:
          print("No GPS data available. Please upload a GPS track first.")
          return
      
      if tags is None:
          print("No tags provided. Please specify OSM tags to filter features.")
          return

      # Calculate default distance if not provided
      if distance is None:
        min_lon = min(lon for _, lon in self.track_coordinates)
        max_lon = max(lon for _, lon in self.track_coordinates)
        min_lat = min(lat for lat, _ in self.track_coordinates)
        max_lat = max(lat for lat, _ in self.track_coordinates)
        
        # Approximate distance as half the diagonal of the bounding box
        distance = int(
            math.sqrt((max_lon - min_lon)**2 + (max_lat - min_lat)**2) * 111000 / 2
        )
        print(f"Automatically calculated distance: {distance} meters.")

      if distance < 100 or distance > 10000:
        raise ValueError("Distance must be between 100 and 10,000 meters.")
      
      try:
          api = overpy.Overpass()
          lat, lon = self.map_center
          
          # Build the Overpass query dynamically based on the provided tags
          tag_filters = "".join([f'["{key}"="{value}"]' for key, value in tags.items()])
          query = f"""
          [out:json];
          (
              node{tag_filters}(around:{distance},{lat},{lon});
              way{tag_filters}(around:{distance},{lat},{lon});
              relation{tag_filters}(around:{distance},{lat},{lon});
          );
          out body;
          >;
          out skel qt;
          """

          print(f"Executing Overpass query:\n{query}")  # Debugging: Print the query
    
          # Fetch data
          result = api.query(query)
          
          # Debugging: Print the number of nodes, ways, and relations returned
          print(f"Nodes: {len(result.nodes)}, Ways: {len(result.ways)}, Relations: {len(result.relations)}")
          
          # Calculate the bounding box for clipping
          buffer = 0.01  # Adjust this value as needed (in degrees)
          bbox_polygon = box(lon - buffer, lat - buffer, lon + buffer, lat + buffer)
          
          # Convert nodes, ways, and relations to geometries
          geometries = []
          for node in result.nodes:
              point = Point(float(node.lon), float(node.lat))
              if bbox_polygon.contains(point):  # Only include points within the bounding box
                geometries.append(point)

          for way in result.ways:
            coords = [(float(node.lon), float(node.lat)) for node in way.nodes]
            if len(coords) > 1:  # Ensure there are at least two points to form a line
                line = LineString(coords)
                clipped_line = line.intersection(bbox_polygon)  # Clip the line to the bounding box
                if clipped_line.is_empty:
                    continue
                if clipped_line.geom_type == "MultiLineString":
                    geometries.extend([LineString(segment) for segment in clipped_line.geoms])
                else:
                    geometries.append(clipped_line)

          for relation in result.relations:
              outer_coords = []
              for member in relation.members:
                  if member.role == "outer" and hasattr(member, "way"):
                      coords = [(float(node.lon), float(node.lat)) for node in member.way.nodes]
                      outer_coords.append(coords)
              
              if outer_coords:
                  polygons = [Polygon(coords) for coords in outer_coords]
                  geometries.append(MultiPolygon(polygons) if len(polygons) > 1 else polygons[0])
          
          # Create a GeoDataFrame
          if not geometries:
              print(f"No {feature_type} found within the specified radius.")
              return
          
          gdf = gpd.GeoDataFrame(geometry=geometries, crs="EPSG:4326")
          
          # Store the loaded features
          self.map_objects[feature_type] = gdf
          print(f"Loaded {len(gdf)} {feature_type} features from OSM.")
      except Exception as e:
          print(f"Error loading OSM features: {e}")

    def load_osm_features(self, feature_type="trees", distance=1000):
        """Load OpenStreetMap features near the track."""
        if not self.map_center:
            print("No GPS data available. Please upload a GPS track first.")
            return
        
        try:
            # Fetch OSM data within a specified radius around the map center
            tags = {"natural": "tree"} if feature_type == "trees" else {"building": True}
            gdf = ox.geometries_from_point(self.map_center, tags, dist=distance)
            
            # Store the loaded features
            self.map_objects[feature_type] = gdf
            print(f"Loaded {len(gdf)} {feature_type} features from OSM.")
        except Exception as e:
            print(f"Error loading OSM features: {e}")

    def style_map_objects(self, feature_type="trees", color="green", size=5, linewidth=None):
        """
        Apply custom styles to map objects.
        
        Parameters:
            feature_type (str): The type of feature to style (e.g., "trees", "streets").
            color (str): The color to apply to the feature.
            size (int): The size for point features (e.g., trees).
            linewidth (int): The line width for line or polygon features (e.g., streets, rivers, buildings).
        """
        if feature_type not in self.map_objects:
            print(f"No {feature_type} data available. Load the data first using `load_osm_features`.")
            return
        
        gdf = self.map_objects[feature_type]
        
        # Apply custom styles
        gdf["color"] = color
        if size is not None:
            gdf["size"] = size  # For point features
        if linewidth is not None:
            gdf["linewidth"] = linewidth  # For line or polygon features
        
        print(f"Styled {len(gdf)} {feature_type} with color={color}, size={size}, linewidth={linewidth}.")

    def load_gps_file(self, file_path):
        """Load and parse GPS track file (GPX, KML, GeoJSON)."""
        try:
            if file_path.endswith('.gpx'):
                with open(file_path, 'r') as gpx_file:
                    gpx = gpxpy.parse(gpx_file)
                    for track in gpx.tracks:
                        for segment in track.segments:
                            for point in segment.points:
                                self.track_coordinates.append((point.latitude, point.longitude))
            elif file_path.endswith('.geojson'):
                gdf = gpd.read_file(file_path)
                self.track_coordinates = list(zip(gdf.geometry.y, gdf.geometry.x))
            elif file_path.endswith('.kml'):
                gdf = gpd.read_file(file_path, driver='KML')
                self.track_coordinates = list(zip(gdf.geometry.y, gdf.geometry.x))
            else:
                raise ValueError("Unsupported file format. Please use GPX, KML, or GeoJSON.")
            
            self.map_center = (np.mean([coord[0] for coord in self.track_coordinates]),
                               np.mean([coord[1] for coord in self.track_coordinates]))
            print("GPS track loaded successfully!")
        except Exception as e:
            print(f"Error loading GPS file: {e}")

    def generate_static_map(self, output_file="map.png", color_scheme="viridis"):
      """Generate a static map with customizable styles."""
      if not self.track_coordinates:
          print("No GPS data available. Please upload a GPS track first.")
          return
      
      # Create a GeoDataFrame for the track
      track_geom = LineString([(lon, lat) for lat, lon in self.track_coordinates])
      gdf_track = gpd.GeoDataFrame(geometry=[track_geom], crs="EPSG:4326")
      
      # Calculate the bounding box of the GPS track
      min_lon = min(lon for _, lon in self.track_coordinates)
      max_lon = max(lon for _, lon in self.track_coordinates)
      min_lat = min(lat for lat, _ in self.track_coordinates)
      max_lat = max(lat for lat, _ in self.track_coordinates)
      
      # Add a small buffer around the bounding box to ensure features are visible
      buffer = 0.01  # Adjust this value as needed (in degrees)
      bbox = (min_lon - buffer, min_lat - buffer, max_lon + buffer, max_lat + buffer)

      # Create a bounding box polygon for clipping
      bbox_polygon = box(bbox[0], bbox[1], bbox[2], bbox[3])
      
      # Plot the map
      fig, ax = plt.subplots(figsize=(10, 10))
      gdf_track.plot(ax=ax, color="red", linewidth=2, label="Hiking Track")
      
      # Plot styled map objects
      for feature_type, gdf in self.map_objects.items():
          if "color" in gdf.columns:
              # Clip geometries to the bounding box
              clipped_gdf = gpd.clip(gdf, bbox_polygon)
              
              if feature_type == "streets":
                  clipped_gdf.plot(ax=ax, color=clipped_gdf["color"], linewidth=clipped_gdf.get("linewidth", 1), label=feature_type)
              elif feature_type == "rivers":
                  clipped_gdf.plot(ax=ax, color=clipped_gdf["color"], linewidth=clipped_gdf.get("linewidth", 2), label=feature_type)
              elif feature_type == "buildings":
                  clipped_gdf.plot(ax=ax, facecolor=clipped_gdf["color"], edgecolor="black", alpha=0.7, linewidth=clipped_gdf.get("linewidth", 1), label=feature_type)
              else:  # Points (e.g., trees)
                  clipped_gdf.plot(ax=ax, color=clipped_gdf["color"], markersize=clipped_gdf.get("size", 5), label=feature_type)
      
      # Set the map's extent to the bounding box of the GPS track
      #ax.set_xlim(bbox[0], bbox[2])
      #ax.set_ylim(bbox[1], bbox[3])
        
      # Save the map
      plt.savefig(output_file, dpi=300)
      print(f"Static map saved as {output_file}")

    def fetch_elevation_data(self, coordinates):
        """Fetch elevation data using the Open-Elevation API."""
        url = "https://api.open-elevation.com/api/v1/lookup"
        locations = [{"latitude": lat, "longitude": lon} for lat, lon in coordinates]
        payload = {"locations": locations}
        
        try:
            response = requests.post(url, json=payload)
            response.raise_for_status()
            data = response.json()
            elevations = [result["elevation"] for result in data["results"]]
            return elevations
        except Exception as e:
            print(f"Error fetching elevation data: {e}")
            return None

    def generate_interactive_map(self, output_file="interactive_map.html"):
      """Generate an interactive map using Folium."""
      if not self.track_coordinates:
          print("No GPS data available. Please upload a GPS track first.")
          return
      
      # Create a Folium map centered on the track
      m = folium.Map(location=self.map_center, zoom_start=13, tiles="OpenStreetMap")
      
      # Add the hiking track
      folium.PolyLine(locations=self.track_coordinates, color="red", weight=5, opacity=0.7).add_to(m)
      
      # Add styled map objects
      for feature_type, gdf in self.map_objects.items():
          if "color" in gdf.columns:
              # Limit to 100 markers for performance
              for _, row in gdf.iloc[:100].iterrows():
                  folium.CircleMarker(
                      location=(row.geometry.y, row.geometry.x),
                      radius=row.get("size", 5),
                      color=row["color"],
                      fill=True,
                      fill_color=row["color"]
                  ).add_to(m)

      # Add styled map objects
      for feature_type, gdf in self.map_objects.items():
          if "color" in gdf.columns:
              # Create a heatmap layer
              heat_data = [[row.geometry.y, row.geometry.x] for _, row in gdf.iterrows()]
              HeatMap(heat_data).add_to(m)
            
      # Save the map
      m.save(output_file)
      print(f"Interactive map saved as {output_file}")


if __name__ == "__main__":
    app = HikeMapper()
    
    # Load GPS file
    app.load_gps_file('/Users/werner/Downloads/Course_à_pied_dans_l_après_midi.gpx')
    
    # Load trees
    app.load_osm_features_overpass(
        feature_type="trees",
        tags={"natural": "tree"},
        distance=1000
    )
    app.style_map_objects(feature_type="trees", color="green", size=5)
    
    # Load streets
    app.load_osm_features_overpass(
        feature_type="streets",
        tags={"highway": "residential"},  # Fetch residential streets
        distance=1000
    )
    app.style_map_objects(feature_type="streets", color="black", linewidth=1)
    
    # Load buildings
    app.load_osm_features_overpass(
        feature_type="buildings",
        tags={"building": "*"},  # Match any building type
        distance=1000
    )
    app.style_map_objects(feature_type="buildings", color="gray", linewidth=1)
    
    # Load rivers
    app.load_osm_features_overpass(
        feature_type="rivers",
        tags={"waterway": "river"},  # Fetch rivers
        distance=1000
    )
    app.style_map_objects(feature_type="rivers", color="blue", linewidth=2)
    
    # Generate static map
    static_output = "static_map.png"
    app.generate_static_map(output_file=static_output)
