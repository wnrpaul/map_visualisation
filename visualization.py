# visualization.py
import contextily as ctx
import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.geometry import LineString, Point, Polygon, MultiLineString, box


class StaticMapRenderer:
  def set_coordinates(self, coordinates):
    """
    Set the coordinates.
    
    Parameters:
        map_objects (dict): A dictionary containing GeoDataFrames for different
        feature types.
    """
    self.coordinates = coordinates

  def generate_osm_map(self, output_file="osm_map.png"):
    """Generate the OpenStreetMap-style map with GPS track overlay."""
    if not self.coordinates:
        print("No GPS data available. Please upload a GPS track first.")
        return

    # Create track geometry
    track_geom = LineString([(lon, lat) for lat, lon in self.coordinates])
    gdf_track = gpd.GeoDataFrame(geometry=[track_geom], crs="EPSG:4326")

    # Re-project to Web Mercator
    gdf_track = gdf_track.to_crs(epsg=3857)

    # Calculate bounds with proper buffer (500 meters)
    min_x, min_y, max_x, max_y = gdf_track.total_bounds
    buffer = 500  # Meters
    bbox = (
        min_x - buffer,
        min_y - buffer,
        max_x + buffer,
        max_y + buffer
    )

    # Create plot
    fig, ax = plt.subplots(figsize=(10, 10))

    # Plot invisible track to set bounds
    gdf_track.plot(ax=ax, alpha=0)

    # Add basemap with automatic zoom detection
    ctx.add_basemap(
        ax,
        source=ctx.providers.OpenStreetMap.Mapnik,
        crs=gdf_track.crs.to_string()
    )

    # Set axis limits and remove borders
    ax.set_xlim(bbox[0], bbox[2])
    ax.set_ylim(bbox[1], bbox[3])
    ax.set_axis_off()

    # Add visible track overlay
    gdf_track.plot(ax=ax, color='red', linewidth=2)

    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"OSM Reference Map saved to {output_file}")


# class InteractiveMapRenderer:

#    def generate_interactive_map(self, output_file="interactive_map.html"):
#        """Generate an interactive map using Folium."""
#        if not self.track_coordinates:
#            print("No GPS data available. Please upload a GPS track first.")
#            return

#        # Create a Folium map centered on the track
#        m = folium.Map(location=self.map_center, zoom_start=13, tiles="OpenStreetMap")

#        # Add the hiking track
#        folium.PolyLine(locations=self.track_coordinates, color="red", weight=5, opacity=0.7).add_to(m)

#        # Add styled map objects
#        for feature_type, gdf in self.map_objects.items():
#            if "color" in gdf.columns:
#                # Limit to 100 markers for performance
#                for _, row in gdf.iloc[:100].iterrows():
#                    folium.CircleMarker(
#                        location=(row.geometry.y, row.geometry.x),
#                        radius=row.get("size", 5),
#                        color=row["color"],
#                        fill=True,
#                        fill_color=row["color"]
#                    ).add_to(m)

#        # Add styled map objects
#        for feature_type, gdf in self.map_objects.items():
#            if "color" in gdf.columns:
#                # Create a heatmap layer
#                heat_data = [[row.geometry.y, row.geometry.x] for _, row in gdf.iterrows()]
#                HeatMap(heat_data).add_to(m)

#        # Save the map
#        m.save(output_file)
#        print(f"Interactive map saved as {output_file}")
