# plotter.py
import matplotlib.pyplot as plt
import geopandas as gpd
from shapely import LineString


class ModulePlotter:
  def __init__(self):
    self.coordinates = []
    self.map_objects = {}

  def set_map_objects(self, map_objects):
    """
    Set the map_objects dictionary.
    
    Parameters:
        map_objects (dict): A dictionary containing GeoDataFrames for different feature types.
    """
    self.map_objects = map_objects

  def set_coordinates(self, coordinates):
    """
    Set the coordinates.
    
    Parameters:
        map_objects (dict): A dictionary containing GeoDataFrames for different feature types.
    """
    self.coordinates = coordinates

  def plot_gps_trace(self, ax):
    """
    Plot the GPS trace on the given axis.
    
    Parameters:
        ax (matplotlib.axes.Axes): The axis on which to plot the GPS track.
    """
    if not self.coordinates:
      print("No GPS data available. Please upload a GPS track first.")
      return

    # Create a GeoDataFrame for the track
    track_geom = LineString([(lon, lat)
                            for lat, lon in self.coordinates])
    gdf_track = gpd.GeoDataFrame(geometry=[track_geom], crs="EPSG:4326")

    # Calculate the bounding box of the GPS track
    min_lon = min(lon for _, lon in self.coordinates)
    max_lon = max(lon for _, lon in self.coordinates)
    min_lat = min(lat for lat, _ in self.coordinates)
    max_lat = max(lat for lat, _ in self.coordinates)
    buffer = 0.01  # Adjust this value as needed (in degrees)
    bbox = (min_lon - buffer, min_lat - buffer,
            max_lon + buffer, max_lat + buffer)

    # Plot the track
    gdf_track.plot(ax=ax, color="red", linewidth=2)

    # Set map extent
    # ax.set_xlim(bbox[0], bbox[2])
    # ax.set_ylim(bbox[1], bbox[3])

  def plot_map_features(self, ax):
    """
    Plot map features on the given axis.
    
    Parameters:
        ax (matplotlib.axes.Axes): The axis on which to plot the map features.
    """
    if not self.map_objects:
        print("No map objects available. Please load OSM features first.")
        return

    # Plot other features
    for feature_type, gdf in self.map_objects.items():
      if "color" in gdf.columns:
        if feature_type == "streets":
          gdf.plot(ax=ax, color=gdf["color"],
                   linewidth=gdf.get("linewidth", 1))
        elif feature_type == "rivers":
          gdf.plot(ax=ax, color=gdf["color"],
                   linewidth=gdf.get("linewidth", 2))
        elif feature_type == "buildings":
          gdf.plot(ax=ax, facecolor=gdf["color"], edgecolor="black",
                   alpha=0.7, linewidth=gdf.get("linewidth", 1))
        else:  # Points (e.g., trees)
          gdf.plot(ax=ax, color=gdf["color"], markersize=gdf.get("size", 5))

  def generate_map_plot(self, output_file="map_plot.png"):
    """
    Generate a clean map plot without frames, axes, or titles.
    """
    if not self.coordinates:
      print("No GPS data available. Please upload a GPS track first.")
      return

    # Create a figure and axis
    fig, ax = plt.subplots(figsize=(10, 10))

    # Plot the GPS track
    self.plot_gps_trace(ax)

    # Plot map features
    self.plot_map_features(ax)

    # Remove frame, axes, and labels
    ax.set_frame_on(False)  # Remove the frame around the plot
    ax.set_xticks([])       # Remove x-axis ticks
    ax.set_yticks([])       # Remove y-axis ticks
    ax.set_xticklabels([])  # Remove x-axis labels
    ax.set_yticklabels([])  # Remove y-axis labels
    ax.grid(False)          # Disable grid lines

    # Save the clean map plot
    plt.tight_layout(pad=0)  # Remove padding
    plt.savefig(output_file, dpi=300, bbox_inches='tight', pad_inches=0)
    print(f"Clean map plot saved as {output_file}")

  def plot_elevation_profile(self, output_file="elevation_profile.png"):
    """Plot elevation profile alongside the static map."""
    if not self.coordinates:
      print("No GPS data available. Please upload a GPS track first.")
      return

    # Fetch elevation data
    elevations = self.fetch_elevation_data(self.coordinates)
    if elevations is None:
        print("Failed to fetch elevation data. Skipping elevation profile.")
        return

    # Calculate cumulative distance
    cumulative_distance = self.calculate_cumulative_distance(
        self.coordinates)

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
