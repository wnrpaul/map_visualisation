from data_loader import DataLoader
from styler import ModuleStyler
from visualization import StaticMapRenderer
from plotter import ModulePlotter
#from .styles import STYLE_CONFIG
#from utils import calculate_zoom_level, calculate_cumulative_distance


class Mapper:
  def __init__(self):
    # Initialize components
    self.data_loader = DataLoader()
    self.map_renderer = StaticMapRenderer()
    self.plotter = ModulePlotter()
    self.styler = ModuleStyler()

    # Shared state
    self.gps_data = None
    self.map_center = None
    self.coordinates = []
    self.elevations = None
    self.map_objects = {}  # Store map objects (e.g., trees, buildings)

  def load_gps_file(self, file_path):
    """Load GPS data from a file."""
    self.gps_data = self.data_loader.load_gps_file(file_path)
    self.coordinates = self.data_loader.coordinates
    self.map_center = self.data_loader.map_center

  def load_osm_features(self, distance=None):
    """Load OSM features using Overpass API."""
    self.data_loader.load_osm_features(distance)
    self.map_objects.update(self.data_loader.map_objects)

  def export_map_objects_to_geojson(self, file_path="osm_features.geojson"):
    """Export OSM features to a GeoJSON file."""
    self.data_loader.export_map_objects_to_geojson(file_path)

  def style_map_objects(self, file_path):
    """Apply custom styles to map objects."""
    self.styler.load_styles_from_yaml(file_path)
    self.styler.set_map_objects(self.map_objects)
    self.styler.style_map_objects()

  def generate_map_plot(self, output_file="static_map.png"):
    self.plotter.set_map_objects(self.map_objects)
    self.plotter.set_coordinates(self.coordinates)
    self.plotter.generate_map_plot(output_file)

  def generate_osm_map(self, output_file="osm_classical_map.png"):
    self.map_renderer.set_coordinates(self.coordinates)
    self.map_renderer.generate_osm_map(output_file)
  
