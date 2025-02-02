# styler.py
import yaml


class ModuleStyler:
  def __init__(self):
    self.map_objects = {}
    self.styles = {}

  def load_styles_from_yaml(self, file_path):
    """
    Load styles from a configuration file or dictionary.
    
    Parameters:
        styles_config (dict): A dictionary containing style definitions.
    """
    with open(file_path, 'r') as f:
        self.styles = yaml.safe_load(f)

  def set_map_objects(self, map_objects):
    """
    Set the map_objects dictionary.

    Parameters:
      map_objects (dict): A dictionary containing GeoDataFrames for different feature types.
    """
    self.map_objects = map_objects

  def style_map_objects(self):
    """
    Apply custom styles to all map objects based on the loaded styles.
    """
    if not self.map_objects:
      print("No map objects available. Load the data first using `load_osm_features`.")
      return
        
    for feature_type, gdf in self.map_objects.items():
      if feature_type in self.styles:
      
        style = self.styles[feature_type]
        gdf["color"] = style.get("color", "black")  # Default to black if no color is specified
        if "size" in style:
            gdf["size"] = style["size"]
        if "linewidth" in style:
            gdf["linewidth"] = style["linewidth"]
        if "alpha" in style:
            gdf["alpha"] = style["alpha"]
        print(f"Styled {len(gdf)} features with tag '{feature_type}' using style: {style}")
      else:
        #print(f"No style defined for feature type '{feature_type}'. Using default style.")
        gdf["color"] = "black"
        gdf["linewidth"] = 1
