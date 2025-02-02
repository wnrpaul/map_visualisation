# hike_mapper/__init__.py

# Import key classes and functions from your modules
from .data_loader import DataLoader
from .visualization import StaticMapRenderer
from .utils import calculate_zoom_level, calculate_cumulative_distance
from .plotter import ModulePlotter
from .styler import ModuleStyler
from .styles import STYLE_CONFIG

# Define what gets imported when someone uses `from hike_mapper import *`
__all__ = [
    'DataLoader',
    'StaticMapRenderer',
    'ModulePlotter',
    'ModuleStyler',
    'STYLE_CONFIG',
]

# Optional: Package metadata
__version__ = "1.0.0"
__author__ = "Paul Werner"
__email__ = "wnrpaul@gmail.com"
__licence__ = "MIT"
__url__ = "https://github.com/wnrpaul/hike-mapper"
__copyright__ = "Copyright 2025, Werner"
__status__ = "Development"
__keyword__ = "hiking, mapping, visualization"
__classifiers__ = [
    "Development Status :: 4 - Beta",
    "Intended Audience :: Developers",
    "License :: OSI Approved :: MIT License",
    "Programming Language :: Python :: 3.12",
]
__requires__ = [
    "geopandas",
    "matplotlib",
    "numpy",
    "overpy",
    "requests",
    "shapely",
    "folium",
    "gpxpy",
    "contextily",
    "osmnx",
    "geopy",
]
__python_requires__ = ">=3.12"
__description__ = "A tool for visualizing hiking tracks with OSM data and elevation profiles."


# hike_mapper/__init__.py
#If your package has many modules and you want to optimize performance, you can use lazy loading in __init__.py to only import modules when they are actually used. For example:
#def __getattr__(name):
#    """Lazy load modules only when needed."""
#    if name == "DataLoader":
#        from .data_loader import DataLoader
#        return DataLoader
#    elif name == "MapVisualizer":
#        from .visualization import MapVisualizer
#        return MapVisualizer
#    raise AttributeError(f"module 'hike_mapper' has no attribute '{name}'")
#This approach is useful for large packages where importing everything at once might slow down startup time.
