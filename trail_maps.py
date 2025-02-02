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
import contextily as ctx

class HikeMapper:
    def __init__(self):
        self.gps_data = None
        self.map_center = None
        self.track_coordinates = []
        self.map_objects = {}  # Store map objects (e.g., trees, buildings)


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
    
    # Generate map plot
    app.generate_map_plot(output_file="static_map.png")

    # Generate classical OSM map
    app.generate_osm_map(output_file="osm_classical_map.png")
