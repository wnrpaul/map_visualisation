from geopy.distance import geodesic
from shapely.geometry import box

from math import log2


def calculate_zoom_level(bbox):
      """Calculate a valid zoom level based on the bounding box."""
      
      # Approximate Earth's circumference in meters
      earth_circumference = 40075016.686
      
      # Calculate the width and height of the bounding box in meters
      width = bbox[2] - bbox[0]
      height = bbox[3] - bbox[1]
      
      # Calculate the zoom level based on the larger dimension
      max_dimension = max(width, height)
      zoom = int(log2(earth_circumference / max_dimension))
      
      # Ensure the zoom level is within valid bounds (0-19)
      return max(0, min(zoom, 19))


def calculate_cumulative_distance(self, coordinates):
      """Calculate cumulative distance along a GPS track."""
      cumulative_distance = [0]  # Start with 0 distance
      for i in range(1, len(coordinates)):
          prev_coord = coordinates[i - 1]
          curr_coord = coordinates[i]
          segment_distance = geodesic(prev_coord, curr_coord).meters  # Distance in meters
          cumulative_distance.append(cumulative_distance[-1] + segment_distance)
      return cumulative_distance


def create_bounding_box(lon, lat, buffer=0.01):
    return box(lon - buffer, lat - buffer, lon + buffer, lat + buffer)
