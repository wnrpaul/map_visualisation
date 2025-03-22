from mapper import Mapper

if __name__ == "__main__":
    # Initialize the Mapper
    app = Mapper()

    # Path to your GPS file
    gps_file_path = "/Users/werner/Downloads/Course_à_pied_dans_l_après_midi.gpx"

    app.load_gps_file(gps_file_path)

    app.load_osm_features(distance=1000)

    #app.export_map_objects_to_geojson(file_path="osm_features.geojson")

    app.style_map_objects(file_path="styles/default.yaml")

    # Generate map plot
    app.generate_map_plot(output_file="static_map.png") 

    # Generate classical OSM map
    app.generate_osm_map(output_file="osm_classical_map.png")
