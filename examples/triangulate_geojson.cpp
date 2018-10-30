#include "rapidjson/document.h"
#include "rapidjson/prettywriter.h"
#include <delaunator.hpp>
#include "./utils.hpp"
#include <cstdio>
#include <fstream>
#include <vector>
#include <initializer_list>
#include <iostream>

std::string serialize_to_json(delaunator::DelaunatorResult const& delaunator, 
                              std::vector<utils::point> const& coords) {
    rapidjson::StringBuffer sb;
    rapidjson::PrettyWriter<rapidjson::StringBuffer> writer(sb);
    writer.StartObject();
        writer.String("type"); writer.String("FeatureCollection");
        writer.String("crs");
            writer.StartObject();
                writer.String("type"); writer.String("name");
                writer.String("properties");
                    writer.StartObject();
                        writer.String("name"); writer.String("urn:ogc:def:crs:EPSG::900913");
                    writer.EndObject();
            writer.EndObject();
        writer.String("features");
            writer.StartArray();
                for(std::size_t i = 0; i < delaunator.triangles.size(); i+=3) {
                    writer.StartObject();
                        writer.String("type"); writer.String("Feature");
                        writer.String("properties"); writer.StartObject(); writer.EndObject();
                        writer.String("geometry");
                            writer.StartObject();
                                writer.String("type"); writer.String("Polygon");
                                writer.String("coordinates");
                                    writer.StartArray();
                                    writer.StartArray();
                                        writer.StartArray();
                                            writer.Double(coords[delaunator.triangles[i]].x);
                                            writer.Double(coords[delaunator.triangles[i]].y);
                                        writer.EndArray();
                                        writer.StartArray();
                                            writer.Double(coords[delaunator.triangles[i + 1]].x);
                                            writer.Double(coords[delaunator.triangles[i + 1]].y);
                                        writer.EndArray();
                                        writer.StartArray();
                                            writer.Double(coords[delaunator.triangles[i + 2]].x);
                                            writer.Double(coords[delaunator.triangles[i + 2]].y);
                                        writer.EndArray();
                                        writer.StartArray();
                                            writer.Double(coords[delaunator.triangles[i]].x);
                                            writer.Double(coords[delaunator.triangles[i]].y);
                                        writer.EndArray();
                                    writer.EndArray();
                                    writer.EndArray();
                            writer.EndObject();
                    writer.EndObject();
                }
            writer.EndArray();
    writer.EndObject();
    return std::string(sb.GetString());
}

int main(int, char* argv[]) {
    const char* filename = argv[1];
    const char* output = argv[2];
    std::string json = utils::read_file(filename);
    auto coords = utils::get_geo_json_points(json);
    auto result = delaunator::delaunator(coords);
    const char* out_json = serialize_to_json(result, coords).c_str();

    if (output) {
        printf("Writing to file %s", output);
        std::ofstream stream;
        stream.open(output);
        stream << out_json;
        stream.close();
    } else {
        std::puts(out_json);
    }
    return 0;
}
