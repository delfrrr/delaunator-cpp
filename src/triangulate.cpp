#include "rapidjson/document.h"
#include "rapidjson/prettywriter.h"
#include "delaunator.h"
#include "json-helpers.h"
#include <cstdio>
#include <fstream>
#include <vector>
#include <initializer_list>
// #include "prettyprint.hpp"
#include <iostream>
using namespace std;

namespace {

    const string serialize_to_json(Delaunator &delaunator) {
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
                    for(long int i = 0; i < delaunator.triangles.size(); i+=3) {
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
                                                writer.Uint(delaunator.coords[2 * delaunator.triangles[i]]);
                                                writer.Uint(delaunator.coords[2 * delaunator.triangles[i] + 1]);
                                            writer.EndArray();
                                            writer.StartArray();
                                                writer.Uint(delaunator.coords[2 * delaunator.triangles[i + 1]]);
                                                writer.Uint(delaunator.coords[2 * delaunator.triangles[i + 1] + 1]);
                                            writer.EndArray();
                                            writer.StartArray();
                                                writer.Uint(delaunator.coords[2 * delaunator.triangles[i + 2]]);
                                                writer.Uint(delaunator.coords[2 * delaunator.triangles[i + 2] + 1]);
                                            writer.EndArray();
                                            writer.StartArray();
                                                writer.Uint(delaunator.coords[2 * delaunator.triangles[i]]);
                                                writer.Uint(delaunator.coords[2 * delaunator.triangles[i] + 1]);
                                            writer.EndArray();
                                        writer.EndArray();
                                        writer.EndArray();
                                writer.EndObject();
                        writer.EndObject();
                    }
                writer.EndArray();
        writer.EndObject();
        return string(sb.GetString());
    }
}


int main(int, char* argv[]) {
    const char* filename = argv[1];
    const char* output = argv[2];
    string json = json_helpers::read_file(filename);
    const vector<double> coords = json_helpers::get_geo_json_points(json);
    Delaunator delaunator(move(coords));
    const char* out_json = serialize_to_json(delaunator).c_str();

    if (output) {
        printf("Writing to file %s", output);
        ofstream stream;
        stream.open(output);
        stream << out_json;
        stream.close();
    } else {
        puts(out_json);
    }
    return 0;
}
