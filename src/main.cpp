// #include <iostream>
#include "rapidjson/document.h"
#include "delaunator.h"
#include <cstdio>
#include <fstream>
#include <string>
#include <exception>
#include <vector>
#include <initializer_list>
using namespace std;

namespace {
    string read_file(const char* filename) {
        ifstream input_file(filename);
        if(input_file.good()) {
            string json_str(
                (istreambuf_iterator<char>(input_file)),
                istreambuf_iterator<char>()
            );
            return json_str;
        } else {
            printf("Error reading file %s", filename);
            throw exception();
        }
    }

    vector<double> get_geo_json_points(const string& json) {
        rapidjson::Document document;
        if(document.Parse(json.c_str()).HasParseError()) {
            fprintf(stderr, "Cannot parse JSON");
            throw exception();
        }
        const rapidjson::Value& features = document["features"];
        vector<double> coords;
        // vector<double> y_vector;
        for(rapidjson::SizeType i = 0; i < features.Size(); i++) {
            const rapidjson::Value& coordinates = features[i]["geometry"]["coordinates"];
            const double x = coordinates[0].GetDouble();
            const double y = coordinates[1].GetDouble();
            coords.push_back(x);
            coords.push_back(y);
            // printf("coordinates %f %f \n", x, y);
        }
        // Points points = {.x_vector = x_vector, .y_vector = y_vector};
        return coords;
    }
}


int main(int, char* argv[]) {
    const char* filename = argv[1];
    string json = read_file(filename);
    const vector<double> coords = get_geo_json_points(json);
    // vector<double> xy_points;
    // xy_points.reserve(points.x_vector.size() + points.y_vector.size());
    // xy_points.insert(xy_points.end(), points.x_vector.begin(), points.x_vector.end());
    // xy_points.insert(xy_points.end(), points.y_vector.begin(), points.y_vector.end());
    // const int size = points.x_vector.size();
    // double coords[size * 2];
    // for(int i = 0; i < size; i++)
    // {
    //     coords[i] = points.x_vector[i];
    //     coords[2 * i] = points.y_vector[i];
    // }
    Delaunator delaunator(coords);
}
