// #include <iostream>
#include "rapidjson/document.h"
#include <cstdio>
#include <fstream>
#include <string>
#include <exception>
#include <vector>
#include <initializer_list>

namespace {
    std::string read_file(const char* filename) {
        std::ifstream input_file(filename);
        if(input_file.good()) {
            std::string json_str(
                (std::istreambuf_iterator<char>(input_file)),
                std::istreambuf_iterator<char>()
            );
            return json_str;
        } else {
            std::printf("Error reading file %s", filename);
            throw std::exception();
        }
    }

    struct Points {
        std::vector<double> x_vector;
        std::vector<double> y_vector;
    };

    const Points get_geo_json_points(const std::string& json) {
        rapidjson::Document document;
        if(document.Parse(json.c_str()).HasParseError()) {
            std::fprintf(stderr, "Cannot parse JSON");
            throw std::exception();
        }
        const rapidjson::Value& features = document["features"];
        std::vector<double> x_vector;
        std::vector<double> y_vector;
        for(rapidjson::SizeType i = 0; i < features.Size(); i++) {
            const rapidjson::Value& coordinates = features[i]["geometry"]["coordinates"];
            const double x = coordinates[0].GetDouble();
            const double y = coordinates[1].GetDouble();
            x_vector.push_back(std::move(x));
            y_vector.push_back(std::move(y));
            // std::printf("coordinates %f %f \n", x, y);
        }
        Points points = {.x_vector = x_vector, .y_vector = y_vector};
        return points;
    }
}


int main(int, char* argv[]) {
    const char* filename = argv[1];
    std::string json = read_file(filename);
    const Points points = get_geo_json_points(json);
    const int size = points.x_vector.size();
    double coords[size * 2];
    for(int i = 0; i < size; i++)
    {
        coords[i] = points.x_vector[i];
        coords[2 * i] = points.y_vector[i];
    }

}
