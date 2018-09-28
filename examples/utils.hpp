#pragma once

#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>
#include "rapidjson/document.h"

namespace utils {

inline std::string read_file(const char* filename) {
    std::ifstream input_file(filename);
    if(input_file.good()) {
        std::string json_str(
            (std::istreambuf_iterator<char>(input_file)),
            std::istreambuf_iterator<char>()
        );
        return json_str;
    } else {
        printf("Error reading file %s", filename);
        throw std::runtime_error("Error reading file");
    }
}

inline std::vector< double> get_geo_json_points(std::string const& json) {
    rapidjson::Document document;
    if(document.Parse(json.c_str()).HasParseError()) {
        throw std::runtime_error("Cannot parse JSON");
    }
    const rapidjson::Value& features = document["features"];
    std::vector<double> coords;
    // vector<double> y_vector;
    for(rapidjson::SizeType i = 0; i < features.Size(); i++) {
        const rapidjson::Value& coordinates = features[i]["geometry"]["coordinates"];
        const double x = coordinates[0].GetDouble();
        const double y = coordinates[1].GetDouble();
        coords.push_back(x);
        coords.push_back(y);
    }
    return coords;
}

template <typename T>
inline T get_json_value(const rapidjson::Value &v);

template <>
inline double get_json_value(const rapidjson::Value &v) {
    return v.GetDouble();
}

template <>
inline size_t get_json_value(const rapidjson::Value &v) {
    return static_cast<size_t>(v.GetUint64());
}

template <typename T>
inline std::vector<T> get_array_points(std::string const& json) {
    std::vector<T> points;
    rapidjson::Document document;
    if(document.Parse(json.c_str()).HasParseError()) {
        throw std::runtime_error("Cannot parse JSON");
    }
    if(!document.IsArray()) {
        throw std::runtime_error("It's not JSON Array");
    }
    points.reserve(static_cast<std::size_t>(document.Size()));
    for(rapidjson::SizeType i = 0; i < document.Size(); i++) {
        // points.push_back(document[i].GetDouble());
        points.push_back(get_json_value<T>(document[i]));
    }
    return points;
}

} // end ns utils
