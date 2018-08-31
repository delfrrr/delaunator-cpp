#pragma once

#include <string>
#include <vector>

namespace json_helpers {
    std::string read_file(const char* filename);
    std::vector<double> get_geo_json_points(const std::string& json);
    std::vector<double> get_array_points(const std::string& json);
}
