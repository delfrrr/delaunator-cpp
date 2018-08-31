
#include "json-helpers.h"
#include <fstream>
#include <exception>
#include "rapidjson/document.h"

using namespace std;

string json_helpers::read_file(const char* filename) {
    ifstream input_file(filename);
    if(input_file.good()) {
        string json_str(
            (istreambuf_iterator<char>(input_file)),
            istreambuf_iterator<char>()
        );
        return json_str;
    } else {
        printf("Error reading file %s", filename);
        throw runtime_error("Error reading file");
    }
}

vector<double> json_helpers::get_geo_json_points(const string& json) {
    rapidjson::Document document;
    if(document.Parse(json.c_str()).HasParseError()) {
        throw runtime_error("Cannot parse JSON");
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
    }
    return coords;
}

vector<double> json_helpers::get_array_points(const string& json) {
    vector<double> points;
    rapidjson::Document document;
    if(document.Parse(json.c_str()).HasParseError()) {
        throw runtime_error("Cannot parse JSON");
    }
    if(!document.IsArray()) {
        throw runtime_error("It's not JSON Array");
    }
    points.reserve(static_cast<int>(document.Size()));
    for(rapidjson::SizeType i = 0; i < document.Size(); i++) {
        points.push_back(document[i].GetDouble());
    }
    return points;

}
