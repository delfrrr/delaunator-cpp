#define CATCH_CONFIG_MAIN
#include "delaunator.h"
#include "catch.hpp"
#include "json-helpers.h"
#include <string>
#include <vector>

using namespace std;

TEST_CASE("triangles match JS version ouput", "[Delaunator]") {
    string points_str = json_helpers::read_file("./test-files/playgrounds-1356-epsg-3857.geojson");
    string triangles_str = json_helpers::read_file("./test-files/playgrounds-1356-triangles.json");
    const vector<double> coords = json_helpers::get_geo_json_points(points_str);
    const vector<double> triangles = json_helpers::get_array_points(triangles_str);
    Delaunator delaunator(move(coords));

    SECTION("length of triangles is the same") {
        REQUIRE(delaunator.triangles.size() == triangles.size());
    }

    SECTION("values are the same") {
        for(size_t i = 0; i < triangles.size(); i++)
        {
            REQUIRE(delaunator.triangles[i] == triangles[i]);
        }
    }
}

