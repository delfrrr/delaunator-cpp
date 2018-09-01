
#include <chrono>
#include "delaunator.h"
#include "json-helpers.h"
#include <string>
#include <vector>

using namespace std;
int main(int, char* argv[]) {
    string points_str = json_helpers::read_file("./test-files/osm-nodes-45331-epsg-3857.geojson");
    const vector<double> coords = json_helpers::get_geo_json_points(points_str);

    auto t_start = chrono::high_resolution_clock::now();
    Delaunator delaunator(move(coords));
    auto t_end = chrono::high_resolution_clock::now();

    auto milliseconds = chrono::duration_cast<chrono::milliseconds>(t_end - t_start).count();

    printf("coords=%lu \n", delaunator.coords.size() / 2);
    printf("milliseconds=%lld \n", milliseconds);
    printf("triangles=%lu \n", delaunator.triangles.size());

    return 0;
}