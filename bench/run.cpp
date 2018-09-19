#include "../examples/utils.hpp"
#include <benchmark/benchmark.h>
#include <delaunator.hpp>
#include <random>
#include <string>
#include <vector>

namespace {
std::vector<double> generate_uniform(size_t n) {
    std::vector<double> coords;
    std::srand(350);
    double norm = static_cast<double>(RAND_MAX) / 1e3;
    for (size_t i = 0; i < n; i++) {
        coords.push_back(static_cast<double>(std::rand()) / norm);
        coords.push_back(static_cast<double>(std::rand()) / norm);
    }

    return coords;
}

void BM_45K_geojson_nodes(benchmark::State& state) {
    std::string points_str = utils::read_file("./test/test-files/osm-nodes-45331-epsg-3857.geojson");
    std::vector<double> coords = utils::get_geo_json_points(points_str);

    while (state.KeepRunning()) {
        delaunator::Delaunator delaunator(coords);
    }
}

void BM_2K_uniform(benchmark::State& state) {
    std::vector<double> coords = generate_uniform(2000);
    while (state.KeepRunning()) {
        delaunator::Delaunator delaunator(coords);
    }
}

void BM_100K_uniform(benchmark::State& state) {
    std::vector<double> coords = generate_uniform(100000);
    while (state.KeepRunning()) {
        delaunator::Delaunator delaunator(coords);
    }
}

void BM_200K_uniform(benchmark::State& state) {
    std::vector<double> coords = generate_uniform(200000);
    while (state.KeepRunning()) {
        delaunator::Delaunator delaunator(coords);
    }
}

void BM_500K_uniform(benchmark::State& state) {
    std::vector<double> coords = generate_uniform(500000);
    while (state.KeepRunning()) {
        delaunator::Delaunator delaunator(coords);
    }
}

void BM_1M_uniform(benchmark::State& state) {
    std::vector<double> coords = generate_uniform(1000000);
    while (state.KeepRunning()) {
        delaunator::Delaunator delaunator(coords);
    }
}
} // namespace

BENCHMARK(BM_45K_geojson_nodes)->Unit(benchmark::kMillisecond);
BENCHMARK(BM_2K_uniform)->Unit(benchmark::kMillisecond);
BENCHMARK(BM_100K_uniform)->Unit(benchmark::kMillisecond);
BENCHMARK(BM_200K_uniform)->Unit(benchmark::kMillisecond);
BENCHMARK(BM_500K_uniform)->Unit(benchmark::kMillisecond);
BENCHMARK(BM_1M_uniform)->Unit(benchmark::kMillisecond);

BENCHMARK_MAIN()
