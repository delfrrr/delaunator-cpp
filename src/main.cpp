// #include <iostream>
#include "rapidjson/document.h"
#include <cstdio>
#include <fstream>
#include <string>
#include <exception>

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

int main(int, char* argv[]) {
    const char* filename = argv[1];
    std::string json = read_file(filename);
    // std::printf("%s\n", json.c_str());
    // std::printf("%s", json.c_str());
}
