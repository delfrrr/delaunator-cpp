#!/usr/bin/env bash

mkdir -p ./includes
rm -rf ./includes/*

git clone --branch v1.1.0 --depth=1 https://github.com/Tencent/rapidjson.git ./includes/rapidjson
git clone --branch master --depth=1 https://github.com/louisdx/cxx-prettyprint.git ./includes/prettyprint
git clone --branch master --depth=1 https://github.com/catchorg/Catch2.git ./includes/catch
