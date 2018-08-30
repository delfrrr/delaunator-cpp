#!/usr/bin/env node

/**
 * generates triangles indices from delaunator;
 */

const {readFileSync, writeFileSync} = require('fs');
const Delaunator = require('delaunator');
const inputFile = process.argv[2];
const outputFile = process.argv[3];
const geoJson = JSON.parse(readFileSync(inputFile));

const n = geoJson.features.length;
const coords = new Float64Array(n * 2);
for(let i = 0; i < n; i++) {
    const f = geoJson.features[i];
    coords[2 * i] = f.geometry.coordinates[0];
    coords[2 * i + 1] = f.geometry.coordinates[1];
}

console.time('Delaunator');
const delaunator = new Delaunator(coords);
console.timeEnd('Delaunator');

const trianglesAr = Array.from(delaunator.triangles);
writeFileSync(outputFile, JSON.stringify(trianglesAr));

