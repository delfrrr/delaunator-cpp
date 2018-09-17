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

const start = Date.now();
const delaunator = new Delaunator(coords);
const end = Date.now();

console.log('points =', coords.length / 2);
console.log('milliseconds =', end - start);
console.log('triangles =', delaunator.triangles.length);

if (outputFile) {
    console.log('Saving to', outputFile)
    const trianglesAr = Array.from(delaunator.triangles);
    writeFileSync(outputFile, JSON.stringify(trianglesAr));
}

