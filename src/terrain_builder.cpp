#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <sys/stat.h>
#include <sys/types.h>

#include <gdal_priv.h>
#include <cpl_conv.h>

#include "quantized_mesh.h"
#include "triangulator.h"
#include "heightmap.h"
#include "cmdline.h"

// Helper to create directories recursively
bool CreateDirectories(const std::string &path) {
    size_t pos = 0;
    while ((pos = path.find('/', pos + 1)) != std::string::npos) {
        std::string dir = path.substr(0, pos);
        mkdir(dir.c_str(), 0755);
    }
    mkdir(path.c_str(), 0755);
    return true;
}

// Extract a tile from a GDAL dataset
std::shared_ptr<Heightmap> ExtractTile(
    GDALDataset *dataset,
    int zoom, int x, int y)
{
    GDALRasterBand *band = dataset->GetRasterBand(1);
    if (!band) {
        std::cerr << "Failed to get raster band" << std::endl;
        return nullptr;
    }

    // Get dataset geotransform
    double geot[6];
    if (dataset->GetGeoTransform(geot) != CE_None) {
        std::cerr << "Failed to get geotransform" << std::endl;
        return nullptr;
    }

    // geot[0] = top left x
    // geot[1] = w-e pixel resolution
    // geot[2] = rotation, 0 if north up
    // geot[3] = top left y
    // geot[4] = rotation, 0 if north up
    // geot[5] = n-s pixel resolution (negative value)

    const TileBounds bounds = GetTileBounds(zoom, x, y);

    // Convert geographic bounds to pixel coordinates
    const int datasetWidth = band->GetXSize();
    const int datasetHeight = band->GetYSize();

    // Calculate pixel coordinates
    const double pixelX1 = (bounds.west - geot[0]) / geot[1];
    const double pixelY1 = (bounds.north - geot[3]) / geot[5];
    const double pixelX2 = (bounds.east - geot[0]) / geot[1];
    const double pixelY2 = (bounds.south - geot[3]) / geot[5];

    int xOff = static_cast<int>(std::floor(pixelX1));
    int yOff = static_cast<int>(std::floor(pixelY1));
    int xSize = static_cast<int>(std::ceil(pixelX2)) - xOff;
    int ySize = static_cast<int>(std::ceil(pixelY2)) - yOff;

    // Clamp to dataset bounds
    if (xOff < 0) { xSize += xOff; xOff = 0; }
    if (yOff < 0) { ySize += yOff; yOff = 0; }
    if (xOff + xSize > datasetWidth) { xSize = datasetWidth - xOff; }
    if (yOff + ySize > datasetHeight) { ySize = datasetHeight - yOff; }

    if (xSize <= 0 || ySize <= 0) {
        std::cerr << "Tile outside dataset bounds" << std::endl;
        return nullptr;
    }

    // Read tile data
    const int bufWidth = std::min(xSize, 256);
    const int bufHeight = std::min(ySize, 256);

    std::vector<float> data(bufWidth * bufHeight);

    CPLErr err = band->RasterIO(
        GF_Read,
        xOff, yOff, xSize, ySize,
        data.data(), bufWidth, bufHeight,
        GDT_Float32, 0, 0);

    if (err != CE_None) {
        std::cerr << "Failed to read raster data" << std::endl;
        return nullptr;
    }

    return std::make_shared<Heightmap>(bufWidth, bufHeight, data);
}

// Build a single terrain tile
bool BuildTile(
    GDALDataset *dataset,
    const std::string &outputDir,
    int zoom, int x, int y,
    float maxError)
{
    // Extract heightmap for this tile
    auto hm = ExtractTile(dataset, zoom, x, y);
    if (!hm || hm->Width() * hm->Height() == 0) {
        std::cerr << "Failed to extract tile " << zoom << "/" << x << "/" << y << std::endl;
        return false;
    }

    // Triangulate the heightmap
    Triangulator tri(hm);
    tri.Run(maxError, 0, 0);

    // Get triangulation results
    const std::vector<glm::ivec2> points2D = tri.Points2D();
    const std::vector<int> triangles = tri.TriangleIndices();

    if (points2D.empty() || triangles.empty()) {
        std::cerr << "Triangulation produced no geometry" << std::endl;
        return false;
    }

    // Convert to geodetic coordinates
    const TileBounds bounds = GetTileBounds(zoom, x, y);
    const int width = hm->Width();
    const int height = hm->Height();

    std::vector<glm::dvec3> vertices;
    vertices.reserve(points2D.size());

    for (const auto &p : points2D) {
        const double u = static_cast<double>(p.x) / (width - 1);
        const double v = static_cast<double>(p.y) / (height - 1);

        const double lon = bounds.west + u * (bounds.east - bounds.west);
        const double lat = bounds.south + v * (bounds.north - bounds.south);
        const float h = hm->At(p.x, p.y);

        vertices.push_back(glm::dvec3(lon, lat, h));
    }

    // Convert triangle indices to uint32_t
    std::vector<uint32_t> indices;
    indices.reserve(triangles.size());
    for (int idx : triangles) {
        indices.push_back(static_cast<uint32_t>(idx));
    }

    // Find edge vertices
    std::vector<uint32_t> westIndices, southIndices, eastIndices, northIndices;

    for (size_t i = 0; i < points2D.size(); ++i) {
        const glm::ivec2 &p = points2D[i];

        if (p.x == 0) westIndices.push_back(i);
        if (p.y == 0) southIndices.push_back(i);
        if (p.x == width - 1) eastIndices.push_back(i);
        if (p.y == height - 1) northIndices.push_back(i);
    }

    // Create quantized mesh tile
    QuantizedMeshTile tile;
    tile.SetBounds(bounds);
    tile.SetVertices(vertices);
    tile.SetIndices(indices);
    tile.SetEdgeIndices(westIndices, southIndices, eastIndices, northIndices);

    // Write to file
    const std::string tileDir = outputDir + "/" + std::to_string(zoom) + "/" + std::to_string(x);
    CreateDirectories(tileDir);

    const std::string tilePath = tileDir + "/" + std::to_string(y) + ".terrain";
    return tile.WriteToFile(tilePath);
}

int main(int argc, char **argv) {
    // Parse command line arguments
    cmdline::parser p;

    p.add<std::string>("input", 'i', "input GDAL data source (VRT or GeoTIFF)", true);
    p.add<std::string>("output", 'o', "output directory for terrain tiles", true);
    p.add<int>("min-zoom", '\0', "minimum zoom level", false, 0);
    p.add<int>("max-zoom", '\0', "maximum zoom level", false, 10);
    p.add<float>("error", 'e', "maximum triangulation error", false, 0.001);

    p.parse_check(argc, argv);

    const std::string inputPath = p.get<std::string>("input");
    const std::string outputDir = p.get<std::string>("output");
    const int minZoom = p.get<int>("min-zoom");
    const int maxZoom = p.get<int>("max-zoom");
    const float maxError = p.get<float>("error");

    // Initialize GDAL
    GDALAllRegister();

    // Open input dataset
    GDALDataset *dataset = static_cast<GDALDataset*>(GDALOpen(inputPath.c_str(), GA_ReadOnly));
    if (!dataset) {
        std::cerr << "Failed to open input file: " << inputPath << std::endl;
        return 1;
    }

    std::cout << "Dataset: " << dataset->GetRasterXSize() << " x " << dataset->GetRasterYSize() << std::endl;

    // Create output directory
    CreateDirectories(outputDir);

    // Generate tiles for each zoom level
    for (int zoom = minZoom; zoom <= maxZoom; ++zoom) {
        const int tilesPerSide = 1 << zoom;  // 2^zoom

        std::cout << "Generating zoom level " << zoom << " (" << tilesPerSide << " x " << tilesPerSide << " tiles)" << std::endl;

        int successCount = 0;
        int totalTiles = tilesPerSide * tilesPerSide;

        for (int y = 0; y < tilesPerSide; ++y) {
            for (int x = 0; x < tilesPerSide; ++x) {
                if (BuildTile(dataset, outputDir, zoom, x, y, maxError)) {
                    ++successCount;
                }

                // Progress indicator
                const int processed = y * tilesPerSide + x + 1;
                if (processed % 10 == 0 || processed == totalTiles) {
                    std::cout << "  Progress: " << processed << "/" << totalTiles << "\r" << std::flush;
                }
            }
        }

        std::cout << std::endl;
        std::cout << "  Generated " << successCount << "/" << totalTiles << " tiles" << std::endl;
    }

    GDALClose(dataset);

    std::cout << "Done!" << std::endl;
    return 0;
}
