#pragma once

#include <cstdint>
#include <vector>
#include <string>
#include <glm/glm.hpp>

// Quantized-mesh-1.0 terrain format structures and utilities

struct QuantizedMeshHeader {
    // The center of the tile in Earth-centered Fixed coordinates
    double CenterX;
    double CenterY;
    double CenterZ;

    // The minimum and maximum heights in the area covered by this tile
    float MinimumHeight;
    float MaximumHeight;

    // The tile's bounding sphere
    double BoundingSphereCenterX;
    double BoundingSphereCenterY;
    double BoundingSphereCenterZ;
    double BoundingSphereRadius;

    // The horizon occlusion point
    double HorizonOcclusionPointX;
    double HorizonOcclusionPointY;
    double HorizonOcclusionPointZ;
};

struct TileBounds {
    double west;   // longitude in degrees
    double south;  // latitude in degrees
    double east;   // longitude in degrees
    double north;  // latitude in degrees
};

class QuantizedMeshTile {
public:
    QuantizedMeshTile();

    // Set vertices in geodetic coordinates (lon, lat, height)
    void SetVertices(const std::vector<glm::dvec3> &vertices);

    // Set triangle indices (counter-clockwise winding)
    void SetIndices(const std::vector<uint32_t> &indices);

    // Set edge indices (for skirt generation)
    void SetEdgeIndices(
        const std::vector<uint32_t> &west,
        const std::vector<uint32_t> &south,
        const std::vector<uint32_t> &east,
        const std::vector<uint32_t> &north);

    // Set the tile bounds for coordinate conversion
    void SetBounds(const TileBounds &bounds);

    // Write the quantized mesh to a file (gzipped)
    bool WriteToFile(const std::string &path) const;

private:
    std::vector<glm::dvec3> m_Vertices;
    std::vector<uint32_t> m_Indices;
    std::vector<uint32_t> m_WestIndices;
    std::vector<uint32_t> m_SouthIndices;
    std::vector<uint32_t> m_EastIndices;
    std::vector<uint32_t> m_NorthIndices;
    TileBounds m_Bounds;

    QuantizedMeshHeader ComputeHeader() const;
    void EncodeVertices(std::vector<uint8_t> &buffer) const;
    void EncodeIndices(std::vector<uint8_t> &buffer) const;
    void EncodeEdgeIndices(std::vector<uint8_t> &buffer) const;
};

// Tile coordinate utilities
TileBounds GetTileBounds(int zoom, int x, int y);

// Convert geodetic (lon, lat, height) to ECEF (Earth-Centered Earth-Fixed)
glm::dvec3 GeodeticToECEF(double lon, double lat, double height);

// ZigZag encoding for signed integers
uint32_t ZigZagEncode(int32_t value);

// High water mark encoding for indices
std::vector<uint32_t> HighWaterMarkEncode(const std::vector<uint32_t> &indices);
