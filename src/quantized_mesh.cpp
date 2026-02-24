#include "quantized_mesh.h"
#include <cmath>
#include <fstream>
#include <algorithm>
#include <zlib.h>
#include <iostream>
#include <cstring>

// WGS84 ellipsoid constants
static const double EARTH_RADIUS_A = 6378137.0;  // semi-major axis in meters
static const double EARTH_RADIUS_B = 6356752.314245;  // semi-minor axis in meters
static const double EARTH_FLATTENING = (EARTH_RADIUS_A - EARTH_RADIUS_B) / EARTH_RADIUS_A;
static const double EARTH_ECCENTRICITY_SQ = 2.0 * EARTH_FLATTENING - EARTH_FLATTENING * EARTH_FLATTENING;

static const double PI = 3.14159265358979323846;
static const double DEG_TO_RAD = PI / 180.0;

// Convert degrees to radians
static inline double toRadians(double degrees) {
    return degrees * DEG_TO_RAD;
}

glm::dvec3 GeodeticToECEF(double lon, double lat, double height) {
    const double lonRad = toRadians(lon);
    const double latRad = toRadians(lat);

    const double cosLat = std::cos(latRad);
    const double sinLat = std::sin(latRad);
    const double cosLon = std::cos(lonRad);
    const double sinLon = std::sin(lonRad);

    const double N = EARTH_RADIUS_A / std::sqrt(1.0 - EARTH_ECCENTRICITY_SQ * sinLat * sinLat);

    const double x = (N + height) * cosLat * cosLon;
    const double y = (N + height) * cosLat * sinLon;
    const double z = (N * (1.0 - EARTH_ECCENTRICITY_SQ) + height) * sinLat;

    return glm::dvec3(x, y, z);
}

TileBounds GetTileBounds(int zoom, int x, int y) {
    // EPSG:4326 (lat/lon) with TMS tiling scheme
    // At zoom level 0, there are 2 tiles covering the world horizontally
    // Each tile is 180 degrees wide and 180 degrees tall

    // For zoom level 0, we have 2 tiles horizontally (not 1)
    const int tilesX = 2 << zoom;  // 2 * 2^zoom for horizontal tiles
    const int tilesY = 1 << zoom;  // 2^zoom for vertical tiles

    const double tileWidth = 360.0 / tilesX;   // 180 degrees at zoom 0
    const double tileHeight = 180.0 / tilesY;  // 180 degrees at zoom 0

    TileBounds bounds;
    bounds.west = -180.0 + x * tileWidth;
    bounds.east = bounds.west + tileWidth;
    bounds.south = -90.0 + y * tileHeight;
    bounds.north = bounds.south + tileHeight;

    return bounds;
}

uint32_t ZigZagEncode(int32_t value) {
    return (value << 1) ^ (value >> 31);
}

std::vector<uint32_t> HighWaterMarkEncode(const std::vector<uint32_t> &indices) {
    std::vector<uint32_t> encoded(indices.size());
    uint32_t highest = 0;

    for (size_t i = 0; i < indices.size(); ++i) {
        const uint32_t idx = indices[i];
        encoded[i] = highest - idx;
        if (idx == highest) {
            ++highest;
        }
    }

    return encoded;
}

QuantizedMeshTile::QuantizedMeshTile() {
}

void QuantizedMeshTile::SetVertices(const std::vector<glm::dvec3> &vertices) {
    m_Vertices = vertices;
}

void QuantizedMeshTile::SetIndices(const std::vector<uint32_t> &indices) {
    m_Indices = indices;
}

void QuantizedMeshTile::SetEdgeIndices(
    const std::vector<uint32_t> &west,
    const std::vector<uint32_t> &south,
    const std::vector<uint32_t> &east,
    const std::vector<uint32_t> &north)
{
    m_WestIndices = west;
    m_SouthIndices = south;
    m_EastIndices = east;
    m_NorthIndices = north;
}

void QuantizedMeshTile::SetBounds(const TileBounds &bounds) {
    m_Bounds = bounds;
}

QuantizedMeshHeader QuantizedMeshTile::ComputeHeader() const {
    QuantizedMeshHeader header;

    // Find min/max heights
    float minHeight = m_Vertices[0].z;
    float maxHeight = m_Vertices[0].z;
    for (const auto &v : m_Vertices) {
        minHeight = std::min(minHeight, static_cast<float>(v.z));
        maxHeight = std::max(maxHeight, static_cast<float>(v.z));
    }

    header.MinimumHeight = minHeight;
    header.MaximumHeight = maxHeight;

    // Compute center of tile in ECEF
    const double centerLon = (m_Bounds.west + m_Bounds.east) / 2.0;
    const double centerLat = (m_Bounds.south + m_Bounds.north) / 2.0;
    const double centerHeight = (minHeight + maxHeight) / 2.0;
    const glm::dvec3 centerECEF = GeodeticToECEF(centerLon, centerLat, centerHeight);

    header.CenterX = centerECEF.x;
    header.CenterY = centerECEF.y;
    header.CenterZ = centerECEF.z;

    // Compute bounding sphere
    // Use center as sphere center and find max distance to any vertex
    double maxDistSq = 0.0;
    for (const auto &v : m_Vertices) {
        const glm::dvec3 ecef = GeodeticToECEF(v.x, v.y, v.z);
        const glm::dvec3 diff = ecef - centerECEF;
        const double distSq = glm::dot(diff, diff);
        maxDistSq = std::max(maxDistSq, distSq);
    }

    header.BoundingSphereCenterX = centerECEF.x;
    header.BoundingSphereCenterY = centerECEF.y;
    header.BoundingSphereCenterZ = centerECEF.z;
    header.BoundingSphereRadius = std::sqrt(maxDistSq);

    // Compute horizon occlusion point
    // This is a simplified computation - proper implementation would be more complex
    const glm::dvec3 scaledCenter = centerECEF / glm::dvec3(EARTH_RADIUS_A, EARTH_RADIUS_A, EARTH_RADIUS_B);
    header.HorizonOcclusionPointX = scaledCenter.x;
    header.HorizonOcclusionPointY = scaledCenter.y;
    header.HorizonOcclusionPointZ = scaledCenter.z;

    return header;
}

void QuantizedMeshTile::EncodeVertices(std::vector<uint8_t> &buffer) const {
    const uint32_t vertexCount = m_Vertices.size();

    // Find min/max heights for quantization
    float minHeight = m_Vertices[0].z;
    float maxHeight = m_Vertices[0].z;
    for (const auto &v : m_Vertices) {
        minHeight = std::min(minHeight, static_cast<float>(v.z));
        maxHeight = std::max(maxHeight, static_cast<float>(v.z));
    }

    // Quantize vertices to u, v, height
    std::vector<uint16_t> u(vertexCount);
    std::vector<uint16_t> v(vertexCount);
    std::vector<uint16_t> height(vertexCount);

    const double lonRange = m_Bounds.east - m_Bounds.west;
    const double latRange = m_Bounds.north - m_Bounds.south;
    const float heightRange = maxHeight - minHeight;

    for (uint32_t i = 0; i < vertexCount; ++i) {
        const double normalizedLon = (m_Vertices[i].x - m_Bounds.west) / lonRange;
        const double normalizedLat = (m_Vertices[i].y - m_Bounds.south) / latRange;
        const float normalizedHeight = heightRange > 0.0f ?
            (m_Vertices[i].z - minHeight) / heightRange : 0.0f;

        u[i] = static_cast<uint16_t>(normalizedLon * 32767.0);
        v[i] = static_cast<uint16_t>(normalizedLat * 32767.0);
        height[i] = static_cast<uint16_t>(normalizedHeight * 32767.0);
    }

    // Delta encode and zigzag encode
    std::vector<uint16_t> uEncoded(vertexCount);
    std::vector<uint16_t> vEncoded(vertexCount);
    std::vector<uint16_t> heightEncoded(vertexCount);

    uEncoded[0] = ZigZagEncode(u[0]);
    vEncoded[0] = ZigZagEncode(v[0]);
    heightEncoded[0] = ZigZagEncode(height[0]);

    for (uint32_t i = 1; i < vertexCount; ++i) {
        uEncoded[i] = ZigZagEncode(static_cast<int16_t>(u[i] - u[i-1]));
        vEncoded[i] = ZigZagEncode(static_cast<int16_t>(v[i] - v[i-1]));
        heightEncoded[i] = ZigZagEncode(static_cast<int16_t>(height[i] - height[i-1]));
    }

    // Write to buffer
    const size_t startSize = buffer.size();
    buffer.resize(startSize + sizeof(uint32_t) + vertexCount * 3 * sizeof(uint16_t));
    uint8_t *ptr = buffer.data() + startSize;

    // Write vertex count
    std::memcpy(ptr, &vertexCount, sizeof(uint32_t));
    ptr += sizeof(uint32_t);

    // Write u array
    std::memcpy(ptr, uEncoded.data(), vertexCount * sizeof(uint16_t));
    ptr += vertexCount * sizeof(uint16_t);

    // Write v array
    std::memcpy(ptr, vEncoded.data(), vertexCount * sizeof(uint16_t));
    ptr += vertexCount * sizeof(uint16_t);

    // Write height array
    std::memcpy(ptr, heightEncoded.data(), vertexCount * sizeof(uint16_t));
}

void QuantizedMeshTile::EncodeIndices(std::vector<uint8_t> &buffer) const {
    const uint32_t triangleCount = m_Indices.size() / 3;
    const bool use32Bit = m_Vertices.size() > 65536;

    // Add padding for alignment
    if (use32Bit) {
        while (buffer.size() % 4 != 0) {
            buffer.push_back(0);
        }
    } else {
        while (buffer.size() % 2 != 0) {
            buffer.push_back(0);
        }
    }

    // High water mark encode indices
    const std::vector<uint32_t> encoded = HighWaterMarkEncode(m_Indices);

    const size_t startSize = buffer.size();

    if (use32Bit) {
        buffer.resize(startSize + sizeof(uint32_t) + encoded.size() * sizeof(uint32_t));
        uint8_t *ptr = buffer.data() + startSize;

        std::memcpy(ptr, &triangleCount, sizeof(uint32_t));
        ptr += sizeof(uint32_t);

        std::memcpy(ptr, encoded.data(), encoded.size() * sizeof(uint32_t));
    } else {
        buffer.resize(startSize + sizeof(uint32_t) + encoded.size() * sizeof(uint16_t));
        uint8_t *ptr = buffer.data() + startSize;

        std::memcpy(ptr, &triangleCount, sizeof(uint32_t));
        ptr += sizeof(uint32_t);

        std::vector<uint16_t> encoded16(encoded.size());
        for (size_t i = 0; i < encoded.size(); ++i) {
            encoded16[i] = static_cast<uint16_t>(encoded[i]);
        }
        std::memcpy(ptr, encoded16.data(), encoded16.size() * sizeof(uint16_t));
    }
}

void QuantizedMeshTile::EncodeEdgeIndices(std::vector<uint8_t> &buffer) const {
    const bool use32Bit = m_Vertices.size() > 65536;

    auto encodeEdge = [&](const std::vector<uint32_t> &indices) {
        const uint32_t count = indices.size();
        const size_t startSize = buffer.size();

        if (use32Bit) {
            buffer.resize(startSize + sizeof(uint32_t) + count * sizeof(uint32_t));
            uint8_t *ptr = buffer.data() + startSize;

            std::memcpy(ptr, &count, sizeof(uint32_t));
            ptr += sizeof(uint32_t);

            std::memcpy(ptr, indices.data(), count * sizeof(uint32_t));
        } else {
            buffer.resize(startSize + sizeof(uint32_t) + count * sizeof(uint16_t));
            uint8_t *ptr = buffer.data() + startSize;

            std::memcpy(ptr, &count, sizeof(uint32_t));
            ptr += sizeof(uint32_t);

            std::vector<uint16_t> indices16(count);
            for (uint32_t i = 0; i < count; ++i) {
                indices16[i] = static_cast<uint16_t>(indices[i]);
            }
            std::memcpy(ptr, indices16.data(), count * sizeof(uint16_t));
        }
    };

    encodeEdge(m_WestIndices);
    encodeEdge(m_SouthIndices);
    encodeEdge(m_EastIndices);
    encodeEdge(m_NorthIndices);
}

bool QuantizedMeshTile::WriteToFile(const std::string &path) const {
    if (m_Vertices.empty() || m_Indices.empty()) {
        std::cerr << "Cannot write empty tile" << std::endl;
        return false;
    }

    // Build uncompressed data
    std::vector<uint8_t> data;

    // Write header
    const QuantizedMeshHeader header = ComputeHeader();
    const size_t headerSize = sizeof(QuantizedMeshHeader);
    data.resize(headerSize);
    std::memcpy(data.data(), &header, headerSize);

    // Write vertex data
    EncodeVertices(data);

    // Write index data
    EncodeIndices(data);

    // Write edge indices
    EncodeEdgeIndices(data);

    // Compress with gzip
    uLongf compressedSize = compressBound(data.size());
    std::vector<uint8_t> compressed(compressedSize);

    const int result = compress(compressed.data(), &compressedSize, data.data(), data.size());
    if (result != Z_OK) {
        std::cerr << "Failed to compress tile data" << std::endl;
        return false;
    }

    compressed.resize(compressedSize);

    // Write to file
    std::ofstream file(path, std::ios::binary);
    if (!file) {
        std::cerr << "Failed to open file for writing: " << path << std::endl;
        return false;
    }

    file.write(reinterpret_cast<const char*>(compressed.data()), compressed.size());
    file.close();

    return true;
}
