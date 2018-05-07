#pragma once

struct BoundaryIndexedTriangle {
	int id;
	Point2D<double> vertices[3];
	int atlasVertexParentEdge[3];
	int atlasVertexIndices[3];
	int atlasEdgeIndices[3];
	HexagonalIndex indices;
};

struct AtlasIndexedPolygon {
	std::vector < Point2D < double >> vertices;
	std::vector < int > indices;
	std::vector < int > atlasVertexIndices;
	std::vector < int > atlasVertexParentEdge;
	std::vector < int > atlasEdgeIndices;
};

struct AtlasIndexedTriangle {
	int id;
	Point2D<double> vertices[3];
	int indices[3];
	int atlasVertexParentEdge[3];
	int atlasVertexIndices[3];
	int atlasEdgeIndices[3];
};

struct IndexedIntersectionPolygon {
	std::vector < Point2D < double >> vertices;
	std::vector < unsigned long long > indices;
	std::vector < int > edgeIndices;
};

struct IndexedIntersectionTriangle {
	Point2D < double > vertices[3];
	unsigned long long indices[3];
	int edgeIndices[3];
};


unsigned long long SetIntersectionKey(const unsigned long i0, const unsigned long i1) {
	return ( ( (static_cast<unsigned long long>(i0) << 32) & 0xFFFFFFFF00000000) | (static_cast<unsigned long long>(i1) & 0x00000000FFFFFFFF));
}

void GetIntersectionKey(unsigned long long key, unsigned long & i0, unsigned long & i1) {
	i1 = static_cast<unsigned long>(key & 0x00000000FFFFFFFF);
	i0 = static_cast<unsigned long>((key >> 32) & 0x00000000FFFFFFFF);
}

struct IntersectionInfo {
	unsigned long long intersectionKey;
	int intersectionIndex;
	Point2D<double> position;
	double time;
};


bool IntersectionComparison(const IntersectionInfo & i0, const IntersectionInfo & i1) {
	return i0.time < i1.time;
};

struct BoundarySegmentInfo {
	double startTime;
	double endTime;
	int halfEdge;
};

