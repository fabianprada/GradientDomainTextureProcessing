#ifndef ATLAS_MESH_INLCUDED
#define ATLAS_MESH_INLCUDED

#include "SimpleMesh.h"
#include "ChartDecomposition.h"
#include <set>

class AtlasMesh{
public:
	std::vector<Point2D<double>> vertices;
	std::vector<TriangleIndex> triangles;
	std::vector<int> triangleIndexInChart;
	std::vector<int> triangleChartIndex;
	std::vector<int> halfEdgeToEdgeIndex;
	std::vector<int> vertexMap;
	//std::vector<SquareMatrix<double, 2>> metric;
	int numCharts;
};

class AtlasChart {
public:
	Point2D<double> minCorner;
	Point2D<double> maxCorner;
	Point2D<double> gridOrigin;
	int originCoords[2];
	std::vector<TriangleIndex> triangles;
	std::vector<Point2D<double>> vertices;
	std::vector<int> boundaryHalfEdges;
	std::vector<int> atlasEdgeIndices;

	std::vector<int> meshVertexIndices;
	std::vector<int> meshTriangleIndices;
	//std::vector<SquareMatrix<double, 2>> metric;
};

class IndexedVector2D{
public:
	IndexedVector2D(Point2D<double> p_p, int p_index, int p_vertex){
		p = p_p;
		index = p_index;
		vertex = p_vertex;
	}
	Point2D<double> p;
	int index;
	int vertex;
};

class IndexedVector2DComparison
{
public:
	bool operator()(const IndexedVector2D & p1, const IndexedVector2D & p2) const{
		for (int i = 0; i < 2; i++){
			if (p1.p[i] < p2.p[i])
				return true;
			else if (p2.p[i] < p1.p[i])
				return false;
			else {
				if (p1.vertex < p2.vertex) return true;
				else if (p2.vertex < p1.vertex) return false;
			}
		}
		return false;
	}
};

#include "AtlasMesh.inl"
#include "AtlasCharts.inl"
#endif// ATLAS_MESH_INLCUDED