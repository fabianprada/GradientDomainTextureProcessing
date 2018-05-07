#ifndef HIERARCHICAL_SYSTEM_INCLUDED
#define HIERARCHICAL_SYSTEM_INCLUDED

#include "IndexedPolygon.h"

template<class Real>
struct CellLineSampleInfo {
	Point2D<Real> pos;
	SquareMatrix<Real, 2> tensor;
	Point2D<Real> v[4];
	int cellOffset;
};

template<class Real>
bool CellLineSampleCompare(const CellLineSampleInfo<Real> &a, const CellLineSampleInfo<Real> &b) {
	return a.cellOffset < b.cellOffset;
}


template<class Real>
struct HexagonalSampleInfo {
	int cellId;
	Point2D<Real> pos;
	SquareMatrix<Real, 2> tensor;
	Point2D<Real> v[6];
	int fineNodes[6];
};

class CellIntersection {
public:
	int chartId;
	int cellId;
	double time;
};

bool CellIntersectionComparsion(const CellIntersection & a, const CellIntersection & b) {
	return a.time < b.time;
}



class RasterLine {
public:
	int lineStartIndex;
	int lineEndIndex;
	int prevLineIndex;
	int nextLineIndex;
	int coeffStartIndex;
};

class DeepLine {
public:
	int coarseLineStartIndex;
	int coarseLineEndIndex;
	int finePrevLineIndex;
	int fineCurrentLineIndex;
	int fineNextLineIndex;
};

class SegmentedRasterLine {
public:
	std::vector<RasterLine> segments;
};

class MultigridBlockInfo {
public:
	MultigridBlockInfo(const int p_blockWidth = 128, const int p_blockHeight = 16, const int p_paddingWidth = 2, const int p_paddingHeight = 2, const int p_boundaryDilation = 0) {
		blockWidth = p_blockWidth;
		blockHeight = p_blockHeight;
		paddingWidth = p_paddingWidth;
		paddingHeight = p_paddingHeight;
		boundaryDilation = p_boundaryDilation;
	}
	int blockWidth;
	int blockHeight;
	int paddingWidth;
	int paddingHeight;
	int boundaryDilation;
};

class BlockDeepSegment {
public:
	int currentStart;
	int currentEnd;
	int previousStart;
	int nextStart;
	int deepStart;
};

class BlockDeepSegmentedLine {
public:
	std::vector<BlockDeepSegment> blockDeepSegments;
};

class BlockTask {
public:
	std::vector<BlockDeepSegmentedLine> blockPaddedSegmentedLines;
	std::vector<BlockDeepSegmentedLine> blockDeepSegmentedLines;
};
class ThreadTask {
public:
	int taskDeepTexels;
	std::vector<BlockTask> blockTasks;
};

bool threadTaskComparison(const ThreadTask & task1, const ThreadTask & task2) {
	return task1.taskDeepTexels < task2.taskDeepTexels;
}

class InteriorTexelToCellLine{
public:
	int texelStartIndex;
	int texelEndIndex;
	int coeffOffset;
	int length;
	int previousCellStartIndex;
	int nextCellStartIndex;
};


class ProlongationLine{
public:
	int startIndex;
	int length;
	int centerLineIndex;
	int prevLineIndex;
	int nextLineIndex;
	bool alignedStart;
};

class RestrictionLine{
public:
	RestrictionLine() {
		startIndex = length = centerLineIndex = prevLineIndex = nextLineIndex = -1;
	}
	int startIndex;
	int length;
	int centerLineIndex;
	int prevLineIndex;
	int nextLineIndex;
	int outputStart;
};

class AuxiliarNode {
public:
	Point2D<double> position;
	int index;
};



class GridChart {
public:
	Point2D<double> corner;
	int cornerCoords[2];
	int centerOffset[2];
	double cellSizeW;
	double cellSizeH;
	int width;
	int height;
	int globalTexelDeepOffset;
	int globalTexelBoundaryOffset;
	int globalIndexTexelOffset;
	int globalIndexInteriorTexelOffset;
	int globalIndexCellOffset;
	int globalIndexInteriorCellOffset;
	int gridIndexOffset;
	Image<int> cellType;
	Image<int> nodeType;
	Image<int> globalInteriorTexelIndex;
	Image<int> globalTexelIndex;
	Image<int> globalTexelDeepIndex;
	Image<int> globalTexelBoundaryIndex;


	Image<int> localCellIndex;
	std::vector<CellIndex> cellIndices;

	Image<int> localInteriorCellIndex;
	std::vector<CellIndex> interiorCellCorners;

	std::vector<CellIndex> interiorCellGlobalCorners;

	int numInteriorCells;

	Image<int> localBoundaryCellIndex;
	int numBoundaryCells;

	std::vector<int> interiorCellIndexToLocalCellIndex;
	std::vector<int> boundaryCellIndexToLocalCellIndex;

	std::vector<std::vector<AtlasIndexedPolygon>> boundaryPolygons;
	std::vector<std::vector<BoundaryIndexedTriangle>> boundaryTriangles;

	int numBoundaryTriangles;

	std::vector<AuxiliarNode> auxiliarNodes;

	Image<int> triangleId;
	Image<Point2D<double>> baricentricCoords;
};


class GridNodeInfo {
public:
	int chartId;
	int ci;
	int cj;
	int nodeType;
};

class TexelLineInfo {
public:
	TexelLineInfo() {
		lineIndex = -1;
		offset = -1;
	}
	int lineIndex;
	int offset;
};

class GridAtlas {
public:

	std::vector<ThreadTask> threadTasks;
	std::vector<int> boundaryAndDeepIndex;
	std::vector<int> boundaryGlobalIndex;
	std::vector<int> deepGlobalIndex;
	std::vector<GridNodeInfo> nodeInfo;
	std::vector<GridChart> gridCharts;
	std::vector<SegmentedRasterLine> segmentedLines;
	std::vector<RasterLine> rasterLines;
	std::vector<RasterLine> restrictionLines;
	std::vector<DeepLine> deepLines;
	std::vector<ProlongationLine> prolongationLines;
	Eigen::SparseMatrix<double> coarseToFineNodeProlongation;
	//int resolution;
	int numTexels;
	int numInteriorTexels;
	int numDeepTexels;
	int numBoundaryTexels;
	int numCells;
	int numBoundaryCells;
	int numInteriorCells;
	int numBoundaryNodes;
	int numMidPoints;
	int numFineNodes;
};

struct BoundaryDeepIndex {
	int boundaryIndex;
	int deepGlobalIndex;
	int deepIndex;
	int offset;
};

struct BoundaryBoundaryIndex {
	int coarsePrincipalBoundaryIndex;
	int coarseSecondaryBoundaryIndex;
	int fineDeepIndex;
	int offset;
	double weight;
};

class HierarchicalSystem {
public:
	std::vector<SparseMatrix<double, int>> boundaryRestriction;
	std::vector<SparseMatrix<double, int>> prolongation;

	std::vector<GridAtlas> gridAtlases;

	std::vector<SparseMatrix<double, int>> boundaryCoarseFineProlongation;
	std::vector<SparseMatrix<double, int>> boundaryFineCoarseRestriction;
	std::vector<std::vector<BoundaryDeepIndex>> boundaryDeepIndices;
	std::vector<std::vector<BoundaryBoundaryIndex>> boundaryBoundaryIndices;
};


template<class Real>
class MultigridLevelCoefficients{
public:
	std::vector<Real> deepCoefficients;
	SparseMatrix<Real, int> boundaryDeepMatrix;
	SparseMatrix<Real, int> boundaryBoundaryMatrix;
};

template<class DataType>
class MultigridLevelVariables {
public:
	std::vector<DataType> x;
	std::vector<DataType> rhs;
	std::vector<DataType> residual;
	std::vector<DataType> boundary_rhs;
	std::vector<DataType> boundary_value;
	std::vector<DataType> variable_boundary_value;
};

template<class Real>
class MultigridLevelIndices {
public:
	std::vector<ThreadTask> threadTasks;
	std::vector<int> boundaryGlobalIndex;
	std::vector<SegmentedRasterLine> segmentedLines;
	std::vector<RasterLine> rasterLines;
	std::vector<RasterLine> restrictionLines;
	std::vector<ProlongationLine> prolongationLines;
	SparseMatrix<Real, int> boundaryRestriction;
};


#include "SparseMatrixParser.h"
#include "Atlas.h"

#include "Metric.inl"
#include "BoundaryQuadraticElements.inl"
#include "ProlongationAndRestriction.inl"
#include "HierarchyConstruction.inl"
#include "UpdateCoefficients.inl"
#include "IterativeSolvers.inl"



#endif//HIERARCHICAL_SYSTEM_INCLUDED
