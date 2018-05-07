#ifndef POLYGON_CLIPPING_INCLUDED
#define POLYGON_CLIPPING_INCLUDED
#include <Misha/Geometry.h>
#include "IndexedPolygon.h"

void ClipConvexPolygon(std::vector<Point2D<double>> & vertices, const Point2D<double> & normal, const double & offset){

	int verticesCount = vertices.size();

	std::vector<Point2D<double>> outputVertices;

	Point2D<double> previousVertex = vertices[verticesCount - 1];
	double previousLevel = Point2D<double>::Dot(previousVertex, normal) - offset;
	bool isPreviousInterior =  previousLevel > 0;
	for (int i = 0; i < vertices.size(); i++){
		Point2D<double> currentVertex = vertices[i];
		double currentLevel = Point2D<double>::Dot(currentVertex, normal) - offset;
		bool isInterior = currentLevel > 0;
		if (isInterior){
			if (!isPreviousInterior){//Entrying edge
				//Add intersection
				double alpha = -currentLevel / (previousLevel - currentLevel);
				Point2D<double> intersection = previousVertex*alpha + currentVertex*(1.0 - alpha);
				outputVertices.push_back(intersection);
			}
			outputVertices.push_back(currentVertex);
		}
		else{
			if (isPreviousInterior){//Exiting edge
				//Add intersection
				double alpha = -currentLevel / (previousLevel - currentLevel);
				Point2D<double> intersection = previousVertex*alpha + currentVertex*(1.0 - alpha);
				outputVertices.push_back(intersection);
			}
		}
		previousVertex = currentVertex;
		previousLevel = currentLevel;
		isPreviousInterior = isInterior;
	}

	vertices = outputVertices;
}

int ClipTriangleToCell(std::vector<Point2D<double>> & vertices, int i, int j, int width, int height, bool  verbose = false){
	if (verbose) for (int v = 0; v < vertices.size(); v++) printf("(%f %f) ", vertices[v][0], vertices[v][1]);
	if (verbose) for (int v = 0; v < vertices.size(); v++) printf("(%f %f) ", vertices[v][0] * double(width) - (double(i) + 0.5), vertices[v][1] * double(height) - (double(j) + 0.5));

	Point2D<double> normal(0, 1);
	double offset = (double(j) + 0.5) / static_cast<double>(height);
	ClipConvexPolygon(vertices, normal, offset);

	if (verbose) printf("\n Clipping N (%f %f) O (%f) \n", normal[0], normal[1], offset * double(height) - (double(j) + 0.5));
	if (verbose) for (int v = 0; v < vertices.size(); v++) printf("(%f %f) ", vertices[v][0] * double(width) - (double(i) + 0.5), vertices[v][1] * double(height) - (double(j) + 0.5));

	if (!vertices.size()) return 0;

	normal = Point2D<double>(0, -1);
	offset = -(double(j + 1) + 0.5) / static_cast<double>(height);
	ClipConvexPolygon(vertices, normal, offset);

	if (verbose) printf("\n Clipping N (%f %f) O (%f) \n", normal[0], normal[1], (offset * double(height) + (double(j) + 0.5)));
	if (verbose) for (int v = 0; v < vertices.size(); v++) printf("(%f %f) ", vertices[v][0] * double(width) - (double(i) + 0.5), vertices[v][1] * double(height) - (double(j) + 0.5));

	if (!vertices.size()) return 0;

	normal = Point2D<double>(1, 0);
	offset = (double(i) + 0.5) / static_cast<double>(width);
	ClipConvexPolygon(vertices, normal, offset);

	if (verbose) printf("\n Clipping N (%f %f) O (%f) \n", normal[0], normal[1], offset* double(width) - (double(i) + 0.5));
	if (verbose) for (int v = 0; v < vertices.size(); v++) printf("(%f %f) ", vertices[v][0] * double(width) - (double(i) + 0.5), vertices[v][1] * double(height) - (double(j) + 0.5));

	if (!vertices.size()) return 0;

	normal = Point2D<double>(-1, 0);
	offset = -(double(i + 1) + 0.5) / static_cast<double>(width);
	ClipConvexPolygon(vertices, normal, offset);

	if (verbose) printf("\n Clipping N (%f %f) O (%f) \n", normal[0], normal[1], -(offset* double(width) + (double(i) + 0.5)));
	if (verbose) for (int v = 0; v < vertices.size(); v++) printf("(%f %f) ", vertices[v][0] * double(width) - (double(i) + 0.5), vertices[v][1] * double(height) - (double(j) + 0.5));

	return vertices.size();
}

int ClipTriangleToPrimalCell(std::vector<Point2D<double>> & vertices, int i, int j, double cellSizeW, double cellSizeH, bool  verbose = false){
	if (verbose) for (int v = 0; v < vertices.size(); v++) printf("(%f %f) ", vertices[v][0], vertices[v][1]);
	//if (verbose) for (int v = 0; v < vertices.size(); v++) printf("(%f %f) ", vertices[v][0] * double(width) - (double(i) + 0.5), vertices[v][1] * double(height) - (double(j) + 0.5));

	Point2D<double> normal(0, 1);
	double offset = double(j) * cellSizeH;
	ClipConvexPolygon(vertices, normal, offset);

	//if (verbose) printf("\n Clipping N (%f %f) O (%f) \n", normal[0], normal[1], offset * double(height) - (double(j) + 0.5));
	//if (verbose) for (int v = 0; v < vertices.size(); v++) printf("(%f %f) ", vertices[v][0] * double(width) - (double(i) + 0.5), vertices[v][1] * double(height) - (double(j) + 0.5));

	if (!vertices.size()) return 0;

	normal = Point2D<double>(0, -1);
	offset = -(double(j + 1)*cellSizeH);
	ClipConvexPolygon(vertices, normal, offset);

	//if (verbose) printf("\n Clipping N (%f %f) O (%f) \n", normal[0], normal[1], (offset * double(height) + (double(j) + 0.5)));
	//if (verbose) for (int v = 0; v < vertices.size(); v++) printf("(%f %f) ", vertices[v][0] * double(width) - (double(i) + 0.5), vertices[v][1] * double(height) - (double(j) + 0.5));

	if (!vertices.size()) return 0;

	normal = Point2D<double>(1, 0);
	offset = double(i) *cellSizeW;
	ClipConvexPolygon(vertices, normal, offset);

	//if (verbose) printf("\n Clipping N (%f %f) O (%f) \n", normal[0], normal[1], offset* double(width) - (double(i) + 0.5));
	//if (verbose) for (int v = 0; v < vertices.size(); v++) printf("(%f %f) ", vertices[v][0] * double(width) - (double(i) + 0.5), vertices[v][1] * double(height) - (double(j) + 0.5));

	if (!vertices.size()) return 0;

	normal = Point2D<double>(-1, 0);
	offset = -(double(i + 1)*cellSizeW);
	ClipConvexPolygon(vertices, normal, offset);

	//if (verbose) printf("\n Clipping N (%f %f) O (%f) \n", normal[0], normal[1], -(offset* double(width) + (double(i) + 0.5)));
    //if (verbose) for (int v = 0; v < vertices.size(); v++) printf("(%f %f) ", vertices[v][0] * double(width) - (double(i) + 0.5), vertices[v][1] * double(height) - (double(j) + 0.5));

	return vertices.size();
}


int ClipTriangleToTriangle(std::vector<Point2D<double>> & vertices, std::vector<Point2D<double>> & boundaryVertices, double colinearError = 1e-10, bool  verbose = false){
	
	std::vector<Point3D<double>> inputLines(vertices.size());
	Point2D<double> inputCenter = (vertices[0] + vertices[1] + vertices[2]) / 3.0;
	for (int v = 0; v < vertices.size(); v++){
		Point2D<double> dir = vertices[(v + 1) % 3] - vertices[v];
		Point2D<double> mid = (vertices[(v + 1) % 3] + vertices[v])/2.0;
		Point2D<double> normal(-dir[1], dir[0]);
		normal /= Point2D<double>::Length(normal);
		double offset = Point2D<double>::Dot(mid, normal);
		if (Point2D<double>::Dot(inputCenter, normal) < offset){
			normal *= -1.0;
			offset *= -1.0;
		}
		inputLines[v][0] = normal[0];
		inputLines[v][1] = normal[1];
		inputLines[v][2] = offset;
	}

	Point2D<double> boundaryCenter = (boundaryVertices[0] + boundaryVertices[1] + boundaryVertices[2]) / 3.0;
	for (int k = 0; k < boundaryVertices.size(); k++){
		Point2D<double> dir = boundaryVertices[(k + 1) % 3] - boundaryVertices[k];
		Point2D<double> mid = (boundaryVertices[(k + 1) % 3] + boundaryVertices[k]) / 2.0;
		Point2D<double> normal(-dir[1], dir[0]);
		normal /= Point2D<double>::Length(normal);
		double offset = Point2D<double>::Dot(mid, normal);

		if (Point2D<double>::Dot(boundaryCenter, normal) < offset){
			normal *= -1.0;
			offset *= -1.0;
		}
		offset -= 1e-12;
		Point3D<double> boundaryLine(normal[0], normal[1], offset);

		//Check if boundary line match inputLine
		bool colinearSegment = false;
		for (int v = 0; v < inputLines.size(); v++){
			if (fabs(inputLines[v][0] - boundaryLine[0]) < colinearError && fabs(inputLines[v][1] - boundaryLine[1]) < colinearError && fabs(inputLines[v][2] - boundaryLine[2]) < colinearError){
				//printf("Found colinearity! \n");
				colinearSegment = true;
			}
		}
		if (!colinearSegment){
			ClipConvexPolygon(vertices, normal, offset);
		}
	}

	return vertices.size();
}

int ClipTriangleToTriangle(std::vector<Point2D<double>> & vertices, Point2D<double> boundaryVertices[3], double colinearError = 1e-10, bool  verbose = false) {

	std::vector<Point3D<double>> inputLines(vertices.size());
	Point2D<double> inputCenter = (vertices[0] + vertices[1] + vertices[2]) / 3.0;
	for (int v = 0; v < vertices.size(); v++) {
		Point2D<double> dir = vertices[(v + 1) % 3] - vertices[v];
		Point2D<double> mid = (vertices[(v + 1) % 3] + vertices[v]) / 2.0;
		Point2D<double> normal(-dir[1], dir[0]);
		normal /= Point2D<double>::Length(normal);
		double offset = Point2D<double>::Dot(mid, normal);
		if (Point2D<double>::Dot(inputCenter, normal) < offset) {
			normal *= -1.0;
			offset *= -1.0;
		}
		inputLines[v][0] = normal[0];
		inputLines[v][1] = normal[1];
		inputLines[v][2] = offset;
	}

	Point2D<double> boundaryCenter = (boundaryVertices[0] + boundaryVertices[1] + boundaryVertices[2]) / 3.0;
	for (int k = 0; k < 3; k++) {
		Point2D<double> dir = boundaryVertices[(k + 1) % 3] - boundaryVertices[k];
		Point2D<double> mid = (boundaryVertices[(k + 1) % 3] + boundaryVertices[k]) / 2.0;
		Point2D<double> normal(-dir[1], dir[0]);
		normal /= Point2D<double>::Length(normal);
		double offset = Point2D<double>::Dot(mid, normal);

		if (Point2D<double>::Dot(boundaryCenter, normal) < offset) {
			normal *= -1.0;
			offset *= -1.0;
		}
		offset -= 1e-12;
		Point3D<double> boundaryLine(normal[0], normal[1], offset);

		//Check if boundary line match inputLine
		bool colinearSegment = false;
		for (int v = 0; v < inputLines.size(); v++) {
			if (fabs(inputLines[v][0] - boundaryLine[0]) < colinearError && fabs(inputLines[v][1] - boundaryLine[1]) < colinearError && fabs(inputLines[v][2] - boundaryLine[2]) < colinearError) {
				//printf("Found colinearity! \n");
				colinearSegment = true;
			}
		}
		if (!colinearSegment) {
			ClipConvexPolygon(vertices, normal, offset);
		}
	}

	return vertices.size();
}

//Vertex type 
// -2 non initialized
// -1 exterior
// 0 on edge
// 1 interior

int  ClipPartiallyIndexedPolygonToIndexedEdge(AtlasIndexedPolygon &  polygon, const Point2D<double> & normal, double & offset, const int edgeIndex, int atlasVertexIndices[2]){

	std::vector<Point2D<double>> outputVertices;
	std::vector<int> outputVertexIndices; 
	std::vector<int> outputEdgeIndices;
	std::vector<int> outputParentVertexEdgeIndices;
	int lastAddedVertexType = -2;

	Point2D<double> previousVertex = polygon.vertices[polygon.vertices.size() - 1];
	int previousVertexIndex = polygon.atlasVertexIndices[polygon.vertices.size() - 1];
	int previousEdgeIndex = polygon.atlasEdgeIndices[polygon.vertices.size() - 2];
	int nextEdgeIndex = polygon.atlasEdgeIndices[polygon.vertices.size() - 1];
	int previousVertexEdgeSupport = polygon.atlasVertexParentEdge[polygon.vertices.size() - 1];

	int previousVertexType = -2;
	double previousLevel;

	bool emptyPolygon = true;

	bool isCurrentVertexACorner = previousVertexIndex != -1 && (previousVertexIndex == atlasVertexIndices[0] || previousVertexIndex == atlasVertexIndices[1]);
	//bool isPreviousEdgeColinear = previousEdgeIndex != -1 && previousEdgeIndex == edgeIndex;
   // bool isNextEdgeColinear = nextEdgeIndex != -1 && nextEdgeIndex == edgeIndex;
	bool isOnTheEdge = previousVertexEdgeSupport != -1 && previousVertexEdgeSupport == edgeIndex;
	if (isOnTheEdge || isCurrentVertexACorner){
		previousVertexType = 0;
		previousLevel = 0;
	}
	else{
		previousLevel = Point2D<double>::Dot(previousVertex, normal) - offset;
		previousVertexType = previousLevel > 0 ? 1 : -1;
	}

	for (int i = 0; i < polygon.vertices.size(); i++){
		Point2D<double> currentVertex = polygon.vertices[i];
		int currentVertexIndex = polygon.atlasVertexIndices[i];
		int currentVertexType = -2;
		int currentVertexEdgeSupport = polygon.atlasVertexParentEdge[i];
		double currentLevel;

		isCurrentVertexACorner = currentVertexIndex != -1 && (currentVertexIndex == atlasVertexIndices[0] || currentVertexIndex == atlasVertexIndices[1]);
		isOnTheEdge = currentVertexEdgeSupport != -1 && currentVertexEdgeSupport == edgeIndex;
		//isPreviousEdgeColinear = isNextEdgeColinear;
		previousEdgeIndex = nextEdgeIndex;

		nextEdgeIndex = polygon.atlasEdgeIndices[i];
		//isNextEdgeColinear = nextEdgeIndex != -1 && nextEdgeIndex == edgeIndex;

		if (isOnTheEdge || isCurrentVertexACorner){
			currentVertexType = 0;
			currentLevel = 0;
		}
		else{
			currentLevel = Point2D<double>::Dot(currentVertex, normal) - offset;
			currentVertexType = currentLevel > 0 ? 1 : -1;
		}

		if (previousVertexType == -1){
			if (currentVertexType == -1){
				//Do nothing
			}
			else if (currentVertexType == -0){
				//Do nothing
			}
			else{//Entrying edge
				double alpha = -currentLevel / (previousLevel - currentLevel);
				Point2D<double> intersection = previousVertex*alpha + currentVertex*(1.0 - alpha);
				outputVertices.push_back(intersection);
				outputVertexIndices.push_back(-1);
				outputEdgeIndices.push_back(previousEdgeIndex);
				outputParentVertexEdgeIndices.push_back(edgeIndex);
				//lastAddedVertexType = 0;
			}
		}
		else if (previousVertexType == 0){
			outputVertices.push_back(previousVertex);
			outputVertexIndices.push_back(previousVertexIndex);
			outputParentVertexEdgeIndices.push_back(previousVertexEdgeSupport);
			//lastAddedVertexType = 0;
			if (currentVertexType < 0){
				outputEdgeIndices.push_back(edgeIndex);
			}
			else if (currentVertexType == 0){
				outputEdgeIndices.push_back(edgeIndex);
			}
			else{
				outputEdgeIndices.push_back(previousEdgeIndex);
			}
		}
		else{
			outputVertices.push_back(previousVertex);
			outputVertexIndices.push_back(previousVertexIndex);
			outputParentVertexEdgeIndices.push_back(previousVertexEdgeSupport);
			emptyPolygon = false;
			//lastAddedVertexType = 1;
			if (currentVertexType < 0){//Exiting edge
				outputEdgeIndices.push_back(previousEdgeIndex);

				double alpha = -currentLevel / (previousLevel - currentLevel);
				Point2D<double> intersection = previousVertex*alpha + currentVertex*(1.0 - alpha);
				outputVertices.push_back(intersection);
				outputVertexIndices.push_back(-1);
				outputEdgeIndices.push_back(edgeIndex);
				outputParentVertexEdgeIndices.push_back(edgeIndex);
				//lastAddedVertexType = 0;
			}
			else if (currentVertexType == 0){
				outputEdgeIndices.push_back(previousEdgeIndex);
			}
			else{
				outputEdgeIndices.push_back(previousEdgeIndex);
			}
		}

		previousVertex = currentVertex;
		previousLevel = currentLevel;
		previousVertexType = currentVertexType;
		previousVertexIndex = currentVertexIndex;
		previousVertexEdgeSupport = currentVertexEdgeSupport;
	}

	if (emptyPolygon){
		polygon.vertices.clear();
		polygon.atlasVertexIndices.clear();
		polygon.atlasEdgeIndices.clear();
		polygon.atlasVertexParentEdge.clear();
	}
	else{
		polygon.vertices = outputVertices;
		polygon.atlasVertexIndices = outputVertexIndices;
		polygon.atlasEdgeIndices = outputEdgeIndices;
		polygon.atlasVertexParentEdge = outputParentVertexEdgeIndices;

		if (polygon.vertices.size() != polygon.atlasVertexIndices.size() || polygon.vertices.size() != polygon.atlasEdgeIndices.size() || polygon.vertices.size() != polygon.atlasVertexParentEdge.size()){
			printf("Polygon array size does not match! \n");
			return 0;
		}

		//Check for non consecutive colinear edges
		for (int i = 0; i < polygon.atlasEdgeIndices.size(); i++){
			if (polygon.atlasEdgeIndices[i] != -1 && polygon.atlasEdgeIndices[i] == polygon.atlasEdgeIndices[(i + 1) % polygon.atlasEdgeIndices.size()]){
				printf("Unexpected consecutive colinear edges! \n");
				return 0;
			}
		}
	}

	return 1;
}

//Only for convex polygons
int ClipPartiallyIndexedPolygonToIndexedTriangle(AtlasIndexedPolygon & polygon, const AtlasIndexedTriangle & triangle, bool verbose = false){

	Point2D<double> triangleCenter = (triangle.vertices[0] + triangle.vertices[1] + triangle.vertices[2]) / 3.0;

	for (int k = 0; k < 3; k++){

		//Compute edge normal (pointing inside)
		Point2D<double> dir = triangle.vertices[(k + 1) % 3] - triangle.vertices[k];
		Point2D<double> mid = (triangle.vertices[(k + 1) % 3] + triangle.vertices[k]) / 2.0;
		Point2D<double> normal(-dir[1], dir[0]);
		normal /= Point2D<double>::Length(normal);
		double offset = Point2D<double>::Dot(mid, normal);
		if (Point2D<double>::Dot(triangleCenter, normal) < offset){
			normal *= -1.0;
			offset *= -1.0;
		}

		int edgeIndex = triangle.atlasEdgeIndices[k];
		int atlasVertexIndices[2] = { triangle.atlasVertexIndices[k], triangle.atlasVertexIndices[(k + 1) % 3] };
		int clippingResult = ClipPartiallyIndexedPolygonToIndexedEdge(polygon, normal, offset, edgeIndex, atlasVertexIndices);
			
		if (verbose){
			printf("polygon after clipping against edge %d, corners %d %d, normal  %f %f, offset %f. \n", edgeIndex, atlasVertexIndices[0], atlasVertexIndices[1],normal[0],normal[1],offset);
			printf("%d \n", polygon.vertices.size());
			for (int v = 0; v < polygon.vertices.size(); v++)printf("%f %f %f \n", polygon.vertices[v][0], polygon.vertices[v][1], 0);
			for (int v = 0; v < polygon.vertices.size(); v++) printf("%d ", polygon.atlasVertexIndices[v]);
			printf("\n");
			for (int v = 0; v < polygon.vertices.size(); v++) printf("%d ", polygon.atlasEdgeIndices[v]);
			printf("\n");
		}

		if (!clippingResult){
			printf("Clipping failed! \n");
			return -1;
		}
	}

	return polygon.vertices.size();
}

void SetAtlasIndexedPolygonFromTriangle(const AtlasIndexedTriangle & triangle, AtlasIndexedPolygon & polygon){
	for (int k = 0; k < 3; k++){
		polygon.vertices.push_back(triangle.vertices[k]);
		polygon.indices.push_back(triangle.indices[k]);
		polygon.atlasVertexIndices.push_back(triangle.atlasVertexIndices[k]);
		polygon.atlasEdgeIndices.push_back(triangle.atlasEdgeIndices[k]);
		polygon.atlasVertexParentEdge.push_back(triangle.atlasVertexParentEdge[k]);
	}
}

//Points in general positions
int  ClipIndexedIntersectionPolygonToIndexedIntersectionEdge(IndexedIntersectionPolygon &  polygon, const Point2D<double> & normal, double & offset, const int edgeIndex){

	std::vector<Point2D<double>> outputVertices;
	std::vector<unsigned long long> outputIndices;
	std::vector<int> outputEdgeIndices;

	Point2D<double> previousVertex = polygon.vertices[polygon.vertices.size() - 1];
	unsigned long long previousVertexIndex = polygon.indices[polygon.vertices.size() - 1];
	int previousEdgeIndex = polygon.edgeIndices[polygon.vertices.size() - 2];
	int nextEdgeIndex = polygon.edgeIndices[polygon.vertices.size() - 1];
	
	double previousLevel = Point2D<double>::Dot(previousVertex, normal) - offset;
	bool isPreviousInterior = previousLevel > 0;

	for (int i = 0; i < polygon.vertices.size(); i++){
		Point2D<double> currentVertex = polygon.vertices[i];
		unsigned long long currentVertexIndex = polygon.indices[i];
		previousEdgeIndex = nextEdgeIndex;
		nextEdgeIndex = polygon.edgeIndices[i];
		
		double currentLevel = Point2D<double>::Dot(currentVertex, normal) - offset;
		bool isCurrentInterior = currentLevel > 0;
		
		if (!isPreviousInterior){
			if (isCurrentInterior){//Entrying edge
				double alpha = -currentLevel / (previousLevel - currentLevel);
				Point2D<double> intersection = previousVertex*alpha + currentVertex*(1.0 - alpha);
				outputVertices.push_back(intersection);

				unsigned long long intersectionVertexKey = SetIntersectionKey(edgeIndex, previousEdgeIndex);
				outputIndices.push_back(intersectionVertexKey);

				outputEdgeIndices.push_back(previousEdgeIndex);
			}
		}
		else{
			outputVertices.push_back(previousVertex);
			outputIndices.push_back(previousVertexIndex);
			outputEdgeIndices.push_back(previousEdgeIndex);

			if (!isCurrentInterior){//Exiting edge
				double alpha = -currentLevel / (previousLevel - currentLevel);
				Point2D<double> intersection = previousVertex*alpha + currentVertex*(1.0 - alpha);
				outputVertices.push_back(intersection);
				unsigned long long intersectionVertexKey = SetIntersectionKey(edgeIndex, previousEdgeIndex);
				outputIndices.push_back(intersectionVertexKey);
				outputEdgeIndices.push_back(edgeIndex);
			}
		}
		previousVertex = currentVertex;
		previousLevel = currentLevel;
		previousVertexIndex = currentVertexIndex;
		isPreviousInterior = isCurrentInterior;
	}

	polygon.vertices = outputVertices;
	polygon.indices = outputIndices;
	polygon.edgeIndices = outputEdgeIndices;

	return 1;
}


int ClipIndexedIntersectionPolygonToIndexedIntersectionTriangle(IndexedIntersectionPolygon & polygon, const IndexedIntersectionTriangle & triangle, bool verbose = false){
	
	Point2D<double> triangleCenter = (triangle.vertices[0] + triangle.vertices[1] + triangle.vertices[2]) / 3.0;
	unsigned long long cornerKeys[6];
	for (int k = 0; k < 3; k++){
		cornerKeys[2 * k] = SetIntersectionKey(triangle.edgeIndices[k], triangle.edgeIndices[(k + 2) % 3]);
		cornerKeys[2 * k + 1] = SetIntersectionKey(triangle.edgeIndices[(k + 2) % 3], triangle.edgeIndices[k]);
	}
	unsigned long long cornerIndices[6];
	for (int k = 0; k < 3; k++){
		cornerIndices[2 * k] = cornerIndices[2 * k + 1] = triangle.indices[k];
	}
	
	if (verbose){
		printf("Triangle vertices\n");
		for (int v = 0; v < 3; v++)printf("%f %f %f \n", triangle.vertices[v][0], triangle.vertices[v][1], 0);
		printf("Triangle corner keys \n");
		printf("%llu %llu %llu\n", cornerKeys[0], cornerKeys[1], cornerKeys[2]);
		printf("Triangle corner indices \n");
		printf("%llu %llu %llu\n", cornerIndices[0], cornerIndices[1], cornerIndices[2]);
	}

	bool reverseOrientation = false;

	for (int k = 0; k < 3; k++){
		//Compute edge normal (pointing inside)
		Point2D<double> dir = triangle.vertices[(k + 1) % 3] - triangle.vertices[k];
		Point2D<double> mid = (triangle.vertices[(k + 1) % 3] + triangle.vertices[k]) / 2.0;
		Point2D<double> normal(-dir[1], dir[0]);
		normal /= Point2D<double>::Length(normal);
		double offset = Point2D<double>::Dot(mid, normal);
		if (Point2D<double>::Dot(triangleCenter, normal) < offset){
			normal *= -1.0;
			offset *= -1.0;
			reverseOrientation = true;
		}

		int edgeIndex = triangle.edgeIndices[k];
		int clippingResult = ClipIndexedIntersectionPolygonToIndexedIntersectionEdge(polygon, normal, offset, edgeIndex);

		if (verbose){
			printf("polygon after clipping against edge %d, normal  %f %f, offset %f. \n",edgeIndex,normal[0], normal[1], offset);
			printf("%d \n", polygon.vertices.size());
			for (int v = 0; v < polygon.vertices.size(); v++)printf("%f %f %f \n", polygon.vertices[v][0], polygon.vertices[v][1], 0);
			printf(" Vertex Indices \n");
			for (int v = 0; v < polygon.vertices.size(); v++)printf("%llu ", polygon.indices[v]);
			printf("\n");
			printf(" Edge Indices \n");
			for (int v = 0; v < polygon.vertices.size(); v++) printf("%d ", polygon.edgeIndices[v]);
			printf("\n");
		}

		//if (!clippingResult){
		//	printf("Clipping failed! \n");
		//	return -1;
		//}
	}

	//Identify triangle corners
	for (int i = 0; i < polygon.indices.size(); i++) for (int k = 0; k < 6; k++) if (polygon.indices[i] == cornerKeys[k]) polygon.indices[i] = cornerIndices[k];

	if (reverseOrientation){
		int n = polygon.vertices.size();
		std::vector < Point2D < double >> reversedVertices(n);
		std::vector < unsigned long long > reversedIndices(n);
		std::vector < int > reversedEdges(n);
		for (int k = 0; k < n; k++){
			reversedVertices[k] = polygon.vertices[n - 1 - k];
			reversedIndices[k] = polygon.indices[n - 1 - k];
			reversedEdges[k] = polygon.edgeIndices[(2*n - k - 2)% n];
		}
		polygon.vertices = reversedVertices;
		polygon.indices = reversedIndices;
		polygon.edgeIndices = reversedEdges;
	}

	return polygon.vertices.size();
}

void GetTriangleIntegerBBox(Point2D< double > tPos[3], const double invCellSizeW, const double invCellSizeH, int minCorner[2], int maxCorner[2]) {
	double fminx = std::min< double >(std::min< double >(tPos[0][0], tPos[1][0]), tPos[2][0]);
	fminx = std::max< double >(fminx, 0.f);
	double fminy = std::min< double >(std::min< double >(tPos[0][1], tPos[1][1]), tPos[2][1]);
	fminy = std::max< double >(fminy, 0.f);
	double fmaxx = std::max< double >(std::max< double >(tPos[0][0], tPos[1][0]), tPos[2][0]);
	fmaxx = std::min< double >(fmaxx, 1.f);
	double fmaxy = std::max< double >(std::max< double >(tPos[0][1], tPos[1][1]), tPos[2][1]);
	fmaxy = std::min< double >(fmaxy, 1.f);

	minCorner[0] = static_cast<int>(floor(double(fminx)*invCellSizeW));
	minCorner[1] = static_cast<int>(floor(double(fminy)*invCellSizeH));
	maxCorner[0] = static_cast<int>(ceil(double(fmaxx)*invCellSizeW));
	maxCorner[1] = static_cast<int>(ceil(double(fmaxy)*invCellSizeH));
}

void GetEdgeIntegerBBox(Point2D< double > tPos[3], const double invCellSizeW, const double invCellSizeH, int minCorner[2], int maxCorner[2]) {
	double fminx = std::min< double >(tPos[0][0], tPos[1][0]);
	fminx = std::max< double >(fminx, 0.f);
	double fminy = std::min< double >(tPos[0][1], tPos[1][1]);
	fminy = std::max< double >(fminy, 0.f);
	double fmaxx = std::max< double >(tPos[0][0], tPos[1][0]);
	fmaxx = std::min< double >(fmaxx, 1.f);
	double fmaxy = std::max< double >(tPos[0][1], tPos[1][1]);
	fmaxy = std::min< double >(fmaxy, 1.f);

	minCorner[0] = static_cast<int>(floor(double(fminx)*invCellSizeW));
	minCorner[1] = static_cast<int>(floor(double(fminy)*invCellSizeH));
	maxCorner[0] = static_cast<int>(ceil(double(fmaxx)*invCellSizeW));
	maxCorner[1] = static_cast<int>(ceil(double(fmaxy)*invCellSizeH));
}

SquareMatrix< double, 2 > GetBaricentricMap(Point2D< double > tPos[3]) {
	SquareMatrix< double, 2 > parametrizationMap;
	parametrizationMap.coords[0][0] = tPos[1][0] - tPos[0][0];
	parametrizationMap.coords[0][1] = tPos[1][1] - tPos[0][1];
	parametrizationMap.coords[1][0] = tPos[2][0] - tPos[0][0];
	parametrizationMap.coords[1][1] = tPos[2][1] - tPos[0][1];
	return parametrizationMap.inverse();
}


#endif //POLYGON_CLIPPING_INCLUDED

