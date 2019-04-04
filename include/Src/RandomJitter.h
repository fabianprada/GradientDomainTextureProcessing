#pragma once
#include <set>
class SimpleIndexedVector2D {
public:
	SimpleIndexedVector2D(Point2D<double> p_p, int p_index) {
		p = p_p;
		index = p_index;
	}
	Point2D<double> p;
	int index;
};

class SimpleIndexedVector2DComparison
{
public:
	bool operator()(const SimpleIndexedVector2D & p1, const SimpleIndexedVector2D & p2) const {
		for (int i = 0; i < 2; i++) {
			if (p1.p[i] < p2.p[i])
				return true;
			else if (p2.p[i] < p1.p[i])
				return false;
		}
		return false;
	}
};

void AddRandomJitter(std::vector<Point2D<double>> & textureCoordinates, double jitterScale) {
	int lastVertexIndex = 0;
	std::set<SimpleIndexedVector2D, SimpleIndexedVector2DComparison> IndexedPointSet;
	std::set<SimpleIndexedVector2D, SimpleIndexedVector2DComparison>::iterator it;
	std::vector<int> vertexIndex(textureCoordinates.size(), -1);
	for (int t = 0; t < textureCoordinates.size(); t++) {
		int currentCorner = -1;
		SimpleIndexedVector2D idxP(textureCoordinates[t], lastVertexIndex);
		it = IndexedPointSet.find(idxP);
		if (it == IndexedPointSet.end()) {
			IndexedPointSet.insert(idxP);
			currentCorner = lastVertexIndex;
			lastVertexIndex++;
		}
		else {
			SimpleIndexedVector2D indexPoint = *it;
			currentCorner = indexPoint.index;
		}
		vertexIndex[t] = currentCorner;
	}
	srand(time(NULL));
	std::vector<Point2D<double>>randomJitter(lastVertexIndex);
	for (int i = 0; i < randomJitter.size(); i++) {
		randomJitter[i] = Point2D < double >(1.0 - 2.0 * double(rand()) / double(RAND_MAX), 1.0 - 2.0 *  double(rand()) / double(RAND_MAX))*jitterScale;
	}
	for (int i = 0; i < textureCoordinates.size(); i++) {
		textureCoordinates[i] += randomJitter[vertexIndex[i]];
	}
}