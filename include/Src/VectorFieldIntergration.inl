#pragma once
#include "Hierarchy.h"

class InteriorCellLine {
public:
	int prevLineIndex;
	int nextLineIndex;
	int length;
};

int InitializeGridChartInteriorCellLines(const AtlasChart & atlasChart, const GridChart & gridChart, std::vector<InteriorCellLine> & interiorCellLines, std::vector<std::pair<int, int>> & interiorCellLineIndex) {

	const Image<int> & cellType = gridChart.cellType;
	int width = cellType.width();
	int height = cellType.height();

	const Image<int> & globalTexelIndex = gridChart.globalTexelIndex;
	const Image<int> & globalTexelInteriorIndex = gridChart.globalInteriorTexelIndex;

	int localInteriorCellIndex = 0;

	for (int j = 0; j < height; j++) {
		int offset = 0;
		bool previousIsInterior = false;
		int rasterStart = -1;
		while (offset < width) {
			bool currentIsInterior = cellType(offset, j) == 1;
			if ((offset == 0 || offset == width - 1) && currentIsInterior) {
				printf("Unexpected interior cell! \n");
				return 0;
			}
			if (currentIsInterior && !previousIsInterior) rasterStart = offset; //Start raster line
			if (!currentIsInterior && previousIsInterior) { //Terminate raster line
				InteriorCellLine newLine;
				newLine.prevLineIndex = globalTexelIndex(rasterStart, j);
				newLine.nextLineIndex = globalTexelIndex(rasterStart, j + 1);
				newLine.length = offset - rasterStart;

				if (newLine.prevLineIndex == -1 || newLine.nextLineIndex == -1) {
					printf("Inavlid Indexing! \n");
					return 0;
				}
				int currentLine = interiorCellLines.size();

				for (int k = 0; k < offset - rasterStart; k++) {
					if (gridChart.interiorCellCorners[localInteriorCellIndex][0] != globalTexelInteriorIndex(rasterStart + k, j)) {
						printf("Unexpected corner id! \n");
						return 0;
					}

					interiorCellLineIndex.push_back(std::pair<int, int>(currentLine, k));
					localInteriorCellIndex++;
				}

				interiorCellLines.push_back(newLine);
			}
			previousIsInterior = currentIsInterior;
			offset++;
		}
	}

	return 1;
}


int InitializeGridAtladInteriorCellLines(const std::vector<AtlasChart> & atlasCharts, const std::vector<GridChart> & gridCharts, std::vector<InteriorCellLine> & interiorCellLines, std::vector<std::pair<int, int>> & interiorCellLineIndex) {
	for (int i = 0; i < gridCharts.size(); i++) {
		InitializeGridChartInteriorCellLines(atlasCharts[i], gridCharts[i], interiorCellLines, interiorCellLineIndex);
	}
	return 1;
}

template<class Real>
int InitializeChartNormalizedVectorField(const std::vector<SquareMatrix<double, 2>> & parameterMetric, const AtlasChart & atlasChart, const GridChart & gridChart, const std::vector<std::pair<int, int>> & interiorCellLineIndex, const std::vector<int> & fineBoundaryIndex, std::vector<std::vector<CellLineSampleInfo<Real>>> & cellLineSamples, std::vector<HexagonalSampleInfo<Real>> & hexSamples) {

	double cumArea = 0;
	double cumBoundaryMass = 0;

	////Rasterize
	double largestPrecisionError = 0;
	int polygonCount = 0;
	double precision_error = 1e-3;
	//double elementCumArea = 0;
	for (int t = 0; t < atlasChart.triangles.size(); t++) {

		Point2D< double > tPos[3];
		for (int i = 0; i < 3; i++) tPos[i] = atlasChart.vertices[atlasChart.triangles[t][i]] - gridChart.corner;

		SquareMatrix<double, 2> parameter_metric = parameterMetric[t];

		SquareMatrix<double, 2> inverse_parameter_metric = parameter_metric.inverse();
		double parameter_scale_factor = sqrt(parameter_metric.determinant());

		//BBox
		int minCorner[2];
		int maxCorner[2];
		GetTriangleIntegerBBox(tPos, 1.0 / gridChart.cellSizeW, 1.0 / gridChart.cellSizeH, minCorner, maxCorner);

		std::vector<Point2D<double>> parametricVertices(3);
		parametricVertices[0] = tPos[0], parametricVertices[1] = tPos[1], parametricVertices[2] = tPos[2];

		AtlasIndexedTriangle atlasTriangle;
		for (int k = 0; k < 3; k++) {
			atlasTriangle.vertices[k] = tPos[k];
			atlasTriangle.atlasEdgeIndices[k] = atlasChart.atlasEdgeIndices[3 * t + k];
			atlasTriangle.atlasVertexIndices[k] = atlasChart.triangles[t][k];
			atlasTriangle.atlasVertexParentEdge[k] = -1;
		}

		for (int j = minCorner[1]; j < maxCorner[1]; j++) {
			for (int i = minCorner[0]; i < maxCorner[0]; i++) {
				int localInteriorIndex = gridChart.localInteriorCellIndex(i, j);
				int localBoundaryIndex = gridChart.localBoundaryCellIndex(i, j);
				if (localInteriorIndex != -1 && localBoundaryIndex != -1) {
					printf("ERROR: Cell simultaneosly interior and boundary!\n");
					return 0;
				}
				if (localInteriorIndex != -1) {

					std::vector<Point2D<double>> polygonVertices = parametricVertices;
					if (ClipTriangleToPrimalCell(polygonVertices, i, j, gridChart.cellSizeW, gridChart.cellSizeH)) {

						int globalInteriorIndex = localInteriorIndex + gridChart.globalIndexInteriorCellOffset;
						int cellLineId = interiorCellLineIndex[globalInteriorIndex].first;
						int cellLineOffset = interiorCellLineIndex[globalInteriorIndex].second;

						CellLineSampleInfo<Real> cellSample;

						cellSample.cellOffset = cellLineOffset;

						int varIndices[4];
						for (int k = 0; k < 4; k++) varIndices[k] = gridChart.interiorCellCorners[localInteriorIndex][k];

						Point2D<double> polygonBaricenter(0, 0);
						double pArea = 0;
						for (int p = 2; p < polygonVertices.size(); p++) {
							Point2D<double> dm[2] = { polygonVertices[p - 1] - polygonVertices[0], polygonVertices[p] - polygonVertices[0] };
							double tArea = fabs(dm[0][0] * dm[1][1] - dm[0][1] * dm[1][0]) / 2.0;
							polygonBaricenter += tArea * (polygonVertices[0] + polygonVertices[p - 1] + polygonVertices[p]) / 3.0;
							pArea += tArea;
						}
						polygonBaricenter /= pArea;

						polygonBaricenter[0] /= gridChart.cellSizeW;
						polygonBaricenter[1] /= gridChart.cellSizeH;
						polygonBaricenter[0] -= (double(i));
						polygonBaricenter[1] -= (double(j));
						if (polygonBaricenter[0] < 0 - precision_error || polygonBaricenter[0] > 1 + precision_error || polygonBaricenter[1] < 0 - precision_error || polygonBaricenter[1] > 1 + precision_error) {
							printf("Sample out of unit box! (%f %f)\n", polygonBaricenter[0], polygonBaricenter[1]);
							return 0;
						}

						cellSample.pos = Point2D<Real>(polygonBaricenter[0], polygonBaricenter[1]);
						SquareMatrix<Real, 2> tensor;
						tensor(0, 0) = inverse_parameter_metric(0, 0);
						tensor(0, 1) = inverse_parameter_metric(0, 1);
						tensor(1, 0) = inverse_parameter_metric(1, 0);
						tensor(1, 1) = inverse_parameter_metric(1, 1);
						cellSample.tensor = tensor;
						for (int k = 0; k < 4; k++) {
							cellSample.v[k] = Point2D<Real>(0, 0);
						}

						for (int p = 2; p < polygonVertices.size(); p++) {

							Point2D<double> dm[2] = { polygonVertices[p - 1] - polygonVertices[0], polygonVertices[p] - polygonVertices[0] };

							SquareMatrix<double, 2> unit_map_differential;
							unit_map_differential(0, 0) = dm[0][0];
							unit_map_differential(0, 1) = dm[0][1];

							unit_map_differential(1, 0) = dm[1][0];
							unit_map_differential(1, 1) = dm[1][1];

							double unit_map_scale_factor = fabs(unit_map_differential.determinant());

							Point2D<double> cellSamples[6]; //Belong to [0,1] x [0,1]
							for (int s = 0; s < 6; s++) {
								cellSamples[s] = polygonVertices[0] + dm[0] * integrator4_samplePos[s][0] + dm[1] * integrator4_samplePos[s][1];
								cellSamples[s][0] /= gridChart.cellSizeW;
								cellSamples[s][1] /= gridChart.cellSizeH;
								cellSamples[s][0] -= (double(i));
								cellSamples[s][1] -= (double(j));
								if (cellSamples[s][0] < 0 - precision_error || cellSamples[s][0] > 1 + precision_error || cellSamples[s][1] < 0 - precision_error || cellSamples[s][1] > 1 + precision_error) {
									printf("Sample out of unit box! (%f %f)\n", cellSamples[s][0], cellSamples[s][1]);
									return 0;
								}
							}

							//elementCumArea += parameter_scale_factor*unit_map_scale_factor / 2.0;
							cumArea += parameter_scale_factor*unit_map_scale_factor / 2.0;


							//Intergrate gradient field
							Point2D<double> sampleGradients[4][6];
							for (int k = 0; k < 4; k++) {
								for (int s = 0; s < 6; s++) {
									sampleGradients[k][s] = BilinearElementGradient(k, cellSamples[s]);
									sampleGradients[k][s][0] /= gridChart.cellSizeW;
									sampleGradients[k][s][1] /= gridChart.cellSizeH;
									sampleGradients[k][s] = inverse_parameter_metric * sampleGradients[k][s];
									//sampleGradients[k][s] = BilinearElementGradient_Scaled(k, cellSamples[s],width,height);
								}
							}


							for (int k = 0; k < 4; k++) {
								double integral = 0;
								for (int s = 0; s < 6; s++) {
									Point2D<double> v_d = (parameter_metric * sampleGradients[k][s]) * integrator4_sampleWeight[s] * parameter_scale_factor * unit_map_scale_factor / 2.0;
									cellSample.v[k] += Point2D<Real>(v_d[0], v_d[1]);
								}
							}
						}
						cellLineSamples[cellLineId].push_back(cellSample);
					}
				}
				else if (localBoundaryIndex != -1) {
					double fine_precision_error = 1e-4;

					std::vector<Point2D<double>> cellClippedVertices = parametricVertices;
					double cellClippedArea = 0;
					if (ClipTriangleToPrimalCell(cellClippedVertices, i, j, gridChart.cellSizeW, gridChart.cellSizeH)) {
						for (int p = 2; p < cellClippedVertices.size(); p++) {
							Point2D<double> dm[2] = { cellClippedVertices[p - 1] - cellClippedVertices[0], cellClippedVertices[p] - cellClippedVertices[0] };

							SquareMatrix<double, 2> unit_map_differential;
							unit_map_differential(0, 0) = dm[0][0];
							unit_map_differential(0, 1) = dm[0][1];

							unit_map_differential(1, 0) = dm[1][0];
							unit_map_differential(1, 1) = dm[1][1];

							double unit_map_scale_factor = fabs(unit_map_differential.determinant());

							cellClippedArea += parameter_scale_factor*unit_map_scale_factor / 2.0;
						}
					}

					double cumPolygonArea = 0;
					std::vector<BoundaryIndexedTriangle> cellBoundaryTriangles = gridChart.boundaryTriangles[localBoundaryIndex];
					for (int bt = 0; bt < cellBoundaryTriangles.size(); bt++) {
						BoundaryIndexedTriangle atlasBoundaryTriangle = cellBoundaryTriangles[bt];
						std::vector<Point2D<double>> boundaryVertices(3);
						boundaryVertices[0] = atlasBoundaryTriangle.vertices[0];
						boundaryVertices[1] = atlasBoundaryTriangle.vertices[1];
						boundaryVertices[2] = atlasBoundaryTriangle.vertices[2];
						int boundaryTriangleId = atlasBoundaryTriangle.id;

						AtlasIndexedPolygon atlasBoundaryPolygon;
						SetAtlasIndexedPolygonFromBoundaryTriangle(atlasBoundaryTriangle, atlasBoundaryPolygon);
						int clippingResult = ClipPartiallyIndexedPolygonToIndexedTriangle(atlasBoundaryPolygon, atlasTriangle);
						if (clippingResult > 0) {

							HexagonalSampleInfo<Real> hexagonalSample;

							int varIndices[6];
							const HexagonalIndex & hexIndices = gridChart.boundaryTriangles[localBoundaryIndex][bt].indices;
							for (int k = 0; k < 6; k++) {
								int _fineBoundaryIndex = fineBoundaryIndex[hexIndices[k]];
								if (_fineBoundaryIndex != -1) {
									hexagonalSample.fineNodes[k] = _fineBoundaryIndex;
									hexagonalSample.v[k] = Point2D<Real>(0, 0);
								}
								else {
									printf("ERROR:Invalid fine boundary index! \n");
									return 0;
								}
							}


							SquareMatrix< double, 2 > boundaryMapDifferential;
							boundaryMapDifferential(0, 0) = boundaryVertices[1][0] - boundaryVertices[0][0];
							boundaryMapDifferential(0, 1) = boundaryVertices[1][1] - boundaryVertices[0][1];
							boundaryMapDifferential(1, 0) = boundaryVertices[2][0] - boundaryVertices[0][0];
							boundaryMapDifferential(1, 1) = boundaryVertices[2][1] - boundaryVertices[0][1];

							SquareMatrix< double, 2 > barycentricMap = boundaryMapDifferential.inverse();

							SquareMatrix< double, 2 > boundaryMetric = boundaryMapDifferential.transpose() * parameter_metric * boundaryMapDifferential;
							SquareMatrix< double, 2 > boundaryMetricInverse = boundaryMetric.inverse();
							SquareMatrix<double, 6> polygonStiffness;
							SquareMatrix<double, 6> polygonMass;

							double polygonArea = 0;

							Point2D<double> polygonBaricenter(0, 0);
							double pArea = 0;
							for (int p = 2; p < atlasBoundaryPolygon.vertices.size(); p++) {
								Point2D<double> dm[2] = { atlasBoundaryPolygon.vertices[p - 1] - atlasBoundaryPolygon.vertices[0], atlasBoundaryPolygon.vertices[p] - atlasBoundaryPolygon.vertices[0] };
								double tArea = fabs(dm[0][0] * dm[1][1] - dm[0][1] * dm[1][0]) / 2.0;
								polygonBaricenter += tArea * (atlasBoundaryPolygon.vertices[0] + atlasBoundaryPolygon.vertices[p - 1] + atlasBoundaryPolygon.vertices[p]) / 3.0;
								pArea += tArea;
							}
							polygonBaricenter /= pArea;
							polygonBaricenter = barycentricMap*(polygonBaricenter - boundaryVertices[0]);

							if (polygonBaricenter[0] < 0 - precision_error || polygonBaricenter[1] < 0 - precision_error || (polygonBaricenter[0] + polygonBaricenter[1]) > 1 + precision_error) {
								printf("Sample out of unit right triangle! (%f %f)\n", polygonBaricenter[0], polygonBaricenter[1]);
								return 0;
							}

							hexagonalSample.pos = Point2D<Real>(polygonBaricenter[0], polygonBaricenter[1]);
							SquareMatrix<Real, 2> tensor;
							tensor(0, 0) = boundaryMetricInverse(0, 0);
							tensor(0, 1) = boundaryMetricInverse(0, 1);
							tensor(1, 0) = boundaryMetricInverse(1, 0);
							tensor(1, 1) = boundaryMetricInverse(1, 1);
							hexagonalSample.tensor = tensor;

							for (int p = 2; p < atlasBoundaryPolygon.vertices.size(); p++) {
								Point2D<double> dm[2] = { atlasBoundaryPolygon.vertices[p - 1] - atlasBoundaryPolygon.vertices[0], atlasBoundaryPolygon.vertices[p] - atlasBoundaryPolygon.vertices[0] };

								SquareMatrix<double, 2> unit_map_differential;
								unit_map_differential(0, 0) = dm[0][0];
								unit_map_differential(0, 1) = dm[0][1];

								unit_map_differential(1, 0) = dm[1][0];
								unit_map_differential(1, 1) = dm[1][1];

								//double unit_map_scale_factor = fabs(unit_map_differential.determinant());

								SquareMatrix< double, 2 > unitMapMetric = unit_map_differential.transpose() * parameter_metric * unit_map_differential;

								double _triangleMass = sqrt(fabs(unitMapMetric.determinant())) / 2.0;

								if (_triangleMass > 0) {

									polygonArea += _triangleMass;

									SquareMatrix< double, 2 > boundaryBaseToUnitBase = unit_map_differential.inverse() * boundaryMapDifferential;
									Point2D<double> triangleSamples[6]; //Belong to unit right triangles

									for (int s = 0; s < 6; s++) {
										triangleSamples[s] = atlasBoundaryPolygon.vertices[0] + dm[0] * integrator4_samplePos[s][0] + dm[1] * integrator4_samplePos[s][1];
										triangleSamples[s] = barycentricMap*(triangleSamples[s] - boundaryVertices[0]);
										if (triangleSamples[s][0] < 0 - precision_error || triangleSamples[s][1] < 0 - precision_error || (triangleSamples[s][0] + triangleSamples[s][1]) > 1 + precision_error) {
											printf("Sample out of unit right triangle! (%f %f)\n", triangleSamples[s][0], triangleSamples[s][1]);

											//printf("Boundary triangle \n");
											//printf("%d \n", 3);
											//for (int v = 0; v < 3; v++)printf("%f %f %f \n", boundaryVertices[v][0], boundaryVertices[v][1], 0);
											//printf("Clipped polygon \n");
											//printf("%d \n", atlasBoundaryPolygon.vertices.size());
											//for (int v = 0; v < atlasBoundaryPolygon.vertices.size(); v++)printf("%f %f %f \n", atlasBoundaryPolygon.vertices[v][0], atlasBoundaryPolygon.vertices[v][1], 0);
											//std::vector<Point2D<double>> _polygonVertices = parametricVertices;
											//ClipTriangleToTriangle(_polygonVertices, boundaryVertices);
											return 0;
										}
										else {
											triangleSamples[s][0] = std::max<double>(triangleSamples[s][0], 0);
											triangleSamples[s][1] = std::max<double>(triangleSamples[s][1], 0);
											double excess = (triangleSamples[s][0] + triangleSamples[s][1]) - 1;
											if (excess > 0) {
												triangleSamples[s][0] -= (excess / 2);
												triangleSamples[s][1] -= (excess / 2);
											}
										}
									}


									Point2D<double> sampleGradients[6][6];
									for (int k = 0; k < 6; k++) {
										for (int s = 0; s < 6; s++) {
											sampleGradients[k][s] = QuadraticElementGradient(k, triangleSamples[s]);
										}
									}

									for (int k = 0; k < 6; k++) {
										double integral = 0;
										for (int s = 0; s < 6; s++) {
											Point2D<double> v_d = sampleGradients[k][s] * integrator4_sampleWeight[s] * _triangleMass;
											hexagonalSample.v[k] += Point2D<Real>(v_d[0], v_d[1]);
										}
									}

									polygonCount++;
								}
								else {
									printf("WARNING: Element discarded due to zero mass. Triangle %d. Boundary cell %d . Element %d. Sub triangle %d\n", t, localBoundaryIndex, bt, p - 2);

									printf("Atlas triangle \n");
									printf("%d \n", 3);
									for (int v = 0; v < 3; v++)printf("%f %f %f \n", atlasTriangle.vertices[v][0], atlasTriangle.vertices[v][1], 0);
									printf("%d %d %d \n", atlasTriangle.atlasVertexIndices[0], atlasTriangle.atlasVertexIndices[1], atlasTriangle.atlasVertexIndices[2]);
									printf("%d %d %d \n", atlasTriangle.atlasEdgeIndices[0], atlasTriangle.atlasEdgeIndices[1], atlasTriangle.atlasEdgeIndices[2]);
									for (int _bt = 0; _bt < cellBoundaryTriangles.size(); _bt++) {
										printf("ELEMENT %d \n", _bt);

										BoundaryIndexedTriangle _atlasBoundaryTriangle = cellBoundaryTriangles[_bt];
										AtlasIndexedPolygon _atlasBoundaryPolygon;
										SetAtlasIndexedPolygonFromBoundaryTriangle(_atlasBoundaryTriangle, _atlasBoundaryPolygon);

										printf("Boundary triangle \n");
										printf("%d \n", 3);
										for (int v = 0; v < 3; v++)printf("%f %f %f \n", _atlasBoundaryTriangle.vertices[v][0], _atlasBoundaryTriangle.vertices[v][1], 0);
										printf("%d %d %d \n", _atlasBoundaryTriangle.atlasVertexIndices[0], _atlasBoundaryTriangle.atlasVertexIndices[1], _atlasBoundaryTriangle.atlasVertexIndices[2]);
										printf("%d %d %d \n", _atlasBoundaryTriangle.atlasEdgeIndices[0], _atlasBoundaryTriangle.atlasEdgeIndices[1], _atlasBoundaryTriangle.atlasEdgeIndices[2]);

										ClipPartiallyIndexedPolygonToIndexedTriangle(_atlasBoundaryPolygon, atlasTriangle, _bt == bt);
									}
									return 0;
								}
							}

							cumPolygonArea += polygonArea;

							double integratedPolygonMass = 0;
							for (int k = 0; k < 6; k++)for (int l = 0; l < 6; l++)integratedPolygonMass += polygonMass(k, l);
							if (fabs(integratedPolygonMass - polygonArea) > precision_error) {
								printf("out of precision! \n");
								return 0;
							}
							hexSamples.push_back(hexagonalSample);
						}
						else if (clippingResult < 0) {

							return 0;
						}
					}

					cumBoundaryMass += cumPolygonArea;
					cumArea += cumPolygonArea;
					if (fabs(cumPolygonArea - cellClippedArea) > fine_precision_error) {
						printf("Out of precision! %.12g %.12g\n", cumPolygonArea, cellClippedArea);

						return 0;
					}
				}
			}
		}
	}

	if (0) printf("Largest precision error %g.\n", largestPrecisionError);

	return 1;
}



template<class Real>
int InitializeNormalizedVectorFieldIntegration(const std::vector<std::vector<SquareMatrix<double, 2>>> & parameterMetric, const std::vector<AtlasChart> atlasCharts, const std::vector<GridChart> & gridCharts, const std::vector<std::pair<int, int>> & interiorCellLineIndex, const std::vector<int> & fineBoundaryIndex, std::vector<std::vector<CellLineSampleInfo<Real>>> & cellLineSamples, std::vector<HexagonalSampleInfo<Real>> & hexSamples) {
	for (int i = 0; i < gridCharts.size(); i++) {
		if (!InitializeChartNormalizedVectorField(parameterMetric[i], atlasCharts[i], gridCharts[i], interiorCellLineIndex, fineBoundaryIndex, cellLineSamples, hexSamples)) {
			return 0;
		}
	}
	return 1;
}


template <class Real>
int IntergrateNormalizedVectorFieldFromSamples(const Real invCellSizeW, const Real invCellSizeH, const std::vector<InteriorCellLine> & interiorCellLines, const std::vector<std::vector<CellLineSampleInfo<Real>>> & cellLineSamples, const std::vector<HexagonalSampleInfo<Real>> & hexSamples, const std::vector<Real> & potential, std::vector<Real> & rhs, const std::vector<Real> & boundary_potential, std::vector<Real> & boundary_rhs, bool verbose = false) {

	memset(&rhs[0], 0, rhs.size() * sizeof(Real));
	memset(&boundary_rhs[0], 0, boundary_rhs.size() * sizeof(Real));
	clock_t begin = clock();

	auto UpdateRow = [&](int r)
	{
		const Real * _inPrevious = &potential[interiorCellLines[r].prevLineIndex];
		const Real * _inNext = &potential[interiorCellLines[r].nextLineIndex];

		Real * _outPrevious = &rhs[interiorCellLines[r].prevLineIndex];
		Real * _outNext = &rhs[interiorCellLines[r].nextLineIndex];

		Real cornerValues[4];
		cornerValues[0] = *_inPrevious;
		_inPrevious++;
		cornerValues[1] = *_inPrevious;
		cornerValues[3] = *_inNext;
		_inNext++;
		cornerValues[2] = *_inNext;

		Real rhsValues[4] = { 0, 0, 0, 0 };

		int numSamples = cellLineSamples[r].size();
		const CellLineSampleInfo<Real> * sample = &cellLineSamples[r][0];
		int currentOffset = 0;

		for (int j = 0; j < numSamples; j++) {
			if (sample->cellOffset >currentOffset) {

				cornerValues[0] = cornerValues[1];
				_inPrevious++;
				cornerValues[1] = *_inPrevious;
				cornerValues[3] = cornerValues[2];
				_inNext++;
				cornerValues[2] = *_inNext;

				*_outPrevious += rhsValues[0];
				_outPrevious++;
				*_outNext += rhsValues[3];
				_outNext++;
				rhsValues[0] = rhsValues[1];
				rhsValues[1] = 0;
				rhsValues[3] = rhsValues[2];
				rhsValues[2] = 0;

				currentOffset++;
			}
			Point2D<Real> pos = sample->pos;
			Point2D<Real> sampledDiffusedGradient = -(Point2D<Real>((pos[1] - 1.0)*cornerValues[0] + (1.0 - pos[1])*cornerValues[1] + pos[1] * cornerValues[2] - pos[1] * cornerValues[3], (pos[0] - 1.0)*cornerValues[0] - pos[0] * cornerValues[1] + pos[0] * cornerValues[2] + (1.0 - pos[0])*cornerValues[3])); //Negative gradient
			sampledDiffusedGradient[0] *= invCellSizeW;
			sampledDiffusedGradient[1] *= invCellSizeH;
			Point2D<Real> aux = sample->tensor * sampledDiffusedGradient;
			double len = Point2D<Real>::Dot(sampledDiffusedGradient, aux);
			if (len > 0) {
				len = sqrt(len);
				sampledDiffusedGradient = aux / len;
			}
			else {
				if (0)printf("WARNING: Zero vector! \n");
				sampledDiffusedGradient *= 0;
			}

			rhsValues[0] += Point2D<Real>::Dot(sampledDiffusedGradient, sample->v[0]);
			rhsValues[1] += Point2D<Real>::Dot(sampledDiffusedGradient, sample->v[1]);
			rhsValues[2] += Point2D<Real>::Dot(sampledDiffusedGradient, sample->v[2]);
			rhsValues[3] += Point2D<Real>::Dot(sampledDiffusedGradient, sample->v[3]);
			sample++;
		}

		*_outPrevious += rhsValues[0];
		_outPrevious++;
		*_outPrevious += rhsValues[1];
		*_outNext += rhsValues[3];
		_outNext++;
		*_outNext += rhsValues[2];
	};

	//for (int r = 0; r < interiorCellLines.size(); r++) UpdateRow(r);

	int threads = omp_get_max_threads();
	std::vector<int> lineRange(threads + 1);
	int blockSize = interiorCellLines.size() / threads;
	for (int t = 0; t < threads; t++) lineRange[t] = t*blockSize;
	lineRange[threads] = interiorCellLines.size();
#pragma omp parallel for
	for (int t = 0; t < threads; t++) {
		const int tId = omp_get_thread_num();
		const int firstLine = lineRange[tId];
		const int lastLine = lineRange[tId + 1];
		for (int r = firstLine; r < lastLine; r++) UpdateRow(r);
	}

	if (verbose)printf("Integrating cells %.4f \n", double(clock() - begin) / CLOCKS_PER_SEC);

	begin = clock();
	for (int i = 0; i < hexSamples.size(); i++) {
		const HexagonalSampleInfo<Real> & sample = hexSamples[i];
		Real cornerValues[6] = { boundary_potential[sample.fineNodes[0]], boundary_potential[sample.fineNodes[1]], boundary_potential[sample.fineNodes[2]], boundary_potential[sample.fineNodes[3]], boundary_potential[sample.fineNodes[4]], boundary_potential[sample.fineNodes[5]] };

		Point2D<Real> pos = sample.pos;

		Point2D<Real> sampledDiffusedGradient =
			-(Point2D<Real>(4 * pos[0] + 4 * pos[1] - 3.0, 4 * pos[0] + 4 * pos[1] - 3.0) *  cornerValues[0]
			+ Point2D<Real>(4 * pos[0] - 1.0, 0.0) * cornerValues[1]
			+ Point2D<Real>(0, 4 * pos[1] - 1.0) *  cornerValues[2]
			+ Point2D<Real>(4 * pos[1], 4 * pos[0]) *  cornerValues[3]
			+ Point2D<Real>(-4 * pos[1], -4 * pos[0] - 8 * pos[1] + 4.0) * cornerValues[4]
			+ Point2D<Real>(-8 * pos[0] - 4 * pos[1] + 4, -4 * pos[0]) * cornerValues[5]);

		Point2D<Real> aux = sample.tensor * sampledDiffusedGradient;
		double len = Point2D<Real>::Dot(sampledDiffusedGradient, aux);

		if (len > 0) {
			len = sqrt(len);
			sampledDiffusedGradient = aux / len;
		}
		else {
			if (0) printf("WARNING: Zero vector! \n");
			sampledDiffusedGradient *= 0;
		}

		for (int k = 0; k < 6; k++) {
			boundary_rhs[sample.fineNodes[k]] += Point2D<Real>::Dot(sampledDiffusedGradient, sample.v[k]);
		}
	}
	if (verbose) printf("Integrating hexagons %.4f \n", double(clock() - begin) / CLOCKS_PER_SEC);
	return 1;
}
