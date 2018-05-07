#pragma once
#include "Hierarchy.h"
#include "PolygonClipping.h"

int InitializeChartMassAndStiffness(const std::vector<SquareMatrix<double, 2>> & parameterMetric, const AtlasChart & atlasChart, const GridChart & gridChart, const std::vector<int> & boundaryAndDeepIndex, const std::vector<int> & fineBoundaryIndex, std::vector<double> & deepMassCoefficients, std::vector<double> & deepStiffnessCoefficients,
	std::vector<Eigen::Triplet<double>> & boundaryBoundaryMassTriplets, std::vector<Eigen::Triplet<double>> & boundaryBoundaryStiffnessTriplets,
	std::vector<Eigen::Triplet<double>> & boundaryDeepMassTriplets, std::vector<Eigen::Triplet<double>> & boundaryDeepStiffnessTriplets,
	bool computeCellBasedStiffness, const std::vector<Point3D<double>> & inputSignal, const std::vector<Point3D<double>> & boundarySignal, std::vector<double> & texelToCellCoeffs, std::vector<Eigen::Triplet<double>> boundaryCellStiffnessTriplets[3]) {
	double cumArea = 0;
	double cumBoundaryMass = 0;

	std::vector<SquareMatrix<double, 4>> cellStiffness;
	std::vector<SquareMatrix<double, 4>> cellMass;
	cellMass.resize(gridChart.numInteriorCells);
	cellStiffness.resize(gridChart.numInteriorCells);

	std::vector<SquareMatrix<double, 6>> hexagonalStiffness;
	std::vector<SquareMatrix<double, 6>> hexagonalMass;
	hexagonalMass.resize(gridChart.numBoundaryTriangles);
	hexagonalStiffness.resize(gridChart.numBoundaryTriangles);

	////Rasterize
	int zeroAreaElementCount = 0;
	double largestPrecisionError = 0;
	int polygonCount = 0;
	double precision_error = 1e-3;
	//double elementCumArea = 0;

	//#pragma omp parallel for 
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
					//return 0;
				}
				if (localInteriorIndex != -1) {
					std::vector<Point2D<double>> polygonVertices = parametricVertices;
					if (ClipTriangleToPrimalCell(polygonVertices, i, j, gridChart.cellSizeW, gridChart.cellSizeH)) {
						SquareMatrix<double, 4> polygonStiffness;
						SquareMatrix<double, 4> polygonMass;
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
									//return 0;
								}
							}

							//Integrate scalar product
							double sampleValues[4][6];
							for (int k = 0; k < 4; k++) {
								for (int s = 0; s < 6; s++) {
									sampleValues[k][s] = BilinearElementValue(k, cellSamples[s]);
									//sampleValues[k][s] = BilinearElementValue_Scaled(k, cellSamples[s],width,height);
								}
							}

							double integratedScalarProduct[4][4];
							for (int k = 0; k < 4; k++)for (int l = 0; l < 4; l++) {
								double integral = 0;
								for (int s = 0; s < 6; s++) {
									integral += sampleValues[k][s] * sampleValues[l][s] * integrator4_sampleWeight[s];
								}
								polygonMass(l, k) += integral*parameter_scale_factor*unit_map_scale_factor / 2.0;
								//cellMass[localInteriorIndex](l, k) += integral*parameter_scale_factor*unit_map_scale_factor / 2.0;
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
									//sampleGradients[k][s] = BilinearElementGradient_Scaled(k, cellSamples[s],width,height);
								}
							}

							double integratedGradientProduct[4][4];
							for (int k = 0; k < 4; k++)for (int l = 0; l < 4; l++) {
								double integral = 0;
								for (int s = 0; s < 6; s++) {
									integral += Point2D<double>::Dot(sampleGradients[l][s], inverse_parameter_metric * sampleGradients[k][s])* integrator4_sampleWeight[s] * parameter_scale_factor*unit_map_scale_factor / 2.0;
								}
								polygonStiffness(l, k) += integral;
								//cellStiffness[localInteriorIndex](l, k) += integral;
							}
						}
						//#pragma omp critical
						{
							for (int k = 0; k < 4; k++)for (int l = 0; l < 4; l++) {
								cellMass[localInteriorIndex](l, k) += polygonMass(l, k);
								cellStiffness[localInteriorIndex](l, k) += polygonStiffness(l, k);
							}
						}
					}
				}
				else if (localBoundaryIndex != -1) {
					double fine_precision_error = 1e-4;
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

							for (int p = 2; p < atlasBoundaryPolygon.vertices.size(); p++) {
								Point2D<double> dm[2] = { atlasBoundaryPolygon.vertices[p - 1] - atlasBoundaryPolygon.vertices[0], atlasBoundaryPolygon.vertices[p] - atlasBoundaryPolygon.vertices[0] };

								SquareMatrix<double, 2> unit_map_differential;
								unit_map_differential(0, 0) = dm[0][0];
								unit_map_differential(0, 1) = dm[0][1];

								unit_map_differential(1, 0) = dm[1][0];
								unit_map_differential(1, 1) = dm[1][1];

								//double unit_map_scale_factor = abs(unit_map_differential.determinant());

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

									double values[6][6];
									Point2D<double> gradients[6][6];

									for (int k = 0; k < 6; k++) {
										for (int s = 0; s < 6; s++) {
											values[k][s] = QuadraticElementValue(k, triangleSamples[s]);
											gradients[k][s] = barycentricMap.transpose() * QuadraticElementGradient(k, triangleSamples[s]);
										}
									}

									for (int k = 0; k < 6; k++)for (int l = 0; l < 6; l++) {
										double integral = 0;
										for (int s = 0; s < 6; s++) {
											integral += values[k][s] * values[l][s] * integrator4_sampleWeight[s];
										}
										polygonMass(l, k) += integral*_triangleMass;
									}

									SquareMatrix<double, 6> localStiffness;
									for (int k = 0; k < 6; k++)for (int l = 0; l < 6; l++) {
										double integral = 0;
										for (int s = 0; s < 6; s++) {
											integral += Point2D<double>::Dot(gradients[k][s], inverse_parameter_metric * gradients[l][s]) * integrator4_sampleWeight[s];
										}
										localStiffness(l, k) += integral*_triangleMass;
										polygonStiffness(l, k) += integral*_triangleMass;
									}

									for (int dk = 0; dk < 6; dk++)for (int dl = 0; dl < 6; dl++) {
										double diff = fabs(localStiffness(dk, dl) - localStiffness(dl, dk));
										if (diff > largestPrecisionError) {
											largestPrecisionError = diff;
										}
									}
									polygonCount++;
								}
								else {
									zeroAreaElementCount++;
									printf("WARNING: Zero area polygon at cell %d %d \n", gridChart.cornerCoords[0] + i, gridChart.cornerCoords[1] + j);
								}
							}

							cumPolygonArea += polygonArea;

							double integratedPolygonMass = 0;
							for (int k = 0; k < 6; k++)for (int l = 0; l < 6; l++)integratedPolygonMass += polygonMass(k, l);
							if (fabs(integratedPolygonMass - polygonArea) > precision_error) {
								printf("out of precision! \n");
								//return 0;
							}
							//#pragma omp critical
							{
								for (int dk = 0; dk < 6; dk++)for (int dl = 0; dl < 6; dl++) {
									hexagonalStiffness[boundaryTriangleId](dk, dl) += (polygonStiffness(dk, dl) + polygonStiffness(dl, dk)) / 2.0;
									hexagonalMass[boundaryTriangleId](dk, dl) += (polygonMass(dk, dl) + polygonMass(dl, dk)) / 2.0;
								}
							}
						}
					}

					cumBoundaryMass += cumPolygonArea;
					cumArea += cumPolygonArea;
				}
			}
		}
	}

	if (zeroAreaElementCount) printf("WARNING: Element with zero area = %d \n", zeroAreaElementCount);

	int offset_i[4] = { 0, 1 , 1 ,0 };
	int offset_j[4] = { 0, 0 , 1 ,1 };
	int cellOffset[4] = { 3, 2, 0, 1 };

	auto NeighbourOffset = [&](int k, int l) {
		return  (offset_j[l] - offset_j[k] + 1) * 3 + (offset_i[l] - offset_i[k] + 1);
	};

	for (int i = 0; i < gridChart.interiorCellCorners.size(); i++) {
		const CellIndex & indicesGlobal = gridChart.interiorCellGlobalCorners[i];
		const CellIndex & indicesInterior = gridChart.interiorCellCorners[i];

		const int localCellIndex = gridChart.interiorCellIndexToLocalCellIndex[i];
		const int globalCellIndex = localCellIndex + gridChart.globalIndexCellOffset;

		Point<double, 4> prod[3];
		if (computeCellBasedStiffness){
			Point3D < double > values[4] = { inputSignal[indicesGlobal[0]],inputSignal[indicesGlobal[1]], inputSignal[indicesGlobal[2]], inputSignal[indicesGlobal[3]] };
			for (int c = 0; c < 3; c++) {
				Point<double, 4> v;
				v[0] = values[0][c];
				v[1] = values[1][c];
				v[2] = values[2][c];
				v[3] = values[3][c];

				prod[c] = cellStiffness[i] * v;
			}
		}

		for (int k = 0; k < 4; k++) {
			int currentNode = indicesGlobal[k];
			int _currentBoundaryAndDeepIndex = boundaryAndDeepIndex[currentNode];
			if (_currentBoundaryAndDeepIndex < 0) {//Deep
				int deepIndex = -_currentBoundaryAndDeepIndex - 1;
				for (int l = 0; l < 4; l++) {
					//offset_i(k,l);
					//offset_j(k,l);
					//int neighbourIndex = (neighbourNode.cj - currentNode.cj + 1) * 3 + (neighbourNode.ci - currentNode.ci + 1);
					deepMassCoefficients[10 * deepIndex + NeighbourOffset(k, l)] += cellMass[i](k, l);
					deepStiffnessCoefficients[10 * deepIndex + NeighbourOffset(k, l)] += cellStiffness[i](k, l);
				}
				if (computeCellBasedStiffness) {
					//Add cell data
					texelToCellCoeffs[3 * (4 * deepIndex + cellOffset[k]) + 0] = prod[0][k];
					texelToCellCoeffs[3 * (4 * deepIndex + cellOffset[k]) + 1] = prod[1][k];
					texelToCellCoeffs[3 * (4 * deepIndex + cellOffset[k]) + 2] = prod[2][k];
				}
			}
			else {//Bundary
				int boundaryIndex = _currentBoundaryAndDeepIndex - 1;
				for (int l = 0; l < 4; l++) {
					int neighbourNode = indicesGlobal[l];
					int _neighbourBoundaryAndDeepIndex = boundaryAndDeepIndex[neighbourNode];
					if (_neighbourBoundaryAndDeepIndex < 0) {//Deep
						boundaryDeepMassTriplets.push_back(Eigen::Triplet<double>(boundaryIndex, neighbourNode, cellMass[i](k, l)));
						boundaryDeepStiffnessTriplets.push_back(Eigen::Triplet<double>(boundaryIndex, neighbourNode, cellStiffness[i](k, l)));
					}
					else {//Boundary
						boundaryBoundaryMassTriplets.push_back(Eigen::Triplet<double>(fineBoundaryIndex[indicesInterior[k]], fineBoundaryIndex[indicesInterior[l]], cellMass[i](k, l)));
						boundaryBoundaryStiffnessTriplets.push_back(Eigen::Triplet<double>(fineBoundaryIndex[indicesInterior[k]], fineBoundaryIndex[indicesInterior[l]], cellStiffness[i](k, l)));
					}
				}
				if (computeCellBasedStiffness) {
					//Add cell data
					int _fineBoundaryIndex = fineBoundaryIndex[indicesInterior[k]];
					for (int c = 0; c < 3; c++) {
						boundaryCellStiffnessTriplets[c].push_back(Eigen::Triplet<double>(_fineBoundaryIndex, globalCellIndex, prod[c][k]));
					}
				}
			}

		}
	}

	for (int c = 0; c < gridChart.boundaryTriangles.size(); c++) {

		const int localCellIndex = gridChart.boundaryCellIndexToLocalCellIndex[c];
		const int globalCellIndex = localCellIndex + gridChart.globalIndexCellOffset;

		for (int b = 0; b < gridChart.boundaryTriangles[c].size(); b++) {
			const int i = gridChart.boundaryTriangles[c][b].id;
			const HexagonalIndex & indices = gridChart.boundaryTriangles[c][b].indices;
			HexagonalIndex fineHexIndices;
			for (int k = 0; k < 6; k++) fineHexIndices[k] = fineBoundaryIndex[indices[k]];

			//Add cell data
			
			if (computeCellBasedStiffness){
				Point<double, 6> prod[3];
				Point3D < double > values[6] = { boundarySignal[fineHexIndices[0]],boundarySignal[fineHexIndices[1]],boundarySignal[fineHexIndices[2]],boundarySignal[fineHexIndices[3]],boundarySignal[fineHexIndices[4]],boundarySignal[fineHexIndices[5]] };
				for (int c = 0; c < 3; c++) {
					Point<double, 6> v;
					v[0] = values[0][c];
					v[1] = values[1][c];
					v[2] = values[2][c];
					v[3] = values[3][c];
					v[4] = values[4][c];
					v[5] = values[5][c];
					prod[c] = hexagonalStiffness[i] * v;
				}

				for (int k = 0; k < 6; k++) {
					for (int c = 0; c < 3; c++) {
						boundaryCellStiffnessTriplets[c].push_back(Eigen::Triplet<double>(fineHexIndices[k], globalCellIndex, prod[c][k]));
					}
				}
			}
			for (int k = 0; k < 6; k++)for (int l = 0; l < 6; l++) {
				boundaryBoundaryMassTriplets.push_back(Eigen::Triplet<double>(fineHexIndices[k], fineHexIndices[l], hexagonalMass[i](l, k)));
				boundaryBoundaryStiffnessTriplets.push_back(Eigen::Triplet<double>(fineHexIndices[k], fineHexIndices[l], hexagonalStiffness[i](l, k)));
			}
		}
	}

	if (0) printf("Largest precision error %g.\n", largestPrecisionError);

	return 1;
}

int InitializeMassAndStiffness(const std::vector<std::vector<SquareMatrix<double, 2>>> & parameterMetric, const std::vector<AtlasChart> atlasCharts, const GridAtlas & gridAtlas, const std::vector<int> & fineBoundaryIndex, const int numFineBoundaryNodes,
	std::vector<double> & deepMassCoefficients, std::vector<double> & deepStiffnessCoefficients,
	SparseMatrix<double, int> & boundaryBoundaryMassMatrix, SparseMatrix<double, int> & boundaryBoundaryStiffnessMatrix,
	SparseMatrix<double, int> & boundaryDeepMassMatrix, SparseMatrix<double, int> & boundaryDeepStiffnessMatrix,
	bool computeCellBasedStiffness, const std::vector<Point3D<double>> & inputSignal, const std::vector<Point3D<double>> & boundarySignal, std::vector<double> & texelToCellCoeffs, SparseMatrix<double, int> boundaryCellBasedStiffnessRHSMatrix[3]) {

	const std::vector<GridChart> & gridCharts = gridAtlas.gridCharts;
	const std::vector<int> & boundaryAndDeepIndex = gridAtlas.boundaryAndDeepIndex;

	if(computeCellBasedStiffness) texelToCellCoeffs.resize(3 * 4 * gridAtlas.numDeepTexels);

	clock_t m_begin;
	m_begin = clock();

	std::vector<Eigen::Triplet<double>> boundaryBoundaryMassTriplets;
	std::vector<Eigen::Triplet<double>> boundaryBoundaryStiffnessTriplets;
	std::vector<Eigen::Triplet<double>> boundaryDeepMassTriplets;
	std::vector<Eigen::Triplet<double>> boundaryDeepStiffnessTriplets;
	std::vector<Eigen::Triplet<double>> boundaryCellStiffnessTriplets[3];

#pragma omp parallel for
	for (int i = 0; i < gridCharts.size(); i++) {
		std::vector<Eigen::Triplet<double>> chartBoundaryBoundaryMassTriplets;
		std::vector<Eigen::Triplet<double>> chartBoundaryBoundaryStiffnessTriplets;
		std::vector<Eigen::Triplet<double>> chartBoundaryDeepMassTriplets;
		std::vector<Eigen::Triplet<double>> chartBoundaryDeepStiffnessTriplets;
		std::vector<Eigen::Triplet<double>> chartBoundaryCellStiffnessTriplets[3];
		InitializeChartMassAndStiffness(parameterMetric[i], atlasCharts[i], gridCharts[i], boundaryAndDeepIndex, fineBoundaryIndex, deepMassCoefficients, deepStiffnessCoefficients, chartBoundaryBoundaryMassTriplets, chartBoundaryBoundaryStiffnessTriplets,
			chartBoundaryDeepMassTriplets, chartBoundaryDeepStiffnessTriplets,
			computeCellBasedStiffness ,inputSignal, boundarySignal, texelToCellCoeffs, chartBoundaryCellStiffnessTriplets);
#pragma omp critical
		{
			boundaryBoundaryMassTriplets.insert(boundaryBoundaryMassTriplets.end(), chartBoundaryBoundaryMassTriplets.begin(), chartBoundaryBoundaryMassTriplets.end());
			boundaryBoundaryStiffnessTriplets.insert(boundaryBoundaryStiffnessTriplets.end(), chartBoundaryBoundaryStiffnessTriplets.begin(), chartBoundaryBoundaryStiffnessTriplets.end());
			boundaryDeepMassTriplets.insert(boundaryDeepMassTriplets.end(), chartBoundaryDeepMassTriplets.begin(), chartBoundaryDeepMassTriplets.end());
			boundaryDeepStiffnessTriplets.insert(boundaryDeepStiffnessTriplets.end(), chartBoundaryDeepStiffnessTriplets.begin(), chartBoundaryDeepStiffnessTriplets.end());
			if (computeCellBasedStiffness) for (int c = 0; c < 3; c++) boundaryCellStiffnessTriplets[c].insert(boundaryCellStiffnessTriplets[c].end(), chartBoundaryCellStiffnessTriplets[c].begin(), chartBoundaryCellStiffnessTriplets[c].end());
		}
	}

	printf("Computing Triplets =  %.4f \n", double(clock() - m_begin) / CLOCKS_PER_SEC);


	m_begin = clock();

	Eigen::SparseMatrix<double> _boundaryBoundaryMass;
	_boundaryBoundaryMass.resize(numFineBoundaryNodes, numFineBoundaryNodes);
	_boundaryBoundaryMass.setFromTriplets(boundaryBoundaryMassTriplets.begin(), boundaryBoundaryMassTriplets.end());
	SparseMatrixParser(_boundaryBoundaryMass, boundaryBoundaryMassMatrix, true);

	Eigen::SparseMatrix<double> _boundaryBoundaryStiffness;
	_boundaryBoundaryStiffness.resize(numFineBoundaryNodes, numFineBoundaryNodes);
	_boundaryBoundaryStiffness.setFromTriplets(boundaryBoundaryStiffnessTriplets.begin(), boundaryBoundaryStiffnessTriplets.end());
	SparseMatrixParser(_boundaryBoundaryStiffness, boundaryBoundaryStiffnessMatrix, true);

	int numTexels = gridAtlas.numTexels;
	int numBoundaryTexels = gridAtlas.numBoundaryTexels;

	Eigen::SparseMatrix<double> _boundaryDeepMass;
	_boundaryDeepMass.resize(numBoundaryTexels, numTexels);
	_boundaryDeepMass.setFromTriplets(boundaryDeepMassTriplets.begin(), boundaryDeepMassTriplets.end());
	SparseMatrixParser(_boundaryDeepMass, boundaryDeepMassMatrix, false);

	Eigen::SparseMatrix<double> _boundaryDeepStiffness;
	_boundaryDeepStiffness.resize(numBoundaryTexels, numTexels);
	_boundaryDeepStiffness.setFromTriplets(boundaryDeepStiffnessTriplets.begin(), boundaryDeepStiffnessTriplets.end());
	SparseMatrixParser(_boundaryDeepStiffness, boundaryDeepStiffnessMatrix, false);

	printf("Constructing matrices =  %.4f \n", double(clock() - m_begin) / CLOCKS_PER_SEC);

	if (computeCellBasedStiffness){
		for (int c = 0; c < 3; c++) {
			Eigen::SparseMatrix<double> _boundaryCellBasedRHS;
			_boundaryCellBasedRHS.resize(numFineBoundaryNodes, gridAtlas.numCells);
			_boundaryCellBasedRHS.setFromTriplets(boundaryCellStiffnessTriplets[c].begin(), boundaryCellStiffnessTriplets[c].end());
			SparseMatrixParser(_boundaryCellBasedRHS, boundaryCellBasedStiffnessRHSMatrix[c], false);
		}
	}

	return 1;
}


int InitializeMassAndStiffness(std::vector<double> & deepMassCoefficients, std::vector<double> & deepStiffnessCoefficients,
	SparseMatrix<double, int> & boundaryBoundaryMassMatrix, SparseMatrix<double, int> & boundaryBoundaryStiffnessMatrix,
	SparseMatrix<double, int> & boundaryDeepMassMatrix, SparseMatrix<double, int> & boundaryDeepStiffnessMatrix,
	const HierarchicalSystem & hierarchy, const std::vector<std::vector<SquareMatrix<double, 2>>> & parameterMetric, const std::vector<AtlasChart> & atlasCharts,
	const BoundaryProlongationData & boundaryProlongation, bool computeCellBasedStiffness,
	const std::vector<Point3D<double>> & inputSignal, std::vector<double> & texelToCellCoeffs, SparseMatrix<double, int> boundaryCellBasedStiffnessRHSMatrix[3]) {

	//(2) Initialize mass and stiffness
	deepMassCoefficients.resize(10 * hierarchy.gridAtlases[0].numDeepTexels, 0);
	deepStiffnessCoefficients.resize(10 * hierarchy.gridAtlases[0].numDeepTexels, 0);


	SparseMatrix<double, int> fineBoundaryBoundaryMassMatrix;
	SparseMatrix<double, int> fineBoundaryBoundaryStiffnessMatrix;

	SparseMatrix<double, int> fineBoundaryCellStiffnessRHSMatrix[3];
	std::vector<Point3D<double>> fineBoundarySignal;
	
	if (computeCellBasedStiffness){
		const std::vector<int> & boundaryGlobalIndex = hierarchy.gridAtlases[0].boundaryGlobalIndex;
		int numBoundaryTexels = boundaryGlobalIndex.size();
		int numFineBoundaryNodes = boundaryProlongation.numFineBoundarNodes;
		std::vector<Point3D<double>> coarseBoundarySignal;
		coarseBoundarySignal.resize(numBoundaryTexels);
		for (int i = 0; i < numBoundaryTexels; i++)  coarseBoundarySignal[i] = inputSignal[boundaryGlobalIndex[i]];
		fineBoundarySignal.resize(numFineBoundaryNodes);
		boundaryProlongation.coarseBoundaryFineBoundaryProlongation.Multiply(&coarseBoundarySignal[0], &fineBoundarySignal[0]);
	}

	//clock_t m_begin = clock();
	if (!InitializeMassAndStiffness(parameterMetric, atlasCharts, hierarchy.gridAtlases[0], boundaryProlongation.fineBoundaryIndex, boundaryProlongation.numFineBoundarNodes, deepMassCoefficients, deepStiffnessCoefficients, fineBoundaryBoundaryMassMatrix, fineBoundaryBoundaryStiffnessMatrix, boundaryDeepMassMatrix, boundaryDeepStiffnessMatrix,
		computeCellBasedStiffness, inputSignal, fineBoundarySignal, texelToCellCoeffs, fineBoundaryCellStiffnessRHSMatrix)) {
		printf("ERROR: Unable to initialize fine mass and stiffness! \n");
		return 0;
	}

	SparseMatrix<double, int> temp = fineBoundaryBoundaryMassMatrix*boundaryProlongation.coarseBoundaryFineBoundaryProlongation;
	boundaryBoundaryMassMatrix = boundaryProlongation.fineBoundaryCoarseBoundaryRestriction * temp;

	temp = fineBoundaryBoundaryStiffnessMatrix*boundaryProlongation.coarseBoundaryFineBoundaryProlongation;
	boundaryBoundaryStiffnessMatrix = boundaryProlongation.fineBoundaryCoarseBoundaryRestriction * temp;

	if (1) {
		std::vector<double> in(boundaryBoundaryMassMatrix.Rows(), 1.0);
		std::vector<double> out(boundaryBoundaryMassMatrix.Rows(), 0.0);
		boundaryBoundaryMassMatrix.Multiply(GetPointer(in), GetPointer(out));
		for (int i = 0; i < out.size(); i++) {
			if (out[i] == 0.0) printf("WARNING: Zero row at index %d. \n", i);
		}
	}


	if (computeCellBasedStiffness) {
		for (int c = 0; c < 3; c++) {
			boundaryCellBasedStiffnessRHSMatrix[c] = boundaryProlongation.fineBoundaryCoarseBoundaryRestriction * fineBoundaryCellStiffnessRHSMatrix[c];
		}
	}

	return 1;
}

