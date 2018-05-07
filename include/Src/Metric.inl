#pragma once

enum {
	EMBEDDING_METRIC,
	UNIFORM_METRIC,
	UMBILIC_METRIC,
	EDGE_AWARE_METRIC,
};



int InitializeVectorFieldMetric(const std::vector<SquareMatrix<double, 2>> & embeddingMetric, const std::vector<Point2D<double>> & vf, bool normalizeArea, std::vector<SquareMatrix<double, 2>> & outputMetric) {
	int tCount = embeddingMetric.size();

	//Compute vf scale
	double vfScale = 0;
	for (int t = 0; t < tCount; t++) {
		vfScale += Point2D<double>::Dot(vf[t], embeddingMetric[t] * vf[t])* (sqrt(embeddingMetric[t].determinant()) / 2.0);
	}
	double vfNormalization = 1.0 / sqrt(vfScale);

	outputMetric.resize(embeddingMetric.size());
	double totalMass = 0;

	for (int t = 0; t < tCount; t++) {
		SquareMatrix<double, 2> g = embeddingMetric[t];
		SquareMatrix<double, 2> u;
		if (Point2D<double>::Dot(vf[t], embeddingMetric[t] * vf[t]) > 0) {
			Point2D<double> vectors[2];
			vectors[0] = vf[t];
			Point2D<double> temp = g*vectors[0];
			vectors[1] = Point2D<double>(-temp[1], temp[0]);
			for (int e = 0; e < 2; e++) vectors[e] /= sqrt(Point2D<double>::Dot(vectors[e], g*vectors[e]));

			double principalScale = sqrt(Point2D<double>::Dot(vf[t], embeddingMetric[t] * vf[t]))*vfNormalization;

			SquareMatrix<double, 2> B;
			B(0, 0) = vectors[0][0]; B(0, 1) = vectors[0][1];
			B(1, 0) = vectors[1][0]; B(1, 1) = vectors[1][1];


			SquareMatrix<double, 2> D;
			D(0, 0) = 0;
#if 0
			D(1, 1) = principalScale;
#else
			D(1, 1) = 1e-1;
#endif
			D(0, 1) = D(1, 0) = 0.0;

			u = g * B * D * B.inverse();
		}
		else {
			u *= 0;
		}
		outputMetric[t] = u + g*(1e-5);

		totalMass += sqrt(outputMetric[t].determinant()) / 2;
	}

	if (normalizeArea) {
		for (int t = 0; t < tCount; t++) {
			outputMetric[t] /= totalMass;
		}
	}
	return 1;
}

#define NORMALIZE_SURFACE_EMBEDDING 1

int InitializeCurvatureBasedMetric(const TexturedMesh & mesh, bool normalizeArea, std::vector<SquareMatrix<double, 2>> & outputMetric, std::vector<Point3D<double>> & vNormals, int metricMode) {
	outputMetric.resize(mesh.triangles.size());

#if NORMALIZE_SURFACE_EMBEDDING
	double inputMass = GetMeshArea(mesh);
	double edgeScaling = 1.0 / sqrt(inputMass);
#endif

	double totalMass = 0;

	for (int t = 0; t < mesh.triangles.size(); t++) {

		Point3D< double > vPos[3];
		for (int i = 0; i < 3; i++) vPos[i] = mesh.vertices[mesh.triangles[t][i]];

		SquareMatrix<double, 2> g;
#if NORMALIZE_SURFACE_EMBEDDING
		Point3D< double > dv[2] = { (vPos[1] - vPos[0])*edgeScaling, (vPos[2] - vPos[0])*edgeScaling };
#else
		Point3D< double > dv[2] = { (vPos[1] - vPos[0]), (vPos[2] - vPos[0]) };
#endif
		for (int k = 0; k < 2; k++)for (int l = 0; l < 2; l++)g(k, l) = Point3D<double>::Dot(dv[k], dv[l]);

		Point3D< double > vNormal[3];
		for (int i = 0; i < 3; i++) vNormal[i] = vNormals[mesh.triangles[t][i]];
		Point3D<double> dn[2] = { vNormal[1] - vNormal[0], vNormal[2] - vNormal[0] };
		SquareMatrix<double, 2> gg;
		for (int k = 0; k < 2; k++)for (int l = 0; l < 2; l++) gg(k, l) = Point3D<double>::Dot(dn[k], dv[l]);
		gg(0, 1) = gg(1, 0) = (gg(0, 1) + gg(1, 0)) / 2.0;

		SquareMatrix<double, 2> S = g.inverse()*gg;

		if (metricMode == UMBILIC_METRIC) {
			double a = 1.0;
			double b = -(S(0, 0) + S(1, 1));
			double c = S(0, 0)*S(1, 1) - S(0, 1)*S(1, 0);
			double discriminant = b*b - 4.0 *a*c;
			if (discriminant < 0) {
				printf("Unexpected negative discriminant! \n");
				return 0;
			}

			discriminant = sqrt(discriminant);
			double roots[2] = { (-b - discriminant) / (2.0*a), (-b + discriminant) / (2.0*a) };
			Point2D<double> vectors[2] = { Point2D<double>(1, 0), Point2D<double>(0, 1) };
			if (S(1, 0) != 0 || S(0, 1) != 0) {
				if (fabs(S(1, 0)) > fabs(S(0, 1))) {
					vectors[0] = Point2D<double>(-S(1, 0), S(0, 0) - roots[0]);
				}
				else {
					vectors[0] = Point2D<double>(-(S(1, 1) - roots[0]), S(0, 1));
				}
			}

			Point2D<double> temp = g*vectors[0];
			vectors[1] = Point2D<double>(-temp[1], temp[0]);

			for (int e = 0; e < 2; e++) vectors[e] /= sqrt(Point2D<double>::Dot(vectors[e], g*vectors[e]));

			SquareMatrix<double, 2> B;
			B(0, 0) = vectors[0][0]; B(0, 1) = vectors[0][1];
			B(1, 0) = vectors[1][0]; B(1, 1) = vectors[1][1];


			//roots[0] = fabs(roots[0]) * meshScale;
			//roots[1] = fabs(roots[1]) * meshScale;
			////printf("%f %f \n", roots[0], roots[1]);
			//double minEvalue = std::min<double>(roots[0], roots[1]);
			//roots[0] += (epsilon - minEvalue);
			//roots[1] += (epsilon - minEvalue);

			//Major curves
			roots[0] = fabs(roots[1] - roots[0]);
			roots[1] = 0;

			//Minor curves
			//roots[1] = fabs(roots[1] - roots[0]);
			//roots[0] = 0;

			SquareMatrix<double, 2> D;
			D(0, 0) = roots[0]; D(1, 1) = roots[1];
			D(0, 1) = D(1, 0) = 0.0;

			SquareMatrix<double, 2> u = g * B * D * B.inverse();

			//g = u + g*(1e-10);
			g = u + g*(1e-8); //For Camel LIC minor
		}
		else {
			double trace = S.trace();
			double det = S.determinant();
			double squared_curvatures = trace*trace - 2.0*det;
			if (squared_curvatures < 0) {
				printf("Unexpected negative sum of curvatures!\n");
				return 0;
			}
			g *= (1e-10 + squared_curvatures);
		}

		totalMass += sqrt(g.determinant()) / 2;

		outputMetric[t] = g;
	}

	printf("Here Total Mass %g \n", totalMass);

	if (normalizeArea) {
		for (int t = 0; t < mesh.triangles.size(); t++) {
			outputMetric[t] /= totalMass;
		}
	}
	return 1;
}

int InitializePrincipalCurvatureDirection(const TexturedMesh & mesh, std::vector<Point2D<double>> & principalDirection, std::vector<Point3D<double>> & vNormals, bool useSmallestCurvarture = true) {
	principalDirection.resize(mesh.triangles.size());

#if NORMALIZE_SURFACE_EMBEDDING
	double inputMass = GetMeshArea(mesh);
	double edgeScaling = 1.0 / sqrt(inputMass);
#endif
	double totalMass = 0;
	for (int t = 0; t < mesh.triangles.size(); t++) {

		Point3D< double > vPos[3];
		for (int i = 0; i < 3; i++) vPos[i] = mesh.vertices[mesh.triangles[t][i]];

		SquareMatrix<double, 2> g;
#if NORMALIZE_SURFACE_EMBEDDING
		Point3D< double > dv[2] = { (vPos[1] - vPos[0])*edgeScaling, (vPos[2] - vPos[0])*edgeScaling };
#else
		Point3D< double > dv[2] = { (vPos[1] - vPos[0]), (vPos[2] - vPos[0]) };
#endif
		for (int k = 0; k < 2; k++)for (int l = 0; l < 2; l++)g(k, l) = Point3D<double>::Dot(dv[k], dv[l]);

		Point3D< double > vNormal[3];
		for (int i = 0; i < 3; i++) vNormal[i] = vNormals[mesh.triangles[t][i]];
		Point3D<double> dn[2] = { vNormal[1] - vNormal[0], vNormal[2] - vNormal[0] };
		SquareMatrix<double, 2> gg;
		for (int k = 0; k < 2; k++)for (int l = 0; l < 2; l++) gg(k, l) = Point3D<double>::Dot(dn[k], dv[l]);
		gg(0, 1) = gg(1, 0) = (gg(0, 1) + gg(1, 0)) / 2.0;

		SquareMatrix<double, 2> S = g.inverse()*gg;

		double a = 1.0;
		double b = -(S(0, 0) + S(1, 1));
		double c = S(0, 0)*S(1, 1) - S(0, 1)*S(1, 0);
		double discriminant = b*b - 4.0 *a*c;
		if (discriminant < 0) {
			printf("Unexpected negative discriminant! \n");
			return 0;
		}

		discriminant = sqrt(discriminant);
		double root = (-b - discriminant) / (2.0*a);
		Point2D<double> vector =  Point2D<double>(1, 0);
		if (S(1, 0) != 0 || S(0, 1) != 0) {
			if (fabs(S(1, 0)) > fabs(S(0, 1))) {
				vector = Point2D<double>(-S(1, 0), S(0, 0) - root);
			}
			else {
				vector = Point2D<double>(-(S(1, 1) - root), S(0, 1));
			}
		}
		if (useSmallestCurvarture) {
			principalDirection[t] = vector;
		}
		else {
			Point2D<double> temp = g*vector;
			principalDirection[t] = Point2D<double>(-temp[1], temp[0]);
		}
		principalDirection[t] /= sqrt(Point2D<double>::Dot(principalDirection[t], g*principalDirection[t]));
	}
	return 1;
}


int InitializeEmbeddingMetric(const TexturedMesh & mesh, bool normalizeArea, std::vector<SquareMatrix<double, 2>> & embeddingMetric) {
	//double meshArea = GetMeshArea(mesh);
	//printf("Mesh area %f \n", meshArea);
	embeddingMetric.resize(mesh.triangles.size());

	double totalMass = 0;

	for (int t = 0; t < mesh.triangles.size(); t++) {

		Point3D< double > vPos[3];
		for (int i = 0; i < 3; i++) vPos[i] = mesh.vertices[mesh.triangles[t][i]];

		SquareMatrix<double, 2> g;
		Point3D< double > dv[2] = { vPos[1] - vPos[0], vPos[2] - vPos[0] };
		for (int k = 0; k < 2; k++)for (int l = 0; l < 2; l++)g(k, l) = Point3D<double>::Dot(dv[k], dv[l]);

		totalMass += sqrt(g.determinant()) / 2;

		embeddingMetric[t] = g;
	}

	if (normalizeArea) {
		for (int t = 0; t < mesh.triangles.size(); t++) {
			embeddingMetric[t] /= totalMass;
		}
	}

	return 1;
}

int InitializeUniformMetric(const TexturedMesh & mesh, bool normalizeArea, std::vector<SquareMatrix<double, 2>> & embeddingMetric) {
	//double meshArea = GetMeshArea(mesh);
	//printf("Mesh area %f \n", meshArea);
	embeddingMetric.resize(mesh.triangles.size());

	double totalMass = 0;

	for (int t = 0; t < mesh.triangles.size(); t++) {

		Point2D< double > tPos[3];
		for (int i = 0; i < 3; i++) tPos[i] = mesh.textureCoordinates[3 * t + i];

		SquareMatrix<double, 2> g;
		Point2D< double > dt[2] = { tPos[1] - tPos[0], tPos[2] - tPos[0] };
		for (int k = 0; k < 2; k++)for (int l = 0; l < 2; l++)g(k, l) = Point2D<double>::Dot(dt[k], dt[l]);

		totalMass += sqrt(g.determinant()) / 2;

		embeddingMetric[t] = g;
	}

	if (normalizeArea) {
		for (int t = 0; t < mesh.triangles.size(); t++) {
			embeddingMetric[t] /= totalMass;
		}
	}

	return 1;
}


//int InitializeParameterMetricAndDifferentials(const TexturedMesh & mesh, const std::vector<SquareMatrix<double, 2>> & embeddingMetric, std::vector<SquareMatrix<double, 2>> & parameterMetric, std::vector<SquareMatrix<double, 2>> & parameterToEmbeddingDifferential, std::vector<SquareMatrix<double, 2>> & parameterToEmbeddingRotatedDifferential) {
//
//	parameterMetric.resize(mesh.triangles.size());
//	parameterToEmbeddingDifferential.resize(mesh.triangles.size());
//	parameterToEmbeddingRotatedDifferential.resize(mesh.triangles.size());
//	SquareMatrix<double, 2> rot90;
//	rot90(0, 0) = rot90(1, 1) = 0; rot90(0, 1) = 1;	rot90(1, 0) = -1;
//
//	for (int t = 0; t < mesh.triangles.size(); t++) {
//
//		Point3D< double > vPos[3];
//		for (int i = 0; i < 3; i++) vPos[i] = mesh.vertices[mesh.triangles[t][i]];
//
//		SquareMatrix<double, 2> embedding_metric = embeddingMetric[t];
//
//		Point2D< double > tPos[3];
//		for (int i = 0; i < 3; i++) tPos[i] = mesh.textureCoordinates[3 * t + i];
//
//		Point2D< double > dp[2] = { tPos[1] - tPos[0], tPos[2] - tPos[0] };
//
//		//Parametric map
//		SquareMatrix<double, 2> parametric_map_differential;
//		parametric_map_differential(0, 0) = dp[0][0];
//		parametric_map_differential(0, 1) = dp[0][1];
//
//		parametric_map_differential(1, 0) = dp[1][0];
//		parametric_map_differential(1, 1) = dp[1][1];
//
//		SquareMatrix<double, 2> inverse_parametric_map_differential = parametric_map_differential.inverse();
//		SquareMatrix<double, 2> parameter_metric = inverse_parametric_map_differential.transpose() * embedding_metric * inverse_parametric_map_differential;
//
//		parameterToEmbeddingDifferential[t] = inverse_parametric_map_differential*parameter_metric.inverse();
//		parameterToEmbeddingRotatedDifferential[t] = inverse_parametric_map_differential*rot90 / sqrt(parameter_metric.determinant());
//		parameterMetric[t] = parameter_metric;
//	}
//
//	return 1;
//}

//int InitializeParameterMetric(const TexturedMesh & mesh, const std::vector<SquareMatrix<double, 2>> & embeddingMetric, std::vector<SquareMatrix<double, 2>> & parameterMetric) {
//
//	parameterMetric.resize(mesh.triangles.size());
//
//	for (int t = 0; t < mesh.triangles.size(); t++) {
//
//		Point3D< double > vPos[3];
//		for (int i = 0; i < 3; i++) vPos[i] = mesh.vertices[mesh.triangles[t][i]];
//
//		SquareMatrix<double, 2> embedding_metric = embeddingMetric[t];
//
//		Point2D< double > tPos[3];
//		for (int i = 0; i < 3; i++) tPos[i] = mesh.textureCoordinates[3 * t + i];
//
//		Point2D< double > dp[2] = { tPos[1] - tPos[0], tPos[2] - tPos[0] };
//
//		//Parametric map
//		SquareMatrix<double, 2> parametric_map_differential;
//		parametric_map_differential(0, 0) = dp[0][0];
//		parametric_map_differential(0, 1) = dp[0][1];
//
//		parametric_map_differential(1, 0) = dp[1][0];
//		parametric_map_differential(1, 1) = dp[1][1];
//
//		SquareMatrix<double, 2> inverse_parametric_map_differential = parametric_map_differential.inverse();
//		SquareMatrix<double, 2> parameter_metric = inverse_parametric_map_differential.transpose() * embedding_metric * inverse_parametric_map_differential;
//
//		parameterMetric[t] = parameter_metric;
//	}
//
//	return 1;
//}
//


//int InitializeMetric( AtlasMesh & atlasMesh, TexturedMesh & mesh, const int metricMode) {
//	std::vector<SquareMatrix<double, 2>> surfaceMetric;
//	std::vector<SquareMatrix<double, 2>> embeddingMetric;
//
//	InitializeEmbeddingMetric(mesh, true, embeddingMetric);
//
//	if (metricMode == EMBEDDING_METRIC) {
//		surfaceMetric = embeddingMetric;
//	}
//	else if (metricMode == UNIFORM_METRIC) {
//		InitializeUniformMetric(mesh, true, surfaceMetric);
//	}
//	else if (metricMode == UMBILIC_METRIC || metricMode == EDGE_AWARE_METRIC) {
//		UpdateNormals(mesh);
//		Eigen::SparseMatrix<double> meshMassMatrix;
//		Eigen::SparseMatrix<double> meshStiffnessMatrix;
//		InitializeMeshMatrices(mesh, meshMassMatrix, meshStiffnessMatrix);
//#if USE_PARDISO
//		Eigen::PardisoLDLT< Eigen::SparseMatrix< double > > meshSolver(meshMassMatrix + meshStiffnessMatrix*1e-4);
//#else
//		Eigen::SimplicialLDLT< Eigen::SparseMatrix< double > > meshSolver(meshMassMatrix + meshStiffnessMatrix*1e-4);
//#endif
//		SmoothSignal(meshMassMatrix, meshSolver, mesh.normals, true);
//		InitializeCurvatureBasedMetric(mesh, true, surfaceMetric, mesh.normals, metricMode);
//	}
//
//	InitializeParameterMetric(mesh, surfaceMetric, atlasMesh.metric);
//
//	return 1;
//}


int InitializeParameterMetric(const TexturedMesh & mesh, const std::vector<SquareMatrix<double, 2>> & embeddingMetric, const std::vector<AtlasChart> & atlasCharts, std::vector<std::vector<SquareMatrix<double, 2>>> & parameterMetric) {

	parameterMetric.resize(atlasCharts.size());
	for (int i = 0; i < atlasCharts.size(); i++) {
		parameterMetric[i].resize(atlasCharts[i].meshTriangleIndices.size());
		for (int k = 0; k < atlasCharts[i].meshTriangleIndices.size(); k++) {
			int t = atlasCharts[i].meshTriangleIndices[k];
			Point3D< double > vPos[3];
			for (int i = 0; i < 3; i++) vPos[i] = mesh.vertices[mesh.triangles[t][i]];

			SquareMatrix<double, 2> embedding_metric = embeddingMetric[t];

			Point2D< double > tPos[3];
			for (int i = 0; i < 3; i++) tPos[i] = mesh.textureCoordinates[3 * t + i];

			Point2D< double > dp[2] = { tPos[1] - tPos[0], tPos[2] - tPos[0] };

			//Parametric map
			SquareMatrix<double, 2> parametric_map_differential;
			parametric_map_differential(0, 0) = dp[0][0];
			parametric_map_differential(0, 1) = dp[0][1];

			parametric_map_differential(1, 0) = dp[1][0];
			parametric_map_differential(1, 1) = dp[1][1];

			SquareMatrix<double, 2> inverse_parametric_map_differential = parametric_map_differential.inverse();
			SquareMatrix<double, 2> parameter_metric = inverse_parametric_map_differential.transpose() * embedding_metric * inverse_parametric_map_differential;
			parameterMetric[i][k] = parameter_metric;
		}
	}
	return 1;
}

int InitializeMetric(TexturedMesh & mesh, const int metricMode, const std::vector<AtlasChart> & atlasCharts, std::vector<std::vector<SquareMatrix<double, 2>>> & parameterMetric) {
	std::vector<SquareMatrix<double, 2>> surfaceMetric;
	std::vector<SquareMatrix<double, 2>> embeddingMetric;

	InitializeEmbeddingMetric(mesh, true, embeddingMetric);

	if (metricMode == EMBEDDING_METRIC) {
		surfaceMetric = embeddingMetric;
	}
	else if (metricMode == UNIFORM_METRIC) {
		InitializeUniformMetric(mesh, true, surfaceMetric);
	}
	else{
		printf("Not recognized  metric! \n");
		return 0;
	}
	InitializeParameterMetric(mesh, surfaceMetric, atlasCharts, parameterMetric);
	return 1;
}

int InitializeAnisotropicMetric(TexturedMesh & mesh, const std::vector<AtlasChart> & atlasCharts, std::vector<std::vector<SquareMatrix<double, 2>>> & parameterMetric, std::vector<Point2D<double>> & vf,  bool useSmallestCurvarture) {
	std::vector<SquareMatrix<double, 2>> surfaceMetric;
	std::vector<SquareMatrix<double, 2>> embeddingMetric;

	InitializeEmbeddingMetric(mesh, true, embeddingMetric);

	if (vf.size() == 0) {
		UpdateNormals(mesh);
		Eigen::SparseMatrix<double> meshMassMatrix;
		Eigen::SparseMatrix<double> meshStiffnessMatrix;
		InitializeMeshMatrices(mesh, meshMassMatrix, meshStiffnessMatrix);
		Eigen::SimplicialLDLT< Eigen::SparseMatrix< double > > meshSolver(meshMassMatrix + meshStiffnessMatrix*1e-4);
		SmoothSignal(meshMassMatrix, meshSolver, mesh.normals, true);
		InitializePrincipalCurvatureDirection(mesh, vf, mesh.normals, useSmallestCurvarture);
	}
	InitializeVectorFieldMetric(embeddingMetric, vf, true, surfaceMetric);
	InitializeParameterMetric(mesh, surfaceMetric, atlasCharts, parameterMetric);
	return 1;
}