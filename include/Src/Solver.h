#pragma once

#define USE_CHOLMOD 0
#define USE_EIGEN_SIMPLICIAL 1
#define USE_EIGEN_PARDISO 0

#if USE_CHOLMOD
#include <Misha/LinearSolvers.h>

#define CHOLMOD_CHANNELS_IN_PARALLEL 1
#if CHOLMOD_CHANNELS_IN_PARALLEL
template <class Real>
class CholmodCholeskySolver3{
public:
	CholmodSolver<1> solver[3];
	std::vector<Real> out[3];
	std::vector<Real> in[3];
	void init(const SparseMatrix<double, int> & M){
#pragma omp parallel for
		for (int c = 0; c < 3; c++){
			solver[c]._init(M);
		}

		const int numVariables = M.Rows();
		for (int c = 0; c < 3; c++){
			out[c].resize(numVariables);
			in[c].resize(numVariables);
		}
	}

	void update(const SparseMatrix<double, int> & M){
#pragma omp parallel for
		for (int c = 0; c < 3; c++){
			solver[c]._update(M);
		}
	}
};

template <class Real, class DataType>
void solve(CholmodCholeskySolver3<Real> & chol, std::vector<DataType> & x0, const std::vector<DataType> & rhs){
	int numVariables = x0.size();
#pragma omp parallel for
	for (int c = 0; c < 3; c++){
		for (int n = 0; n < numVariables; n++) chol.in[c][n] = rhs[n][c];
		chol.solver[c].solve(&chol.in[c][0], &chol.out[c][0]);
	}
#pragma omp parallel for
	for (int n = 0; n < numVariables; n++) x0[n] = SetData(chol.out[0][n], chol.out[1][n], chol.out[2][n]);
}
#else 

#define CHOLMOD_CHANNELS_IN_BLOCK 1
#if CHOLMOD_CHANNELS_IN_BLOCK
template <class Real>
class CholmodCholeskySolver3{
public:
	CholmodSolver<3> solver;
	std::vector<Real> out;
	std::vector<Real> in;
	void init(const SparseMatrix<double, int> & M) {
		solver._init(M);
		const int numVariables = M.Rows();
		out.resize(3 * numVariables);
		in.resize(3 * numVariables);
	}

	void update(const SparseMatrix<double, int> & M) {
		solver._update(M);
	}
};

template <class Real, class DataType>
void solve(CholmodCholeskySolver3<Real> & chol, std::vector<DataType> & x0, const std::vector<DataType> & rhs){
	int numVariables = x0.size();
	for (int n = 0; n < numVariables; n++) for (int c = 0; c < 3; c++)chol.in[c* numVariables + n] = rhs[n][c];
	chol.solver.solve(&chol.in[0], &chol.out[0]);
	for (int n = 0; n < numVariables; n++) x0[n] = SetData(chol.out[0 * numVariables + n], chol.out[1 * numVariables + n], chol.out[2 * numVariables + n]);
}
#else
template <class Real>
class CholmodCholeskySolver3{
public:
	CholmodSolver<1> solver;
	std::vector<Real> out[3];
	std::vector<Real> in[3];
	void init(const SparseMatrix<double, int> & M) {
		solver._init(M);
		const int numVariables = M.Rows();
		for (int c = 0; c < 3; c++){
			out[c].resize(numVariables);
			in[c].resize(numVariables);
		}
	}

	void update(const SparseMatrix<double, int> & M) {
		solver._update(M);
	}
};

template <class Real, class DataType>
void solve(CholmodCholeskySolver3<Real> & chol, std::vector<DataType> & x0, const std::vector<DataType> & rhs) {
	int numVariables = x0.size();
	for (int c = 0; c < 3; c++) {
		for (int n = 0; n < numVariables; n++) chol.in[c][n] = rhs[n][c];
		chol.solver.solve(&chol.in[c][0], &chol.out[c][0]);
	}
#pragma omp parallel for
	for (int n = 0; n < numVariables; n++) x0[n] = SetData(chol.out[0][n], chol.out[1][n], chol.out[2][n]);
}
#endif
#endif

template <class Real>
class CholmodCholeskySolver1{
public:
	CholmodSolver<1> solver;
	std::vector<Real> out;
	std::vector<Real> in;
	void init(const SparseMatrix<double, int> & M) {
		solver._init(M);
		const int numVariables = M.Rows();
		out.resize(numVariables);
		in.resize(numVariables);
	}
	void update(const SparseMatrix<double, int> & M) {
		solver._update(M);
	}
};

template <class Real, class DataType>
void solve(CholmodCholeskySolver1<Real> & chol, std::vector<DataType> & x0, const std::vector<DataType> & rhs) {
	int numVariables = x0.size();
#pragma omp parallel for
	for (int n = 0; n < numVariables; n++) chol.in[n] = rhs[n];
	chol.solver.solve(&chol.in[0], &chol.out[0]);
#pragma omp parallel for
	for (int n = 0; n < numVariables; n++) x0[n] = chol.out[n];
}
#endif

#include <Eigen/Sparse>
#include <Eigen/Dense>

template <class Real>
class EigenCholeskySolver3{
public:
	typedef typename Eigen::SimplicialLDLT<Eigen::SparseMatrix<Real>> Solver;
	//typedef typename Eigen::SimplicialLLT<Eigen::SparseMatrix<Real>> Solver;
	Solver * solver;
	EigenVector x0_vectors[3];
	EigenVector rhs_vectors[3];
	EigenVector solution_vectors[3];

	void init(const SparseMatrix<double, int> & _M){
		Eigen::SparseMatrix< Real > M;
		SparseMatrixParser(_M, M);

		solver = new Solver();
		solver->analyzePattern(M);

		Eigen::ComputationInfo info = solver->info();
		if (info == Eigen::Success) {
		}
		else if (info == Eigen::NumericalIssue){
			printf("FAILED : Numerical issue! \n");
		}
		else if (info == Eigen::NoConvergence) {
			printf("FAILED : No convergence! \n");
		}
		else if (info == Eigen::InvalidInput) {
			printf("FAILED : Invalid input! \n");
		}
		else {
			printf("FAILED : Undetermined cause! \n");
		}

		const int numVariables = M.rows();
		for (int c = 0; c < 3; c++) {
			x0_vectors[c].resize(numVariables);
			rhs_vectors[c].resize(numVariables);
			solution_vectors[c].resize(numVariables);
		}
	}
	void update(const SparseMatrix<double, int> & _M) {
		Eigen::SparseMatrix< Real > M;
		SparseMatrixParser(_M, M);
		solver->factorize(M);
	}
};

template <class Real,class DataType>
void solve(EigenCholeskySolver3<Real> & chol, std::vector<DataType> & x0, const std::vector<DataType> & rhs){
	int numVariables = x0.size();
#pragma omp parallel for
	for (int n = 0; n < numVariables; n++) for (int c = 0; c < 3; c++) {
		chol.x0_vectors[c][n] = x0[n][c];
		chol.rhs_vectors[c][n] = rhs[n][c];
	}
#pragma omp parallel for
	for (int c = 0; c < 3; c++) chol.solution_vectors[c] = chol.solver->solve(chol.rhs_vectors[c]);
#pragma omp parallel for
	for (int n = 0; n < numVariables; n++) x0[n]= SetData(chol.solution_vectors[0][n], chol.solution_vectors[1][n], chol.solution_vectors[2][n]);
}

template <class Real>
class EigenCholeskySolver1{
public:
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<Real>> * solver;
	EigenVector x0_vector;
	EigenVector rhs_vector;
	EigenVector solution_vector;

	void init(const SparseMatrix<double, int> & _M) {
		Eigen::SparseMatrix< Real > M;
		SparseMatrixParser(_M, M);

		solver = new Eigen::SimplicialLDLT<Eigen::SparseMatrix<Real>>();
		solver->analyzePattern(M);
		Eigen::ComputationInfo info = solver->info();
		if (info == Eigen::Success) {
		}
		else if (info == Eigen::NumericalIssue) {
			printf("FAILED : Numerical issue! \n");
		}
		else if (info == Eigen::NoConvergence) {
			printf("FAILED : No convergence! \n");
		}
		else if (info == Eigen::InvalidInput) {
			printf("FAILED : Invalid input! \n");
		}
		else {
			printf("FAILED : Undetermined cause! \n");
		}

		const int numVariables = M.rows();
		x0_vector.resize(numVariables);
		rhs_vector.resize(numVariables);
		solution_vector.resize(numVariables);
	}

	void update(const SparseMatrix<double, int> & _M){
		Eigen::SparseMatrix< Real > M;
		SparseMatrixParser(_M, M);
		solver->factorize(M);
	}
};

template <class Real, class DataType>
void solve(EigenCholeskySolver1<Real> & chol, std::vector<DataType> & x0, const std::vector<DataType> & rhs){
	int numVariables = x0.size();
#pragma omp parallel for
	for (int n = 0; n < numVariables; n++) {
		chol.x0_vector[n] = x0[n];
		chol.rhs_vector[n] = rhs[n];
	}
	chol.solution_vector = chol.solver->solve(chol.rhs_vector);
#pragma omp parallel for
	for (int n = 0; n < numVariables; n++) x0[n] = chol.solution_vector[n];
}

#if USE_EIGEN_PARDISO
#pragma comment( lib , "mkl_intel_lp64.lib")
#pragma comment( lib , "mkl_intel_thread.lib")
#pragma comment( lib , "mkl_core.lib")
#pragma comment( lib , "libiomp5md.lib")

#include <Eigen/PardisoSupport>
//Pardiso Solver

template <class Real>
class EigenPardisoSolver3 {
public:
	Eigen::PardisoLDLT<Eigen::SparseMatrix<Real>>  * solver;
	EigenVector x0_vectors[3];
	EigenVector rhs_vectors[3];
	EigenVector solution_vectors[3];

	void init(const SparseMatrix<double, int> & _M) {
		Eigen::SparseMatrix< Real > M;
		SparseMatrixParser(_M, M);

		solver = new Eigen::PardisoLDLT<Eigen::SparseMatrix<Real>>();
		solver->analyzePattern(M);
		Eigen::ComputationInfo info = solver->info();
		if (info == Eigen::Success) {
		}
		else if (info == Eigen::NumericalIssue) {
			printf("FAILED : Numerical issue! \n");
		}
		else if (info == Eigen::NoConvergence) {
			printf("FAILED : No convergence! \n");
		}
		else if (info == Eigen::InvalidInput) {
			printf("FAILED : Invalid input! \n");
		}
		else {
			printf("FAILED : Undetermined cause! \n");
		}

		const int numVariables = M.rows();
		for (int c = 0; c < 3; c++) {
			x0_vectors[c].resize(numVariables);
			rhs_vectors[c].resize(numVariables);
			solution_vectors[c].resize(numVariables);
		}
	}
	void update(const SparseMatrix<double, int> & _M) {
		Eigen::SparseMatrix< Real > M;
		SparseMatrixParser(_M, M);
		solver->factorize(M);
	}
};

template <class Real, class DataType>
void solve(EigenPardisoSolver3<Real> & chol, std::vector<DataType> & x0, const std::vector<DataType> & rhs) {
	int numVariables = x0.size();
#pragma omp parallel for
	for (int n = 0; n < numVariables; n++) for (int c = 0; c < 3; c++) {
		chol.x0_vectors[c][n] = x0[n][c];
		chol.rhs_vectors[c][n] = rhs[n][c];
	}
	for (int c = 0; c < 3; c++) chol.solution_vectors[c] = chol.solver->solve(chol.rhs_vectors[c]);
#pragma omp parallel for
	for (int n = 0; n < numVariables; n++) x0[n] = SetData(chol.solution_vectors[0][n], chol.solution_vectors[1][n], chol.solution_vectors[2][n]);
}

template <class Real>
class EigenPardisoSolver1{
public:
	Eigen::PardisoLDLT<Eigen::SparseMatrix<Real>> * solver;
	EigenVector x0_vector;
	EigenVector rhs_vector;
	EigenVector solution_vector;

	void init(const SparseMatrix<double, int> & _M){
		Eigen::SparseMatrix< Real > M;
		SparseMatrixParser(_M, M);

		solver = new Eigen::PardisoLDLT<Eigen::SparseMatrix<Real>>();
		solver->analyzePattern(M);
		Eigen::ComputationInfo info = solver->info();
		if (info == Eigen::Success) {
		}
		else if (info == Eigen::NumericalIssue) {
			printf("FAILED : Numerical issue! \n");
		}
		else if (info == Eigen::NoConvergence) {
			printf("FAILED : No convergence! \n");
		}
		else if (info == Eigen::InvalidInput) {
			printf("FAILED : Invalid input! \n");
		}
		else {
			printf("FAILED : Undetermined cause! \n");
		}

		const int numVariables = M.rows();
		x0_vector.resize(numVariables);
		rhs_vector.resize(numVariables);
		solution_vector.resize(numVariables);
	}
	void update(const SparseMatrix<double, int> & _M) {
		Eigen::SparseMatrix< Real > M;
		SparseMatrixParser(_M, M);
		solver->factorize(M);
	}
};

template <class Real, class DataType>
void solve(EigenPardisoSolver1<Real> & chol, std::vector<DataType> & x0, const std::vector<DataType> & rhs) {
	int numVariables = x0.size();
#pragma omp parallel for
	for (int n = 0; n < numVariables; n++) {
		chol.x0_vector[n] = x0[n];
		chol.rhs_vector[n] = rhs[n];
	}
	chol.solution_vector = chol.solver->solve(chol.rhs_vector);
#pragma omp parallel for
	for (int n = 0; n < numVariables; n++) x0[n] = chol.solution_vector[n];
}
#endif
