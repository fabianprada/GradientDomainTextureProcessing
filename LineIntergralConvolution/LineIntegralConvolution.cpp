#define USE_FLOAT 0
#include <Src/VectorPrecisionType.inl>

#include <Misha/CmdLineParser.h> 
#include <Src/Basis.h>
#include <Misha/FEM.h>
#include <Src/ChartDecomposition.h>
#include <Src/HSV.h>
#include <Src/Solver.h>
#include <Src/Hierarchy.h>
#include <Src/VectorFieldIntergration.inl>
#include <Src/MassAndStiffness.h>
#include <Src/Padding.h>
#include <Src/TexturedMeshVisualization.h>
#include <Src/RandomJitter.h>

cmdLineParameter<  char* > Mesh("mesh");
cmdLineParameter<  char* > OutputTexture("outTexture");
cmdLineReadable MajorCurvature("major");
cmdLineReadable MinorCurvature("minor");
cmdLineParameter<  char* > VectorField("vf");
cmdLineParameter<  int   > Width("width", 2048);
cmdLineParameter<  int   > Height("height", 2048);
cmdLineParameter<  float > ConvolutionWeight("licStiffness", 1e-4);
cmdLineParameter<  float > StiffnessWeight("stiffnessWeight", 1e-4);
cmdLineParameter<  float > SharpenningWeight("sharpnessWeight", 100);
cmdLineParameter<  int   > Levels("levels", 4);


cmdLineParameter<  char* > CameraConfig("camera");
cmdLineParameter<  int   > Threads("threads", omp_get_num_procs());
cmdLineParameter<  int   > DisplayMode("display", TWO_REGION_DISPLAY);


cmdLineParameter<  int   > MultigridBlockHeight("mBlockH", 16);
cmdLineParameter<  int   > MultigridBlockWidth("mBlockW", 128);
cmdLineParameter<  int   > MultigridPaddedHeight("mPadH", 0);
cmdLineParameter<  int   > MultigridPaddedWidth("mPadW", 2);

cmdLineReadable RandomJitter("jitter");
cmdLineReadable Verbose("verbose");
cmdLineReadable DetailVerbose("detail");
cmdLineReadable UseDirectSolver("useDirectSolver");
cmdLineReadable* params[] =
{
	&Mesh,&OutputTexture,&MajorCurvature, &MinorCurvature, &VectorField, &Width,&Height,&ConvolutionWeight,&StiffnessWeight, &SharpenningWeight, &CameraConfig, &Levels,&UseDirectSolver,&Threads,&DisplayMode,&MultigridBlockHeight,&MultigridBlockWidth,&MultigridPaddedHeight,&MultigridPaddedWidth,&Verbose,
	&DetailVerbose,&RandomJitter,
	NULL
};

void ShowUsage(const char* ex)
{
	printf("Usage %s:\n", ex);
	printf("\t --%s <input mesh>\n", Mesh.name);
	printf("\t --%s <output texture>\n", OutputTexture.name);
	printf("\t --%s <use major curvature field> (default)\n", MajorCurvature.name);
	printf("\t --%s <use minor curvature field>\n", MinorCurvature.name);
	printf("\t --%s <vector field file>\n", VectorField.name);
	printf("\t --%s <texture width> [%d]\n", Width.name, Width.value);
	printf("\t --%s <texture height> [%d]\n", Height.name, Height.value);
	printf("\t --%s <LIC stiffness weight> [%f]\n", ConvolutionWeight.name, ConvolutionWeight.value);
	printf("\t --%s <stiffness weight> [%f]\n", StiffnessWeight.name, StiffnessWeight.value);
	printf("\t --%s <sharpness weight> [%f]\n", SharpenningWeight.name, SharpenningWeight.value);
	printf("\t --%s <camera configuration >\n", CameraConfig.name);
	printf("\t --%s <hierarchy levels> [%d]\n", Levels.name, Levels.value);
	printf("\t --%s <threads> [%d]\n", Threads.name, Threads.value);
	printf("\t --%s <disable update>\n", UseDirectSolver.name);
	printf("\t --%s <verbose>\n", Verbose.name);
	printf("\t --%s <detail verbose>\n", DetailVerbose.name);
	printf("\t --%s <display mode> [%d]\n", DisplayMode.name, DisplayMode.value);
	printf("\t \t [%d] One Region \n", ONE_REGION_DISPLAY);
	printf("\t \t [%d] Two Region \n", TWO_REGION_DISPLAY);

	printf("\t --%s <multigrid block width> [%d]\n", MultigridBlockWidth.name, MultigridBlockWidth.value);
	printf("\t --%s <multigrid block height> [%d]\n", MultigridBlockHeight.name, MultigridBlockHeight.value);
	printf("\t --%s <multigrid padded width> [%d]\n", MultigridPaddedWidth.name, MultigridPaddedWidth.value);
	printf("\t --%s <multigrid padded height> [%d]\n", MultigridPaddedHeight.name, MultigridPaddedHeight.value);

	printf("\t --%s <random jitter> \n", RandomJitter.name);

}

template<class Real>
class LineConvolution {
public:
	static Real sharpenningWeight;
	static Real stiffnessWeight;
	static Real convolutionWeight;

	static TexturedMesh mesh;
	static int textureWidth;
	static int textureHeight;
	static int levels;
	static std::vector<Point3D<float>> textureNodePositions;

	static Padding padding;

	static HierarchicalSystem hierarchy;
	static std::vector<CellIndex> cellIndices;

	static std::vector<TextureNodeInfo> textureNodes;
	static Image<int> nodeIndex;

	static SparseMatrix<double, int> anisotropicMass;
	static SparseMatrix<double, int> anisotropicStiffness;

	static SparseMatrix<double, int> mass;
	static SparseMatrix<double, int> stiffness;


	static SparseMatrix< double, int > lineConvolutionMatrix;
	static SparseMatrix< double, int > sharpenningMatrix;

	static int impulseTexel;

	static std::vector<AtlasChart> atlasCharts;
	static std::vector<std::vector<SquareMatrix<double, 2>>> parameterMetric;

	static double lineConvolutionRange;
	static double sharpenningRange;

	static std::vector<VECTOR_TYPE>	randSignal;

	static std::vector<VECTOR_TYPE> exactLineConvolutionSolution;
	static std::vector<VECTOR_TYPE> exactSharpenningSolution;
	static std::vector<VECTOR_TYPE> mass_x0;
	static std::vector<VECTOR_TYPE> stiffness_x0;

	//Impulse Smoothing
	static std::vector<MultigridLevelCoefficients<Real>> multigridLineConvolutionCoefficients;
	static std::vector<MultigridLevelVariables<VECTOR_TYPE>> multigridLineConvolutionVariables;

	//Geodesic Distance
	static std::vector<MultigridLevelCoefficients<Real>> multigridSharpenningCoefficients;
	static std::vector<MultigridLevelVariables<VECTOR_TYPE>> multigridSharpenningVariables;

#if USE_CHOLMOD
	typedef  std::vector<CholmodCholeskySolver3<Real>> BoundarySolverType;
	typedef  CholmodCholeskySolver3<Real>  CoarseSolverType;
	typedef  CholmodCholeskySolver3<Real> DirectSolverType;
#elif USE_EIGEN_SIMPLICIAL
	typedef  std::vector<EigenCholeskySolver3<Real>> BoundarySolverType;
	typedef  EigenCholeskySolver3<Real>  CoarseSolverType;
	typedef  EigenCholeskySolver3<Real> DirectSolverType;
#else
	typedef  std::vector<EigenPardisoSolver3<Real>> BoundarySolverType;
	typedef  EigenPardisoSolver3<Real> CoarseSolverType;
	typedef  EigenPardisoSolver3<Real> DirectSolverType;
#endif

	static BoundarySolverType	boundaryLineConvolutionSolver;
	static BoundarySolverType	boundarySharpenningSolver;

	static CoarseSolverType coarseLineConvolutionSolver;
	static CoarseSolverType coarseSharpenningSolver;

	static DirectSolverType fineLineConvolutionSolver;
	static DirectSolverType fineSharpenningSolver;

	static std::vector<MultigridLevelIndices<Real>> multigridIndices;

	static SparseMatrix<Real, int> coarseBoundaryFineBoundaryProlongation;
	static SparseMatrix<Real, int> fineBoundaryCoarseBoundaryRestriction;
	static std::vector<Real> coarseBoundaryValues;
	static std::vector<Real> coarseBoundaryRHS;
	static std::vector<Real> fineBoundaryValues;
	static std::vector<Real> fineBoundaryRHS;

	//Anisotropic Linear Operators
	static std::vector<double> anisoDeepMassCoefficients;
	static std::vector<double> anisoDeepStiffnessCoefficients;
	static SparseMatrix<double, int> anisoBoundaryBoundaryMassMatrix;
	static SparseMatrix<double, int> anisoBoundaryBoundaryStiffnessMatrix;
	static SparseMatrix<double, int> anisoBoundaryDeepMassMatrix;
	static SparseMatrix<double, int> anisoBoundaryDeepStiffnessMatrix;

	//Isotropic Linear Operators
	static std::vector<double> deepMassCoefficients;
	static std::vector<double> deepStiffnessCoefficients;
	static SparseMatrix<double, int> boundaryBoundaryMassMatrix;
	static SparseMatrix<double, int> boundaryBoundaryStiffnessMatrix;
	static SparseMatrix<double, int> boundaryDeepMassMatrix;
	static SparseMatrix<double, int> boundaryDeepStiffnessMatrix;

	//Samples
	static std::vector<HexagonalSampleInfo<Real>> hexSamples;
	static std::vector<std::vector<CellLineSampleInfo<Real>>> cellLineSamples;
	static std::vector<InteriorCellLine> interiorCellLines;
	static std::vector<std::pair<int, int>> interiorCellLineIndex;

	static unsigned char * outputBuffer;

	//Visulization
	static TexturedMeshVisualization visualization;

	static void UpdateOutputBuffer(const std::vector<VECTOR_TYPE> & solution);

	static void StiffnessWeightCallBack(Visualization* v, const char* prompt);
	static void SharpenningWeightCallBack(Visualization* v, const char* prompt);
	static void ConvolutionWeightCallBack(Visualization* v, const char* prompt);

	static bool stopUpdate;

	static void StopUpdateCallBack(Visualization* v, const char* prompt);
	static void ExportTextureCallBack(Visualization* v, const char* prompt);

	static int Init();
	static void InitializeVisualization(const int width, const int height);
	static void ComputeExactSolution(bool verbose = false);
	static int UpdateSolution(bool verbose = false, bool detailVerbose = false);
	static int InitializeSystem(const int width, const int height);

	static void Display(void) { visualization.Display(); }
	static void MouseFunc(int button, int state, int x, int y);
	static void MotionFunc(int x, int y);
	static void Reshape(int w, int h) { visualization.Reshape(w, h); }
	static void KeyboardFunc(unsigned char key, int x, int y) { visualization.KeyboardFunc(key, x, y); }
	static void Idle();
};

template<class Real> Real														LineConvolution<Real>::sharpenningWeight;
template<class Real> Real														LineConvolution<Real>::stiffnessWeight;
template<class Real> Real														LineConvolution<Real>::convolutionWeight;

template<class Real> TexturedMesh												LineConvolution<Real>::mesh;
template<class Real> int														LineConvolution<Real>::textureWidth;
template<class Real> int														LineConvolution<Real>::textureHeight;

template<class Real> TexturedMeshVisualization									LineConvolution<Real>::visualization;

template<class Real> std::vector<AtlasChart>									LineConvolution<Real>::atlasCharts;
template<class Real> std::vector<std::vector<SquareMatrix<double, 2>>>			LineConvolution<Real>::parameterMetric;

template<class Real> Padding													LineConvolution<Real>::padding;

template<class Real> SparseMatrix<double, int>									LineConvolution<Real>::anisotropicMass;
template<class Real> SparseMatrix<double, int>									LineConvolution<Real>::anisotropicStiffness;


template<class Real> SparseMatrix<double, int>									LineConvolution<Real>::lineConvolutionMatrix;
template<class Real> SparseMatrix<double, int>									LineConvolution<Real>::sharpenningMatrix;

template<class Real> SparseMatrix<double, int>									LineConvolution<Real>::mass;
template<class Real> SparseMatrix<double, int>									LineConvolution<Real>::stiffness;


template<class Real> std::vector<TextureNodeInfo>								LineConvolution<Real>::textureNodes;
template<class Real> Image<int>													LineConvolution<Real>::nodeIndex;
template<class Real> std::vector<CellIndex>										LineConvolution<Real>::cellIndices;

template<class Real> int														LineConvolution<Real>::levels;
template<class Real> HierarchicalSystem											LineConvolution<Real>::hierarchy;

template<class Real> unsigned char *											LineConvolution<Real>::outputBuffer;
template<class Real> std::vector<MultigridLevelIndices<Real>>					LineConvolution<Real>::multigridIndices;

//Impulse Smoothing
template<class Real> std::vector<MultigridLevelCoefficients<Real>>						LineConvolution<Real>::multigridLineConvolutionCoefficients;
template<class Real> std::vector<MultigridLevelVariables<VECTOR_TYPE>>					LineConvolution<Real>::multigridLineConvolutionVariables;
template<class Real> typename LineConvolution<Real>::CoarseSolverType						LineConvolution<Real>::coarseLineConvolutionSolver;

//Geodesic Distance
template<class Real> std::vector<MultigridLevelCoefficients<Real>>						LineConvolution<Real>::multigridSharpenningCoefficients;
template<class Real> std::vector<MultigridLevelVariables<VECTOR_TYPE>>					LineConvolution<Real>::multigridSharpenningVariables;
template<class Real> typename LineConvolution<Real>::CoarseSolverType						LineConvolution<Real>::coarseSharpenningSolver;

template<class Real> typename LineConvolution<Real>::DirectSolverType						LineConvolution<Real>::fineLineConvolutionSolver;
template<class Real> typename LineConvolution<Real>::DirectSolverType						LineConvolution<Real>::fineSharpenningSolver;

template<class Real>  typename LineConvolution<Real>::BoundarySolverType					LineConvolution<Real>::boundaryLineConvolutionSolver;
template<class Real>  typename LineConvolution<Real>::BoundarySolverType					LineConvolution<Real>::boundarySharpenningSolver;

template<class Real> std::vector<VECTOR_TYPE>											LineConvolution<Real>::randSignal;

template<class Real> std::vector<VECTOR_TYPE>											LineConvolution<Real>::exactLineConvolutionSolution;
template<class Real> std::vector<VECTOR_TYPE>											LineConvolution<Real>::exactSharpenningSolution;

template<class Real> std::vector<VECTOR_TYPE>											LineConvolution<Real>::mass_x0;
template<class Real> std::vector<VECTOR_TYPE>											LineConvolution<Real>::stiffness_x0;

static std::vector<VECTOR_TYPE> mass_x0;
static std::vector<VECTOR_TYPE> stiffness_x0;
//Samples
template<class Real> std::vector<HexagonalSampleInfo<Real>>							LineConvolution<Real>::hexSamples;
template<class Real> std::vector<std::vector<CellLineSampleInfo<Real>>>				LineConvolution<Real>::cellLineSamples;
template<class Real> std::vector<InteriorCellLine>									LineConvolution<Real>::interiorCellLines;
template<class Real> std::vector<std::pair<int, int>>								LineConvolution<Real>::interiorCellLineIndex;

template<class Real> int															LineConvolution<Real>::impulseTexel = -1;
template<class Real> std::vector<Point3D<float>>									LineConvolution<Real>::textureNodePositions;

template<class Real> double															LineConvolution<Real>::lineConvolutionRange;
template<class Real> double															LineConvolution<Real>::sharpenningRange;

template<class Real> SparseMatrix<Real, int>										LineConvolution<Real>::coarseBoundaryFineBoundaryProlongation;
template<class Real> SparseMatrix<Real, int>										LineConvolution<Real>::fineBoundaryCoarseBoundaryRestriction;

template<class Real> std::vector<Real>												LineConvolution<Real>::coarseBoundaryValues;
template<class Real> std::vector<Real>												LineConvolution<Real>::coarseBoundaryRHS;
template<class Real> std::vector<Real>												LineConvolution<Real>::fineBoundaryValues;
template<class Real> std::vector<Real>												LineConvolution<Real>::fineBoundaryRHS;

template<class Real> std::vector<double>											LineConvolution<Real>::anisoDeepMassCoefficients;
template<class Real> std::vector<double>											LineConvolution<Real>::anisoDeepStiffnessCoefficients;
template<class Real> SparseMatrix<double, int>										LineConvolution<Real>::anisoBoundaryBoundaryMassMatrix;
template<class Real> SparseMatrix<double, int>										LineConvolution<Real>::anisoBoundaryBoundaryStiffnessMatrix;
template<class Real> SparseMatrix<double, int>										LineConvolution<Real>::anisoBoundaryDeepMassMatrix;
template<class Real> SparseMatrix<double, int>										LineConvolution<Real>::anisoBoundaryDeepStiffnessMatrix;

template<class Real> std::vector<double>											LineConvolution<Real>::deepMassCoefficients;
template<class Real> std::vector<double>											LineConvolution<Real>::deepStiffnessCoefficients;
template<class Real> SparseMatrix<double, int>										LineConvolution<Real>::boundaryBoundaryMassMatrix;
template<class Real> SparseMatrix<double, int>										LineConvolution<Real>::boundaryBoundaryStiffnessMatrix;
template<class Real> SparseMatrix<double, int>										LineConvolution<Real>::boundaryDeepMassMatrix;
template<class Real> SparseMatrix<double, int>										LineConvolution<Real>::boundaryDeepStiffnessMatrix;

template<class Real> bool															LineConvolution<Real>::stopUpdate = false;
template<class Real>
void LineConvolution<Real>::ComputeExactSolution(bool verbose){
	clock_t begin;

	//(1) Line Convolution	
	if (verbose) begin = clock();
	solve(fineLineConvolutionSolver, exactLineConvolutionSolution, multigridLineConvolutionVariables[0].rhs);
	if (verbose) printf("Line convolution %.4f \n", double(clock() - begin) / CLOCKS_PER_SEC);

	//(2) Compute sharpenning RHS
	mass_x0.resize(textureNodes.size());
	MultiplyBySystemMatrix_NoReciprocals(deepMassCoefficients, boundaryDeepMassMatrix, boundaryBoundaryMassMatrix, hierarchy.gridAtlases[0].boundaryGlobalIndex, hierarchy.gridAtlases[0].rasterLines, exactLineConvolutionSolution, mass_x0);

	stiffness_x0.resize(textureNodes.size());
	MultiplyBySystemMatrix_NoReciprocals(deepStiffnessCoefficients, boundaryDeepStiffnessMatrix, boundaryBoundaryStiffnessMatrix, hierarchy.gridAtlases[0].boundaryGlobalIndex, hierarchy.gridAtlases[0].rasterLines, exactLineConvolutionSolution, stiffness_x0);

	for (int i = 0; i < textureNodes.size(); i++) multigridSharpenningVariables[0].rhs[i] = mass_x0[i] + stiffnessWeight  * sharpenningWeight * stiffness_x0[i];

	//(3) Sharpenning
	if (verbose) begin = clock();
	solve(fineSharpenningSolver, exactSharpenningSolution, multigridSharpenningVariables[0].rhs);
	if (verbose) printf("Sharpenning %.4f \n", double(clock() - begin) / CLOCKS_PER_SEC);
}

template<class Real>
void LineConvolution<Real>::UpdateOutputBuffer(const std::vector<VECTOR_TYPE> & solution){

#pragma omp parallel for
	for (int i = 0; i < textureNodes.size(); i++){
		int ci = textureNodes[i].ci;
		int cj = textureNodes[i].cj;
		int offset = 3 * (textureWidth*cj + ci);
		outputBuffer[offset + 0] = (unsigned char)(std::min<double>(std::max<double>(0, solution[i][0]),1.0)*255.0);
		outputBuffer[offset + 1] = (unsigned char)(std::min<double>(std::max<double>(0, solution[i][1]), 1.0)*255.0);
		outputBuffer[offset + 2] = (unsigned char)(std::min<double>(std::max<double>(0, solution[i][2]), 1.0)*255.0);
	}

	glBindTexture(GL_TEXTURE_2D, visualization.textureBuffer);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, textureWidth, textureHeight, 0, GL_RGB, GL_UNSIGNED_BYTE, (GLvoid*)&outputBuffer[0]);
	glBindTexture(GL_TEXTURE_2D, 0);

	glutPostRedisplay();
}

template<class Real>
void LineConvolution<Real>::Idle() {
	
	if (!stopUpdate) {
		if (!UpdateSolution()) {
			printf("Updated solution failed! \n");
		}
		UpdateOutputBuffer(multigridSharpenningVariables[0].x);
	}
}

template<class Real>
void LineConvolution<Real>::MouseFunc(int button, int state, int x, int y) {

	visualization.newX = x; visualization.newY = y;
	visualization.rotating = visualization.scaling = visualization.panning = false;

	if (button == GLUT_LEFT_BUTTON) {
		if (glutGetModifiers() & GLUT_ACTIVE_CTRL) visualization.panning = true;
		else                                        visualization.rotating = true;
	}
	else if (button == GLUT_RIGHT_BUTTON) visualization.scaling = true;
}

template<class Real>
void LineConvolution<Real>::MotionFunc(int x, int y) {

	if (!visualization.showMesh) {
		visualization.oldX = visualization.newX, visualization.oldY = visualization.newY, visualization.newX = x, visualization.newY = y;

		int imageSize = std::min< int >(visualization.screenWidth, visualization.screenHeight);
		if (visualization.panning) visualization.xForm.offset[0] -= (visualization.newX - visualization.oldX) / visualization.imageToScreenScale(), visualization.xForm.offset[1] += (visualization.newY - visualization.oldY) / visualization.imageToScreenScale();
		else
		{
			float dz = float(pow(1.1, double(visualization.newY - visualization.oldY) / 8));
			visualization.xForm.zoom *= dz;
		}

	}
	else {
		visualization.oldX = visualization.newX, visualization.oldY = visualization.newY, visualization.newX = x, visualization.newY = y;
		int screenSize = std::min< int >(visualization.screenWidth, visualization.screenHeight);
		float rel_x = (visualization.newX - visualization.oldX) / (float)screenSize * 2;
		float rel_y = (visualization.newY - visualization.oldY) / (float)screenSize * 2;

		float pRight = rel_x * visualization.zoom, pUp = -rel_y * visualization.zoom;
		float pForward = rel_y * visualization.zoom;
		float rRight = -rel_y, rUp = -rel_x;

		if (visualization.rotating) visualization.camera.rotateUp(visualization.center, rUp), visualization.camera.rotateRight(visualization.center, rRight);
		else if (visualization.scaling) visualization.camera.moveForward(pForward);
		else if (visualization.panning) visualization.camera.moveRight(pRight), visualization.camera.moveUp(pUp);
	}
	glutPostRedisplay();
}


template<class Real>
void LineConvolution<Real>::StiffnessWeightCallBack(Visualization* v, const char* prompt) {
	stopUpdate = false;
	for (int i = 0; i < multigridLineConvolutionVariables[0].x.size(); i++)multigridLineConvolutionVariables[0].x[i] *= 0;
	for (int i = 0; i < multigridSharpenningVariables[0].x.size(); i++)multigridSharpenningVariables[0].x[i] *= 0;

	stiffnessWeight = atof(prompt);

	if (UseDirectSolver.set) {
		sharpenningMatrix = mass + stiffness * stiffnessWeight;
	}

	if (!UpdateLinearSystem(stiffnessWeight, 1.0, hierarchy, multigridSharpenningCoefficients,
		deepMassCoefficients, deepStiffnessCoefficients,
		boundaryBoundaryMassMatrix, boundaryBoundaryStiffnessMatrix,
		boundaryDeepMassMatrix, boundaryDeepStiffnessMatrix,
		coarseSharpenningSolver, boundarySharpenningSolver, fineSharpenningSolver,
		sharpenningMatrix, Verbose.set, false, UseDirectSolver.set)) {
		printf("ERROR : Failed system update! \n");
	}

	if (UseDirectSolver.set) {
		ComputeExactSolution();
		UpdateOutputBuffer(exactSharpenningSolution);
	}
}

template<class Real>
void LineConvolution<Real>::ConvolutionWeightCallBack(Visualization* v, const char* prompt) {
	stopUpdate = false;
	for (int i = 0; i < multigridLineConvolutionVariables[0].x.size(); i++)multigridLineConvolutionVariables[0].x[i] *= 0;
	for (int i = 0; i < multigridSharpenningVariables[0].x.size(); i++)multigridSharpenningVariables[0].x[i] *= 0;

	convolutionWeight = atof(prompt);

	if (UseDirectSolver.set) {
		lineConvolutionMatrix = anisotropicMass + anisotropicStiffness * convolutionWeight;
	}

	if (!UpdateLinearSystem(convolutionWeight, 1.0,hierarchy, multigridLineConvolutionCoefficients,
		anisoDeepMassCoefficients, anisoDeepStiffnessCoefficients,
		anisoBoundaryBoundaryMassMatrix, anisoBoundaryBoundaryStiffnessMatrix,
		anisoBoundaryDeepMassMatrix, anisoBoundaryDeepStiffnessMatrix,
		coarseLineConvolutionSolver, boundaryLineConvolutionSolver, fineLineConvolutionSolver,
		lineConvolutionMatrix, Verbose.set, false, UseDirectSolver.set)) {
		printf("ERROR : Failed system update! \n");
	}

	if (UseDirectSolver.set) {
		ComputeExactSolution();
		UpdateOutputBuffer(exactSharpenningSolution);
	}


}

template<class Real>
void LineConvolution<Real>::SharpenningWeightCallBack(Visualization* v, const char* prompt) {
	stopUpdate = false;
	for (int i = 0; i < multigridLineConvolutionVariables[0].x.size(); i++)multigridLineConvolutionVariables[0].x[i] *= 0;
	for (int i = 0; i < multigridSharpenningVariables[0].x.size(); i++)multigridSharpenningVariables[0].x[i] *= 0;

	sharpenningWeight = atof(prompt);
	if (UseDirectSolver.set) {
		ComputeExactSolution();
		UpdateOutputBuffer(exactSharpenningSolution);
	}
}

template<class Real>
void LineConvolution<Real>::StopUpdateCallBack(Visualization* v, const char* prompt) {
	stopUpdate = !stopUpdate;
}

template<class Real>
void LineConvolution<Real>::ExportTextureCallBack(Visualization* v, const char* prompt) {

	Image<Point3D<float>> outputImage;
	outputImage.resize(textureWidth, textureHeight);
	for (int i = 0; i < outputImage.size(); i++) outputImage[i] = Point3D<float>(outputBuffer[3 * i], outputBuffer[3 * i + 1], outputBuffer[3 * i + 2]) / float(255.0);
	if(padding.nonTrivial) UnpaddImage(padding, outputImage);
	outputImage.write(prompt);
}

template<class Real>
int LineConvolution<Real>::UpdateSolution(bool verbose, bool detailVerbose) {
	clock_t begin;

	//(1)Update smoothed input solution
	if (verbose) begin = clock();
	VCycle(multigridLineConvolutionVariables, multigridLineConvolutionCoefficients, multigridIndices, boundaryLineConvolutionSolver, coarseLineConvolutionSolver, detailVerbose, detailVerbose);
	if (verbose) printf("Smoothing impulse %.4f \n", double(clock() - begin) / CLOCKS_PER_SEC);

	//(2) Compute sharpenning RHS
	mass_x0.resize(textureNodes.size());
	MultiplyBySystemMatrix_NoReciprocals(deepMassCoefficients, boundaryDeepMassMatrix, boundaryBoundaryMassMatrix, hierarchy.gridAtlases[0].boundaryGlobalIndex, hierarchy.gridAtlases[0].rasterLines, multigridLineConvolutionVariables[0].x, mass_x0);

	stiffness_x0.resize(textureNodes.size());
	MultiplyBySystemMatrix_NoReciprocals(deepStiffnessCoefficients, boundaryDeepStiffnessMatrix, boundaryBoundaryStiffnessMatrix, hierarchy.gridAtlases[0].boundaryGlobalIndex, hierarchy.gridAtlases[0].rasterLines, multigridLineConvolutionVariables[0].x, stiffness_x0);

	for (int i = 0; i < textureNodes.size(); i++) multigridSharpenningVariables[0].rhs[i] = mass_x0[i] + stiffnessWeight  * sharpenningWeight * stiffness_x0[i];

	//(3) Update geodesic distance solution	
	if (verbose) begin = clock();
	VCycle(multigridSharpenningVariables, multigridSharpenningCoefficients, multigridIndices, boundarySharpenningSolver, coarseSharpenningSolver, detailVerbose, detailVerbose);
	if (verbose) printf("Solving geodesic distance %.4f \n", double(clock() - begin) / CLOCKS_PER_SEC);

	return 1;
}


template<class Real>
int LineConvolution<Real>::InitializeSystem(const int width, const int height) {

	MultigridBlockInfo multigridBlockInfo(MultigridBlockWidth.value, MultigridBlockHeight.value, MultigridPaddedWidth.value, MultigridPaddedHeight.value, 0);
	if (!InitializeHierarchy(mesh, width, height, levels, textureNodes, cellIndices, hierarchy, atlasCharts,multigridBlockInfo, true, Verbose.set)) {
		printf("ERROR : Failed intialization! \n");
		return 0;
	}

	//Initialize node index
	nodeIndex.resize(width, height);
	for (int i = 0; i < nodeIndex.size(); i++)nodeIndex[i] = -1;
	for (int i = 0; i < textureNodes.size(); i++) {
		if (textureNodes[i].ci < 0 || textureNodes[i].ci > textureWidth - 1 || textureNodes[i].cj < 0 || textureNodes[i].cj > textureHeight - 1) {
			printf("Error: Invalid node! %d %d \n", textureNodes[i].ci, textureNodes[i].cj);
			return 0;
		}
		if (nodeIndex(textureNodes[i].ci, textureNodes[i].cj) != -1) {
			if (0)printf("WARNING: Multiple nodes mapped to pixel %d %d!\n", textureNodes[i].ci, textureNodes[i].cj);
		}
		nodeIndex(textureNodes[i].ci, textureNodes[i].cj) = i;
	}

	BoundaryProlongationData boundaryProlongation;
	if (!InitializeBoundaryProlongationData(hierarchy.gridAtlases[0], boundaryProlongation)) {
		printf("ERROR : Failed boundary prolongation! \n");
		return 0;
	}

	//////////////////////////////////// Initialize multigrid indices

	multigridIndices.resize(levels);
	for (int i = 0; i < levels; i++){
		const GridAtlas & gridAtlas = hierarchy.gridAtlases[i];
		multigridIndices[i].threadTasks = gridAtlas.threadTasks;
		multigridIndices[i].boundaryGlobalIndex = gridAtlas.boundaryGlobalIndex;
		multigridIndices[i].segmentedLines = gridAtlas.segmentedLines;
		multigridIndices[i].rasterLines = gridAtlas.rasterLines;
		multigridIndices[i].restrictionLines = gridAtlas.restrictionLines;
		multigridIndices[i].prolongationLines = gridAtlas.prolongationLines;
		if (i < levels - 1) {
			multigridIndices[i].boundaryRestriction = hierarchy.boundaryRestriction[i];
		}
	}

	//////////////////////////////////// Initialize multigrid coefficients

	//////////////////////////////////// 	Line Convolution coefficients
	{
		std::vector<Point2D<double>> vectorField;
		if (VectorField.set) {
			if (!ReadVector(vectorField, VectorField.value)) {
				printf("ERROR: Unable to read vector field! \n");
				return 0;
			}
		}
		if (!InitializeAnisotropicMetric(mesh, atlasCharts, parameterMetric, vectorField, MinorCurvature.set)) {
			printf("ERROR: Unable to initialize metric \n");
			return 0;
		} 

		std::vector<Point3D<double>> __inputSignal;
		std::vector<double> __texelToCellCoeffs;
		SparseMatrix<double, int> __boundaryCellBasedStiffnessRHSMatrix[3];

		if (!InitializeMassAndStiffness(anisoDeepMassCoefficients, anisoDeepStiffnessCoefficients,
			anisoBoundaryBoundaryMassMatrix, anisoBoundaryBoundaryStiffnessMatrix, anisoBoundaryDeepMassMatrix, anisoBoundaryDeepStiffnessMatrix,
			hierarchy, parameterMetric, atlasCharts, boundaryProlongation, false, __inputSignal, __texelToCellCoeffs, __boundaryCellBasedStiffnessRHSMatrix)) {
			printf("ERROR : Failed intialization! \n");
			return 0;
		}

		if (UseDirectSolver.set) {
			clock_t t_begin;
			t_begin = clock();
			FullMatrixConstruction(hierarchy.gridAtlases[0], anisoDeepMassCoefficients, anisoBoundaryBoundaryMassMatrix, anisoBoundaryDeepMassMatrix, anisotropicMass);
			FullMatrixConstruction(hierarchy.gridAtlases[0], anisoDeepStiffnessCoefficients, anisoBoundaryBoundaryStiffnessMatrix, anisoBoundaryDeepStiffnessMatrix, anisotropicStiffness);
			lineConvolutionMatrix = anisotropicMass + anisotropicStiffness * convolutionWeight;
			printf("Assembling matrices =  %.4f \n", double(clock() - t_begin) / CLOCKS_PER_SEC);
		}

		if(!UpdateLinearSystem(convolutionWeight, 1.0, hierarchy, multigridLineConvolutionCoefficients,
			anisoDeepMassCoefficients, anisoDeepStiffnessCoefficients,
			anisoBoundaryBoundaryMassMatrix, anisoBoundaryBoundaryStiffnessMatrix,
			anisoBoundaryDeepMassMatrix, anisoBoundaryDeepStiffnessMatrix,
			coarseLineConvolutionSolver, boundaryLineConvolutionSolver,fineLineConvolutionSolver,
			lineConvolutionMatrix, Verbose.set, true, UseDirectSolver.set)){
			printf("ERROR : Failed system update! \n");
			return 0;
		}
	}

	//////////////////////////////////// 	Sharpenning coefficients
	SparseMatrix< double, int > sharpenningMatrix;
	{
		if (!InitializeMetric(mesh, EMBEDDING_METRIC, atlasCharts, parameterMetric)) {
			printf("ERROR: Unable to initialize metric \n");
			return 0;
		}

		std::vector<Point3D<double>> __inputSignal;
		std::vector<double> __texelToCellCoeffs;
		SparseMatrix<double, int> __boundaryCellBasedStiffnessRHSMatrix[3];

		if (!InitializeMassAndStiffness(deepMassCoefficients, deepStiffnessCoefficients,
			boundaryBoundaryMassMatrix, boundaryBoundaryStiffnessMatrix, boundaryDeepMassMatrix, boundaryDeepStiffnessMatrix,
			hierarchy, parameterMetric, atlasCharts, boundaryProlongation, false, __inputSignal, __texelToCellCoeffs, __boundaryCellBasedStiffnessRHSMatrix)) {
			printf("ERROR : Failed intialization! \n");
			return 0;
		}

		if (UseDirectSolver.set) {
			clock_t t_begin;
			t_begin = clock();
			FullMatrixConstruction(hierarchy.gridAtlases[0], deepMassCoefficients, boundaryBoundaryMassMatrix, boundaryDeepMassMatrix, mass);
			FullMatrixConstruction(hierarchy.gridAtlases[0], deepStiffnessCoefficients, boundaryBoundaryStiffnessMatrix, boundaryDeepStiffnessMatrix, stiffness);
			sharpenningMatrix = mass + stiffness * stiffnessWeight;
			printf("Assembling matrices =  %.4f \n", double(clock() - t_begin) / CLOCKS_PER_SEC);
		}

		if (!UpdateLinearSystem(stiffnessWeight, 1.0, hierarchy, multigridSharpenningCoefficients,
			deepMassCoefficients, deepStiffnessCoefficients,
			boundaryBoundaryMassMatrix, boundaryBoundaryStiffnessMatrix,
			boundaryDeepMassMatrix, boundaryDeepStiffnessMatrix,
			coarseSharpenningSolver, boundarySharpenningSolver, fineSharpenningSolver,
			sharpenningMatrix, Verbose.set, true, UseDirectSolver.set)){
			printf("ERROR : Failed system update! \n");
			return 0;
		}
	}

	//////////////////////////////////// Initialize multigrid variables

	multigridLineConvolutionVariables.resize(levels);
	for (int i = 0; i < levels; i++) {
		MultigridLevelVariables<VECTOR_TYPE> & variables = multigridLineConvolutionVariables[i];
		variables.x.resize(hierarchy.gridAtlases[i].numTexels);
		variables.rhs.resize(hierarchy.gridAtlases[i].numTexels);
		variables.residual.resize(hierarchy.gridAtlases[i].numTexels);
		variables.boundary_rhs.resize(hierarchy.gridAtlases[i].boundaryGlobalIndex.size());
		variables.boundary_value.resize(hierarchy.gridAtlases[i].boundaryGlobalIndex.size());
		variables.variable_boundary_value.resize(hierarchy.gridAtlases[i].boundaryGlobalIndex.size());
	}

	multigridSharpenningVariables.resize(levels);
	for (int i = 0; i < levels; i++) {
		MultigridLevelVariables<VECTOR_TYPE> & variables = multigridSharpenningVariables[i];
		variables.x.resize(hierarchy.gridAtlases[i].numTexels);
		variables.rhs.resize(hierarchy.gridAtlases[i].numTexels);
		variables.residual.resize(hierarchy.gridAtlases[i].numTexels);
		variables.boundary_rhs.resize(hierarchy.gridAtlases[i].boundaryGlobalIndex.size());
		variables.boundary_value.resize(hierarchy.gridAtlases[i].boundaryGlobalIndex.size());
		variables.variable_boundary_value.resize(hierarchy.gridAtlases[i].boundaryGlobalIndex.size());
	}

	int numTexels = hierarchy.gridAtlases[0].numTexels;
	int numFineNodes = hierarchy.gridAtlases[0].numFineNodes;

	exactLineConvolutionSolution.resize(numTexels);
	exactSharpenningSolution.resize(numTexels);

	randSignal.resize(textureNodes.size());

	for (int i = 0; i < randSignal.size(); i++) {
		Point3D<float> randomColor = HSV2RGB(double(rand()) / double(RAND_MAX), 1, 1);
		randSignal[i] = SetData(randomColor[0], randomColor[1], randomColor[2]);
	}

	MultiplyBySystemMatrix_NoReciprocals(anisoDeepMassCoefficients, anisoBoundaryDeepMassMatrix, anisoBoundaryBoundaryMassMatrix, hierarchy.gridAtlases[0].boundaryGlobalIndex, hierarchy.gridAtlases[0].rasterLines, randSignal, multigridLineConvolutionVariables[0].rhs);
	
	for (int i = 0; i < multigridLineConvolutionVariables[0].x.size(); i++)multigridLineConvolutionVariables[0].x[i] *= 0;
	for (int i = 0; i < multigridSharpenningVariables[0].x.size(); i++)multigridSharpenningVariables[0].x[i] *= 0;

	if (UseDirectSolver.set) {
		ComputeExactSolution();
	}
	else {
		for (int i = 0; i < 4; i++) UpdateSolution();
	}
	return 1;
}

template<class Real>
void LineConvolution<Real>::InitializeVisualization(const int width, const int height) {

	outputBuffer = new unsigned char[height*width * 3];
	memset(outputBuffer, 204, height * width * 3 * sizeof(unsigned char));

	int tCount = mesh.triangles.size();

	visualization.triangles.resize(tCount);
	visualization.vertices.resize(3 * tCount);
	visualization.colors.resize(3 * tCount, Point3D<double>(0.75, 0.75, 0.75));
	visualization.textureCoordinates.resize(3 * tCount);
	visualization.normals.resize(3 * tCount);


	for (int i = 0; i < tCount; i++) for (int k = 0; k < 3; k++) visualization.triangles[i][k] = 3 * i + k;

	for (int i = 0; i<tCount; i++) {
		for (int j = 0; j < 3; j++) {
			visualization.vertices[3 * i + j] = mesh.vertices[mesh.triangles[i][j]];
			visualization.normals[3 * i + j] = mesh.normals[mesh.triangles[i][j]];
			visualization.textureCoordinates[3 * i + j] = mesh.textureCoordinates[3 * i + j];
		}
	}



	std::vector<int> boundaryEdges;
	if (!InitializeBoundaryEdges(mesh, boundaryEdges)) {
		printf("Unable to initialize boundary edges! \n");
	}

	for (int e = 0; e < boundaryEdges.size(); e++) {
		int tIndex = boundaryEdges[e] / 3;
		int kIndex = boundaryEdges[e] % 3;
		for (int c = 0; c < 2; c++) {
			Point3D<double> v = mesh.vertices[mesh.triangles[tIndex][(kIndex + c) % 3]];
			visualization.boundaryEdgeVertices.push_back(v);
		}
	}

	visualization.callBacks.push_back(Visualization::KeyboardCallBack(&visualization, 'f', "convolution weight", "Convolution Weight", ConvolutionWeightCallBack));
	visualization.callBacks.push_back(Visualization::KeyboardCallBack(&visualization, 'g', "sharpeness weight", "Sharpness Weight", SharpenningWeightCallBack));
	visualization.callBacks.push_back(Visualization::KeyboardCallBack(&visualization, 'y', "stiffness weight", "Stiffness Weight", StiffnessWeightCallBack));
	visualization.callBacks.push_back(Visualization::KeyboardCallBack(&visualization, 's', "export texture", "Output Texture", ExportTextureCallBack));
	visualization.callBacks.push_back(Visualization::KeyboardCallBack(&visualization, 'j', "toggle update", StopUpdateCallBack));


	visualization.UpdateVertexBuffer();
	visualization.UpdateFaceBuffer();
	visualization.UpdateTextureBuffer();

	if (UseDirectSolver.set) {
		UpdateOutputBuffer(exactSharpenningSolution);
	}
	else{
		UpdateOutputBuffer(multigridSharpenningVariables[0].x);
	}
}

template<class Real>
int LineConvolution<Real>::Init() {
	levels = Levels.value;
	textureWidth = Width.value;
	textureHeight = Height.value;
	sharpenningWeight = SharpenningWeight.value;
	stiffnessWeight = StiffnessWeight.value;
	convolutionWeight = ConvolutionWeight.value;

	if (!ReadTexturedMesh(mesh, Mesh.value, NULL)) {
		printf("Unable to read mesh data\n");
		return 0;
	}
	if (1) for (int i = 0; i < mesh.textureCoordinates.size(); i++)mesh.textureCoordinates[i][1] = 1.0 - mesh.textureCoordinates[i][1];

	if (RandomJitter.set) {
		printf("Jittering \n");
		double jitterScale = 1e-3 / double(std::max<int>(textureWidth, textureHeight));
		AddRandomJitter(mesh.textureCoordinates, jitterScale);
	}

	ComputePadding(padding, textureWidth, textureHeight, mesh.textureCoordinates);
	if (padding.nonTrivial) {
		PaddTextureCoordinates(padding, textureWidth, textureHeight, mesh.textureCoordinates);
		textureWidth += (padding.left + padding.right);
		textureHeight += (padding.bottom + padding.top);
	}

	//Define centroid and scale for visualization
	Point3D<double> centroid(0.f, 0.f, 0.f);
	for (int i = 0; i < mesh.vertices.size(); i++) centroid += mesh.vertices[i];
	centroid /= double(mesh.vertices.size());
	double radius = 0.f;
	for (int i = 0; i < mesh.vertices.size(); i++) radius = std::max<double>(radius, Point3D<double>::Length(mesh.vertices[i] - centroid));
	for (int i = 0; i < mesh.vertices.size(); i++) mesh.vertices[i] = (mesh.vertices[i] - centroid) / radius;


	printf("Resolution %d x %d \n", textureWidth, textureHeight);
	if (!InitializeSystem(textureWidth, textureHeight)) {
		printf("Unable to initialize system\n");
		return 0;
	}

	//Assign position to exterior nodes using baricentric-exponential map
	FEM::RiemannianMesh< double > * Rmesh = new FEM::RiemannianMesh< double >(GetPointer(mesh.triangles), mesh.triangles.size());
	Rmesh->setMetricFromEmbedding(GetPointer(mesh.vertices));
	Rmesh->makeUnitArea();
	Pointer(FEM::CoordinateXForm< double >) xForms = Rmesh->getCoordinateXForms();

	for (int i = 0; i<textureNodes.size(); i++) {
		if (textureNodes[i].tId != -1 && !textureNodes[i].isInterior) {
			FEM::HermiteSamplePoint< double > _p;
			_p.tIdx = textureNodes[i].tId;
			_p.p = Point2D< double >((double)1. / 3, (double)1. / 3);
			_p.v = textureNodes[i].baricentricCoords - _p.p;

			Rmesh->exp(xForms, _p);

			textureNodes[i].tId = _p.tIdx;
			textureNodes[i].baricentricCoords = _p.p;
		}
	}
	textureNodePositions.resize(textureNodes.size());
	return 1;
}

template<class Real>
int _main(int argc, char* argv[])
{
	if (!LineConvolution<Real>::Init()) return 0;
	if (!OutputTexture.set) {
		glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
		LineConvolution<Real>::visualization.displayMode = DisplayMode.value;
		if (DisplayMode.value == ONE_REGION_DISPLAY) {
			LineConvolution<Real>::visualization.screenWidth = 800;
			LineConvolution<Real>::visualization.screenHeight = 800;
		}
		else if (DisplayMode.value == TWO_REGION_DISPLAY) {
			LineConvolution<Real>::visualization.screenWidth = 1440;
			LineConvolution<Real>::visualization.screenHeight = 720;
		}
		LineConvolution<Real>::visualization.UpdateMainFrameSize();
		glutInitWindowSize(LineConvolution<Real>::visualization.screenWidth, LineConvolution<Real>::visualization.screenHeight);
		glutInit(&argc, argv);
		char windowName[1024];
		sprintf(windowName, "Line Integral Convolution");
		glutCreateWindow(windowName);
		if (glewInit() != GLEW_OK) fprintf(stderr, "[ERROR] glewInit failed\n"), exit(0);
		glutDisplayFunc(LineConvolution<Real>::Display);
		glutReshapeFunc(LineConvolution<Real>::Reshape);
		glutMouseFunc(LineConvolution<Real>::MouseFunc);
		glutMotionFunc(LineConvolution<Real>::MotionFunc);
		glutKeyboardFunc(LineConvolution<Real>::KeyboardFunc);
		if (!UseDirectSolver.set)glutIdleFunc(LineConvolution<Real>::Idle);
		if (CameraConfig.set) LineConvolution<Real>::visualization.ReadSceneConfigurationCallBack(&LineConvolution<Real>::visualization, CameraConfig.value);
		LineConvolution<Real>::InitializeVisualization(LineConvolution<Real>::textureWidth, LineConvolution<Real>::textureHeight);
		glutMainLoop();
	}
	else {
		LineConvolution<Real>::ExportTextureCallBack(&LineConvolution<Real>::visualization, OutputTexture.value);
	}
	return 1;
}

int main(int argc, char* argv[])
{
	cmdLineParse(argc - 1, argv + 1, params);
	if (!Mesh.set) { ShowUsage(argv[0]); return EXIT_FAILURE; }
	omp_set_num_threads(Threads.value);
	printf("Number of threads %d. Number of processors %d. \n", omp_get_max_threads(), omp_get_num_procs());
	_main<REAL_TYPE>(argc, argv);
	return 0;
}
