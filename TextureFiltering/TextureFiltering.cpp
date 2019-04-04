#define USE_FLOAT 1
#include <Src/VectorPrecisionType.inl>

#include <Misha/CmdLineParser.h> 
#include <Src/SimpleMesh.h>
#include <Src/Basis.h>
#include <Misha/FEM.h>
#include <Src/Solver.h>
#include <Src/Hierarchy.h>
#include <Src/MassAndStiffness.h>
#include <Src/InteriorTexelToCellLines.inl>
#include <Src/Padding.h>
#include <Src/TextureFilteringVisualization.h>
#include <Src/RandomJitter.h>

cmdLineParameter<  char* > Mesh("mesh");
cmdLineParameter<  char* > InputTexture("inTexture");
cmdLineParameter<  char* > OutputTexture("outTexture");
cmdLineParameter<  float > StiffnessWeight("stiffness", 1e-3);
cmdLineParameter<  float > SharpnessWeight("sharpness", 1.0);
cmdLineParameter<  int   > DisplayMode("display",FOUR_REGION_DISPLAY);
cmdLineParameter<  int   > Threads("threads", omp_get_num_procs());
cmdLineParameter<  int   > Levels("levels", 4);

cmdLineParameter<  int   > MultigridBlockHeight("mBlockH", 16);
cmdLineParameter<  int   > MultigridBlockWidth ("mBlockW", 128);
cmdLineParameter<  int   > MultigridPaddedHeight("mPadH", 0);
cmdLineParameter<  int   > MultigridPaddedWidth("mPadW", 2);
cmdLineParameter<  int   > MultigridUpdateVcycles("mVCycles", 6);

cmdLineReadable RandomJitter("jitter");
cmdLineParameter<  char* > CameraConfig("camera");
cmdLineReadable UseDirectSolver("useDirectSolver");
cmdLineReadable Verbose("verbose");
cmdLineReadable DetailVerbose("detail");
cmdLineReadable* params[] =
{
	&Mesh, &InputTexture, &OutputTexture, &StiffnessWeight, &SharpnessWeight, &CameraConfig, &Levels,&UseDirectSolver,&Threads,&DisplayMode,&Verbose,
	&DetailVerbose,&MultigridBlockHeight,&MultigridBlockWidth,&MultigridPaddedHeight,&MultigridPaddedWidth,&MultigridUpdateVcycles,&RandomJitter,
	NULL
};

void ShowUsage(const char* ex)
{
	printf("Usage %s:\n", ex);


	printf("\t --%s <input mesh>\n", Mesh.name);
	printf("\t --%s <input texture>\n", InputTexture.name);
	printf("\t --%s <output texture>\n", OutputTexture.name);

	printf("\t --%s <stiffness weight> [%f]\n", StiffnessWeight.name, StiffnessWeight.value);
	printf("\t --%s <sharpness weight> [%f]\n", SharpnessWeight.name, SharpnessWeight.value);
	printf("\t --%s <camera configuration >\n", CameraConfig.name);
	printf("\t --%s <hierarchy levels> [%d]\n", Levels.name,Levels.value);
	printf("\t --%s <threads> [%d]\n", Threads.name, Threads.value);
	printf("\t --%s <disable update>\n", UseDirectSolver.name);
	printf("\t --%s <verbose>\n", Verbose.name);
	printf("\t --%s <detail verbose>\n", DetailVerbose.name);
	printf("\t --%s <display mode> [%d]\n", DisplayMode.name, DisplayMode.value);
	printf("\t \t [%d] One Region \n", ONE_REGION_DISPLAY);
	printf("\t \t [%d] Two Region \n", TWO_REGION_DISPLAY);
	printf("\t \t [%d] Three Region \n", THREE_REGION_DISPLAY);
	printf("\t \t [%d] Four Region \n", FOUR_REGION_DISPLAY);

	printf("\t --%s <multigrid block width> [%d]\n", MultigridBlockWidth.name, MultigridBlockWidth.value);
	printf("\t --%s <multigrid block height> [%d]\n", MultigridBlockHeight.name, MultigridBlockHeight.value);
	printf("\t --%s <multigrid padded width> [%d]\n", MultigridPaddedWidth.name, MultigridPaddedWidth.value);
	printf("\t --%s <multigrid padded height> [%d]\n", MultigridPaddedHeight.name, MultigridPaddedHeight.value);
	printf("\t --%s <multigrid update VCycles> [%d]\n", MultigridUpdateVcycles.name, MultigridUpdateVcycles.value);

	printf("\t --%s <random jitter>\n", RandomJitter.name);
}

enum {
	INPUT_TEXTURE,
	OUTPUT_TEXTURE,
	TEXTURE_COUNT
};

template<class Real>
class TextureFilter{
public:
	static TexturedMesh mesh;
	static int textureWidth;
	static int textureHeight;
	static double stiffnessWeight;
	static double sharpnessWeight;
	static int levels;

	static HierarchicalSystem hierarchy;
	static bool sharpeningWeightUpdated;
	static bool positiveSharpenning;

	static Image<Point3D<float>> filteredTexture;
	//UI
	static std::vector<Real> unifromTexelSharpenningMask;
	static std::vector<Real> cellSharpenningMask;
	static std::vector<Real> uniformCellSharpenningMask;
	static std::vector<Real> texelStiffness[3];
	static std::vector<Point3D<float>> cellCenterPositions;

	static std::vector<Point3D<float>>textureNodePositions;

	static std::vector<AtlasChart> atlasCharts;

	static std::vector<CellIndex> cellIndices;
	
	static std::vector<TextureNodeInfo> textureNodes;

	static SparseMatrix<double,int> mass;
	static SparseMatrix<double,int> stiffness;
	static SparseMatrix< double, int > filteringMatrix;
	//RHS computation
	static std::vector<InteriorTexelToCellLine> interiorTexelToCellLines;
	static std::vector<VECTOR_TYPE> interiorTexelToCellCoeffs;
	static SparseMatrix<Real, int> boundaryCellBasedStiffnessRHSMatrix[3];
	static std::vector<Real> boundaryTexelStiffness[3];
	static std::vector<VECTOR_TYPE> texelModulatedStiffness;

	static std::vector<VECTOR_TYPE> mass_x0;
	static std::vector<VECTOR_TYPE> stiffness_x0;
	static std::vector<VECTOR_TYPE> exactSolution;


	static std::vector<MultigridLevelCoefficients<Real>> multigridFilteringCoefficients;
	static std::vector<MultigridLevelVariables<VECTOR_TYPE>> multigridFilteringVariables;
	static std::vector<MultigridLevelIndices<Real>> multigridIndices;

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
	static BoundarySolverType boundarySolver;
	static CoarseSolverType coarseSolver;
	static DirectSolverType directSolver;


	//Linear Operators
	static std::vector<double> deepMassCoefficients;
	static std::vector<double> deepStiffnessCoefficients;
	static SparseMatrix<double, int> boundaryBoundaryMassMatrix;
	static SparseMatrix<double, int> boundaryBoundaryStiffnessMatrix;
	static SparseMatrix<double, int> boundaryDeepMassMatrix;
	static SparseMatrix<double, int> boundaryDeepStiffnessMatrix;

	static int updateVCycles;
	static int cycleCount;

	static Padding padding;

	//Visulization
	static TextureFilteringVisualization visualization;
	static void ExportTextureCallBack(Visualization* v, const char* prompt);

	static void SharpnessWeightCallBack(Visualization* v, const char* prompt);
	static void StiffnessWeightCallBack(Visualization* v, const char* prompt);

	static int Init();
	static void InitializeVisualization();
	static int UpdateSolution(bool verbose = false, bool detailVerbose = false);
	static void ComputeExactSolution(bool verbose = false);
	static int InitializeSystem(const int width, const int height);

	static void UpdateFilteredColorTexture(const std::vector<VECTOR_TYPE> & solution);

	static void UpdateFilteredTexture(const std::vector<VECTOR_TYPE> & solution);
	static void UpdateMaskTexture();


	static void Display(void){ visualization.Display(); }
	static void MouseFunc(int button, int state, int x, int y);
	static void MotionFunc(int x, int y);
	static void Reshape(int w, int h){ visualization.Reshape(w, h); }
	static void KeyboardFunc(unsigned char key, int x, int y){ visualization.KeyboardFunc(key, x, y); }
	static void Idle();
};

template<class Real> TexturedMesh												TextureFilter<Real>::mesh;
template<class Real> int														TextureFilter<Real>::textureWidth;
template<class Real> int														TextureFilter<Real>::textureHeight;
template<class Real> TextureFilteringVisualization								TextureFilter<Real>::visualization;
template<class Real> SparseMatrix<double,int>									TextureFilter<Real>::mass;
template<class Real> SparseMatrix<double,int>									TextureFilter<Real>::stiffness;
template<class Real> SparseMatrix<double, int>									TextureFilter<Real>::filteringMatrix;

template<class Real> double														TextureFilter<Real>::stiffnessWeight;
template<class Real> double														TextureFilter<Real>::sharpnessWeight;

template<class Real> std::vector<TextureNodeInfo>								TextureFilter<Real>::textureNodes;
template<class Real> std::vector<CellIndex>										TextureFilter<Real>::cellIndices;

template<class Real> int														TextureFilter<Real>::levels;
template<class Real> HierarchicalSystem											TextureFilter<Real>::hierarchy;


template<class Real> bool TextureFilter<Real>::sharpeningWeightUpdated = true;
template<class Real> bool TextureFilter<Real>::positiveSharpenning = true;

//UI
template<class Real> std::vector<Point3D<float>>										TextureFilter<Real>::cellCenterPositions;
template<class Real> std::vector<Real>													TextureFilter<Real>::texelStiffness[3];
template<class Real> std::vector<Real>													TextureFilter<Real>::cellSharpenningMask;
template<class Real> std::vector<Real>													TextureFilter<Real>::uniformCellSharpenningMask;

template<class Real> Image<Point3D<float>>												TextureFilter<Real>::filteredTexture;
template<class Real> std::vector<VECTOR_TYPE>											TextureFilter<Real>::stiffness_x0;
template<class Real> std::vector<VECTOR_TYPE>											TextureFilter<Real>::mass_x0;
template<class Real> std::vector<VECTOR_TYPE>											TextureFilter<Real>::exactSolution;

template<class Real> std::vector<MultigridLevelCoefficients<Real>>						TextureFilter<Real>::multigridFilteringCoefficients;
template<class Real> std::vector<MultigridLevelVariables<VECTOR_TYPE>>					TextureFilter<Real>::multigridFilteringVariables;
template<class Real> std::vector<MultigridLevelIndices<Real>>							TextureFilter<Real>::multigridIndices;


template<class Real> typename TextureFilter<Real>::BoundarySolverType					TextureFilter<Real>::boundarySolver;
template<class Real> typename TextureFilter<Real>::CoarseSolverType						TextureFilter<Real>::coarseSolver;
template<class Real> typename TextureFilter<Real>::DirectSolverType						TextureFilter<Real>::directSolver;

template<class Real> std::vector<AtlasChart>											TextureFilter<Real>::atlasCharts;

template<class Real> std::vector<InteriorTexelToCellLine>								TextureFilter<Real>::interiorTexelToCellLines;
template<class Real> std::vector<VECTOR_TYPE>											TextureFilter<Real>::interiorTexelToCellCoeffs;
template<class Real> SparseMatrix<Real, int>											TextureFilter<Real>::boundaryCellBasedStiffnessRHSMatrix[3];
template<class Real> std::vector<Real>													TextureFilter<Real>::boundaryTexelStiffness[3];
template<class Real> std::vector<VECTOR_TYPE>											TextureFilter<Real>::texelModulatedStiffness;

template<class Real> std::vector<Point3D<float>>										TextureFilter<Real>::textureNodePositions;
template<class Real> std::vector<Real>													TextureFilter<Real>::unifromTexelSharpenningMask;

template<class Real> Padding															TextureFilter<Real>::padding;


template<class Real>  std::vector<double>												TextureFilter<Real>::deepMassCoefficients;
template<class Real>  std::vector<double>												TextureFilter<Real>::deepStiffnessCoefficients;
template<class Real>  SparseMatrix<double, int>											TextureFilter<Real>::boundaryBoundaryMassMatrix;
template<class Real>  SparseMatrix<double, int>											TextureFilter<Real>::boundaryBoundaryStiffnessMatrix;
template<class Real>  SparseMatrix<double, int>											TextureFilter<Real>::boundaryDeepMassMatrix;
template<class Real>  SparseMatrix<double, int>											TextureFilter<Real>::boundaryDeepStiffnessMatrix;

template<class Real> int																TextureFilter<Real>::updateVCycles;
template<class Real> int																TextureFilter<Real>::cycleCount;

template<class Real>
void TextureFilter<Real>::UpdateMaskTexture(){

#pragma omp parallel for
	for (int i = 0; i < textureNodes.size(); i++){
		float texelSharpenningValue = unifromTexelSharpenningMask[i];
		if (texelSharpenningValue != 0.5){
			Point3D<float> newColor;
			if (texelSharpenningValue > 0.5) {
				texelSharpenningValue = 2.0 * texelSharpenningValue - 1.0;
				newColor = Point3D<float>(1.f, 0.f, 0.f) * texelSharpenningValue + Point3D<float>(0.8f, 0.8f, 0.8f) * (1.0 - texelSharpenningValue);
			}
			else {
				texelSharpenningValue = 2.0 *texelSharpenningValue;
				newColor = Point3D<float>(0.f, 0.f, 1.f) *(1.0 - texelSharpenningValue) + Point3D<float>(0.8f, 0.8f, 0.8f) * texelSharpenningValue;
			}
			int ci = textureNodes[i].ci;
			int cj = textureNodes[i].cj;
			int offset = 3 * (textureWidth*cj + ci);
			visualization.maskBufferValues[offset + 0] = (unsigned char)(newColor[0] * 255.0);
			visualization.maskBufferValues[offset + 1] = (unsigned char)(newColor[1] * 255.0);
			visualization.maskBufferValues[offset + 2] = (unsigned char)(newColor[2] * 255.0);
		}
	}

	visualization.UpdateMaskTextureBuffer();
}

template<class Real>
void TextureFilter<Real>::UpdateFilteredColorTexture(const std::vector<VECTOR_TYPE> & solution) {
#pragma omp parallel for
	for (int i = 0; i < textureNodes.size(); i++) {
		int ci = textureNodes[i].ci;
		int cj = textureNodes[i].cj;
		int offset = 3 * (textureWidth*cj + ci);
		for (int c = 0; c < 3; c++) {
			double value = std::min<double>(1.0, std::max<double>(0, solution[i][c]));
			visualization.colorTextureBuffer[offset + c] = (unsigned char)(value*255.0);
		}
	}
}

template<class Real>
void TextureFilter<Real>::UpdateFilteredTexture(const std::vector<VECTOR_TYPE> & solution){
#pragma omp parallel for
	for (int i = 0; i < textureNodes.size(); i++) {
		int ci = textureNodes[i].ci;
		int cj = textureNodes[i].cj;
		filteredTexture(ci,cj) = Point3D<float>(solution[i][0], solution[i][1], solution[i][2]);
	}
}

template<class Real>
void TextureFilter<Real>::Idle() {

	auto RescaleFunction = [](double x) {
		//return x*(2.0 *x + 1);
		return 2.0 *x;
	};

	float radius = 0.1;
	float sharpenningVartiation = 0.2;
	if (visualization.isBrushActive) {
		Point3D<float> selectedPoint;
		bool validSelection = false;
		if (visualization.showMesh) {
			validSelection = visualization.select(visualization.diskX, visualization.diskY, selectedPoint);
		}
		if (validSelection) {
			
			float sharpenningSign = positiveSharpenning ? 1.f : -1.f;
			sharpenningVartiation *= sharpenningSign;
			for (int i = 0; i < cellIndices.size(); i++) {
				float distanceRatio = Point3D<float>::Length(cellCenterPositions[i] - selectedPoint) / radius;
				float factor = 1.0 - distanceRatio;
				factor = factor < 0 ? 0 : factor*factor*(-2.0*factor + 3.0);
				Real uniformSharpenningMaskValue = std::max<Real>(0, std::min<Real>(1.0, uniformCellSharpenningMask[i] + sharpenningVartiation * factor));
				uniformCellSharpenningMask[i] = uniformSharpenningMaskValue;
				cellSharpenningMask[i] = RescaleFunction(uniformSharpenningMaskValue);
			}
			if (1){
				for (int i = 0; i < textureNodePositions.size(); i++) {
					float distanceRatio = Point3D<float>::Length(textureNodePositions[i] - selectedPoint) / radius;
					float factor = 1.0 - distanceRatio;
					factor = factor < 0 ? 0 : factor*factor*(-2.0*factor + 3.0);
					Real sharpenningMaskValue = std::max<Real>(0, std::min<Real>(1.0, unifromTexelSharpenningMask[i] + sharpenningVartiation * factor));
					unifromTexelSharpenningMask[i] = sharpenningMaskValue;
				}
				UpdateMaskTexture();
			}

			cycleCount = 0;
		}
	}
	else if (visualization.isSlideBarActive) {
		if (visualization.slideBarCursorOldPosition != visualization.slideBarCursorPosition) {

			Real diff = (Real)(visualization.slideBarCursorPosition - visualization.slideBarCursorOldPosition);
			visualization.slideBarCursorOldPosition = visualization.slideBarCursorPosition;

			for (int i = 0; i < cellIndices.size(); i++) {
				Real uniformSharpenningMaskValue = std::max<Real>(0, std::min<Real>(1.0, uniformCellSharpenningMask[i] + diff));
				uniformCellSharpenningMask[i] = uniformSharpenningMaskValue;
				cellSharpenningMask[i] = RescaleFunction(uniformSharpenningMaskValue);
			}

			if (1) {
				for (int i = 0; i < textureNodePositions.size(); i++) {
					Real sharpenningMaskValue = std::max<Real>(0, std::min<Real>(1.0, unifromTexelSharpenningMask[i] + diff));
					unifromTexelSharpenningMask[i] = sharpenningMaskValue;
				}
				UpdateMaskTexture();
			}
			cycleCount = 0;
		}
	}

	if (cycleCount < updateVCycles) {
		if (!UpdateSolution()) {
			printf("Updated solution failed! \n");
		}

		if (visualization.textureType == COLOR_TEXTURE) {
			UpdateFilteredColorTexture(multigridFilteringVariables[0].x);
			visualization.UpdateColorTextureBuffer();
		}
		else {
			UpdateFilteredTexture(multigridFilteringVariables[0].x);
			visualization.UpdateTextureBuffer(filteredTexture);
		}
		cycleCount++;
	}
	
}

template<class Real>
void TextureFilter<Real>::MouseFunc(int button, int state, int x, int y) {

	visualization.newX = x; visualization.newY = y;
	visualization.rotating = visualization.scaling = visualization.panning = false;
	visualization.isSlideBarActive = false;
	visualization.isBrushActive = false;

	if (state == GLUT_DOWN && glutGetModifiers() & GLUT_ACTIVE_SHIFT) {
		visualization.isBrushActive = true;
		sharpeningWeightUpdated = false;
		visualization.diskX = x;
		visualization.diskY = y;

		if (button == GLUT_LEFT_BUTTON) {
			positiveSharpenning = true;
		}
		else if (button == GLUT_RIGHT_BUTTON) positiveSharpenning = false;
	}
	else if (visualization.showSlideBar && x > 10 && x < visualization._screenWidth - 10 && (visualization._screenHeight - y) > 18 && (visualization._screenHeight - y) < 32) {//Slide bar update
		visualization.isSlideBarActive = true;
		sharpeningWeightUpdated = false;
		double slideBarCursorPosition = double(x - 20) / double(visualization._screenWidth - 40);
		slideBarCursorPosition = std::min<double>(std::max<double>(slideBarCursorPosition, 0), 1.0);
		visualization.slideBarCursorPosition = slideBarCursorPosition;
		
	}
	else {
		if (visualization.showMesh) {
			visualization.newX = x; visualization.newY = y;

			visualization.rotating = visualization.scaling = visualization.panning = false;
			if (button == GLUT_LEFT_BUTTON) {
				if (glutGetModifiers() & GLUT_ACTIVE_CTRL) visualization.panning = true;
				else                                        visualization.rotating = true;
			}
			else if (button == GLUT_RIGHT_BUTTON) visualization.scaling = true;
		}
		else {

		}
	}

	glutPostRedisplay();
}

template<class Real>
void TextureFilter<Real>::MotionFunc(int x, int y) {

	if (visualization.isBrushActive) {
		sharpeningWeightUpdated = false;
		visualization.diskX = x;
		visualization.diskY = y;
	}
	else if (visualization.showSlideBar && visualization.isSlideBarActive) {
		sharpeningWeightUpdated = false;
		double slideBarCursorPosition = double(x - 20) / double(visualization._screenWidth - 40);
		slideBarCursorPosition = std::min<double>(std::max<double>(slideBarCursorPosition, 0), 1.0);
		visualization.slideBarCursorPosition = slideBarCursorPosition;
	}
	else {
		if (visualization.showMesh) {
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
		else {
			visualization.oldX = visualization.newX, visualization.oldY = visualization.newY, visualization.newX = x, visualization.newY = y;

			int imageSize = std::min< int >(visualization.screenWidth, visualization.screenHeight);
			if (visualization.panning) visualization.xForm.offset[0] -= (visualization.newX - visualization.oldX) / visualization.imageToScreenScale(), visualization.xForm.offset[1] += (visualization.newY - visualization.oldY) / visualization.imageToScreenScale();
			else
			{
				float dz = float(pow(1.1, double(visualization.newY - visualization.oldY) / 8));
				visualization.xForm.zoom *= dz;
			}
		}
	}
	glutPostRedisplay();
}

template<class Real>
void TextureFilter<Real>::ExportTextureCallBack(Visualization* v, const char* prompt){
	UpdateFilteredTexture(multigridFilteringVariables[0].x);
	Image<Point3D<float>> outputTexture = filteredTexture;
	if (padding.nonTrivial) UnpaddImage(padding, outputTexture);

	char* ext = GetFileExtension(prompt);
	if (!strcasecmp(ext, "png")) {
		if (visualization.textureType == COLOR_TEXTURE) {
			outputTexture.write(prompt);
		}
		else if (visualization.textureType == NORMAL_TEXTURE) {
			for (int i = 0; i < outputTexture.size(); i++) {
				outputTexture[i] /= Point3D<float>::Length(outputTexture[i]);
				outputTexture[i] = outputTexture[i] * 0.5f + Point3D<float>(0.5, 0.5, 0.5);
			}
			outputTexture.write(prompt);
		}
	}
	else {
		WriteBinaryImage(outputTexture, prompt);
	}
}

template<class Real>
void  TextureFilter<Real>::SharpnessWeightCallBack(Visualization* v, const char* prompt){
	cycleCount = 0;

	sharpnessWeight = atof(prompt);
	for (int i = 0; i <multigridFilteringVariables[0].rhs.size(); i++) multigridFilteringVariables[0].rhs[i] = mass_x0[i] + stiffness_x0[i] * sharpnessWeight*stiffnessWeight;

	if (UseDirectSolver.set) {
		ComputeExactSolution();
		if (visualization.textureType == COLOR_TEXTURE) {
			UpdateFilteredColorTexture(exactSolution);
			visualization.UpdateColorTextureBuffer();
		}
		else {
			UpdateFilteredTexture(exactSolution);
			visualization.UpdateTextureBuffer(filteredTexture);
		}
	}
}

template<class Real>
void  TextureFilter<Real>::StiffnessWeightCallBack(Visualization* v, const char* prompt) {
	cycleCount = 0;
	
	stiffnessWeight = atof(prompt);
	if (UseDirectSolver.set) {
		filteringMatrix = mass + stiffness * stiffnessWeight;
	}
	if (!UpdateLinearSystem(stiffnessWeight, 1.0, hierarchy, multigridFilteringCoefficients,
		deepMassCoefficients, deepStiffnessCoefficients,
		boundaryBoundaryMassMatrix, boundaryBoundaryStiffnessMatrix,
		boundaryDeepMassMatrix, boundaryDeepStiffnessMatrix,
		coarseSolver, boundarySolver, directSolver,
		filteringMatrix, Verbose.set, false, UseDirectSolver.set)) {
		printf("ERROR : Failed system update! \n");
	}
	for (int i = 0; i <multigridFilteringVariables[0].rhs.size(); i++) multigridFilteringVariables[0].rhs[i] = mass_x0[i] + stiffness_x0[i] * sharpnessWeight*stiffnessWeight;
	
	if (UseDirectSolver.set){
		ComputeExactSolution();
		if (visualization.textureType == COLOR_TEXTURE) {
			UpdateFilteredColorTexture(exactSolution);
			visualization.UpdateColorTextureBuffer();
		}
		else {
			UpdateFilteredTexture(exactSolution);
			visualization.UpdateTextureBuffer(filteredTexture);
		}
	}
}


template<class Real>
void TextureFilter<Real>::ComputeExactSolution(bool verbose) {
	clock_t t_begin;
	if(verbose)t_begin = clock();
	solve(directSolver, exactSolution, multigridFilteringVariables[0].rhs);
	if (verbose) printf("Solving time =  %.4f \n", double(clock() - t_begin) / CLOCKS_PER_SEC);
}

template<class Real>
int TextureFilter<Real>::UpdateSolution(bool verbose, bool detailVerbose){
	if (!sharpeningWeightUpdated) {
		int numTexels = multigridFilteringVariables[0].rhs.size();
		clock_t p_begin;
		if (verbose) p_begin = clock();

		CellStiffnessToTexelStiffness(cellSharpenningMask, interiorTexelToCellLines, interiorTexelToCellCoeffs, boundaryCellBasedStiffnessRHSMatrix, boundaryTexelStiffness, hierarchy.gridAtlases[0].boundaryGlobalIndex, texelModulatedStiffness);
		for (int i = 0; i < numTexels; i++)  multigridFilteringVariables[0].rhs[i] = mass_x0[i] + texelModulatedStiffness[i] * stiffnessWeight;

		if (verbose) printf("RHS update time %.4f  \n", double(clock() - p_begin) / CLOCKS_PER_SEC);	
		sharpeningWeightUpdated = true;
	}

	VCycle(multigridFilteringVariables, multigridFilteringCoefficients, multigridIndices, boundarySolver, coarseSolver, verbose, detailVerbose);

	return 1;
}

template<class Real>
int TextureFilter<Real>::InitializeSystem(const int width, const int height) {

	clock_t t_begin;

	MultigridBlockInfo multigridBlockInfo(MultigridBlockWidth.value, MultigridBlockHeight.value,MultigridPaddedWidth.value,MultigridPaddedHeight.value, 0);
	if (!InitializeHierarchy(mesh, width, height, levels, textureNodes,cellIndices, hierarchy, atlasCharts, multigridBlockInfo,true,Verbose.set)){
		printf("ERROR : Failed intialization! \n");
		return 0;
	}

	BoundaryProlongationData boundaryProlongation;
	if (!InitializeBoundaryProlongationData(hierarchy.gridAtlases[0], boundaryProlongation)){
		printf("ERROR : Failed boundary prolongation! \n");
		return 0;
	}

	std::vector<VECTOR_TYPE> _x0;
	_x0.resize(textureNodes.size());
	
	for (int i = 0; i < textureNodes.size(); i++) {
		Point3D<float> texelValue = mesh.texture(textureNodes[i].ci, textureNodes[i].cj);
		_x0[i] = SetData(texelValue[0], texelValue[1], texelValue[2]);
	}

	std::vector<Point3D<double>> inputSignal(textureNodes.size());
	for (int i = 0; i < textureNodes.size(); i++) inputSignal[i] = mesh.texture(textureNodes[i].ci, textureNodes[i].cj);
	std::vector<double> texelToCellCoeffs;


	printf("Fine system construction ... \n");
	t_begin = clock();

	std::vector<std::vector<SquareMatrix<double, 2>>> parameterMetric;
	if (!InitializeMetric(mesh, EMBEDDING_METRIC, atlasCharts, parameterMetric)) {
		printf("ERROR: Unable to initialize metric \n");
		return 0;
	}

#if USE_FLOAT
	SparseMatrix<double, int> _boundaryCellBasedStiffnessRHSMatrix[3];
	if (!InitializeMassAndStiffness(deepMassCoefficients, deepStiffnessCoefficients,
		boundaryBoundaryMassMatrix, boundaryBoundaryStiffnessMatrix, boundaryDeepMassMatrix, boundaryDeepStiffnessMatrix,
		hierarchy, parameterMetric, atlasCharts,boundaryProlongation, true, inputSignal, texelToCellCoeffs, _boundaryCellBasedStiffnessRHSMatrix)){
		printf("ERROR : Failed intialization! \n");
		return 0;
	}
	for (int c = 0; c < 3; c++) boundaryCellBasedStiffnessRHSMatrix[c] = _boundaryCellBasedStiffnessRHSMatrix[c];
#else
	if (!InitializeMassAndStiffness(deepMassCoefficients, deepStiffnessCoefficients,
		boundaryBoundaryMassMatrix, boundaryBoundaryStiffnessMatrix, boundaryDeepMassMatrix, boundaryDeepStiffnessMatrix,
		hierarchy, parameterMetric, atlasCharts, boundaryProlongation, true, inputSignal, texelToCellCoeffs, boundaryCellBasedStiffnessRHSMatrix)) {
		printf("ERROR : Failed intialization! \n");
		return 0;
	}
#endif

	interiorTexelToCellCoeffs.resize(4 * hierarchy.gridAtlases[0].numDeepTexels);
	for (int i = 0; i < 4 * hierarchy.gridAtlases[0].numDeepTexels; i++)interiorTexelToCellCoeffs[i] = SetData(Real(texelToCellCoeffs[3*i + 0]), Real(texelToCellCoeffs[3*i + 1]), Real(texelToCellCoeffs[3*i + 2]));
	
	printf("Fine system construction =  %.4f \n", double(clock() - t_begin) / CLOCKS_PER_SEC);

	if (!InitializeInteriorTexelToCellLines(interiorTexelToCellLines, hierarchy.gridAtlases[0])) {
		printf("ERROR: Interior texel to cell not initialized! \n");
		return 0;
	}

	for (int c = 0; c < 3; c++) boundaryTexelStiffness[c].resize(hierarchy.gridAtlases[0].boundaryGlobalIndex.size());
	texelModulatedStiffness.resize(hierarchy.gridAtlases[0].numTexels);

	if (UseDirectSolver.set){
		FullMatrixConstruction(hierarchy.gridAtlases[0], deepMassCoefficients, boundaryBoundaryMassMatrix, boundaryDeepMassMatrix, mass);
		FullMatrixConstruction(hierarchy.gridAtlases[0], deepStiffnessCoefficients, boundaryBoundaryStiffnessMatrix, boundaryDeepStiffnessMatrix, stiffness);
		filteringMatrix  = mass + stiffness * stiffnessWeight;
	}

	multigridIndices.resize(levels);
	for (int i = 0; i < levels; i++) {
		const GridAtlas & gridAtlas = hierarchy.gridAtlases[i];
		multigridIndices[i].threadTasks = gridAtlas.threadTasks;
		multigridIndices[i].boundaryGlobalIndex = gridAtlas.boundaryGlobalIndex;
		multigridIndices[i].segmentedLines = gridAtlas.segmentedLines;
		multigridIndices[i].rasterLines = gridAtlas.rasterLines;
		multigridIndices[i].restrictionLines = gridAtlas.restrictionLines;
		multigridIndices[i].prolongationLines = gridAtlas.prolongationLines;
		if (i < levels - 1){
			multigridIndices[i].boundaryRestriction = hierarchy.boundaryRestriction[i];
		}
	}

	if (!UpdateLinearSystem(stiffnessWeight, 1.0, hierarchy, multigridFilteringCoefficients,
		deepMassCoefficients, deepStiffnessCoefficients,
		boundaryBoundaryMassMatrix, boundaryBoundaryStiffnessMatrix,
		boundaryDeepMassMatrix, boundaryDeepStiffnessMatrix,
		coarseSolver, boundarySolver, directSolver,
		filteringMatrix, Verbose.set, true, UseDirectSolver.set)) {
		printf("ERROR : Failed system update! \n");
		return 0;
	}

	multigridFilteringVariables.resize(levels);
	for (int i = 0; i < levels; i++){
		MultigridLevelVariables<VECTOR_TYPE> & variables = multigridFilteringVariables[i];
		variables.x.resize(hierarchy.gridAtlases[i].numTexels);
		variables.rhs.resize(hierarchy.gridAtlases[i].numTexels);
		variables.residual.resize(hierarchy.gridAtlases[i].numTexels);
		variables.boundary_rhs.resize(hierarchy.gridAtlases[i].boundaryGlobalIndex.size());
		variables.boundary_value.resize(hierarchy.gridAtlases[i].boundaryGlobalIndex.size());
		variables.variable_boundary_value.resize(hierarchy.gridAtlases[i].boundaryGlobalIndex.size());
	}

	mass_x0.resize(textureNodes.size());
	MultiplyBySystemMatrix_NoReciprocals(deepMassCoefficients, boundaryDeepMassMatrix, boundaryBoundaryMassMatrix, hierarchy.gridAtlases[0].boundaryGlobalIndex, hierarchy.gridAtlases[0].rasterLines, _x0, mass_x0);

	stiffness_x0.resize(textureNodes.size());
	MultiplyBySystemMatrix_NoReciprocals(deepStiffnessCoefficients, boundaryDeepStiffnessMatrix, boundaryBoundaryStiffnessMatrix, hierarchy.gridAtlases[0].boundaryGlobalIndex, hierarchy.gridAtlases[0].rasterLines, _x0, stiffness_x0);

	for (int i = 0; i <_x0.size(); i++){
		multigridFilteringVariables[0].x[i] = _x0[i];
		multigridFilteringVariables[0].rhs[i] = mass_x0[i] + stiffness_x0[i]*sharpnessWeight*stiffnessWeight;
	}

	
	if (UseDirectSolver.set) {
		exactSolution.resize(textureNodes.size());
		ComputeExactSolution(Verbose.set);
	}
	else {
		for (int i = 0; i < updateVCycles; i++) UpdateSolution();
	}

	filteredTexture.resize(width, height);
	for (int i = 0; i < filteredTexture.size(); i++) filteredTexture[i] = Point3D<float>(0.5f, 0.5f, 0.5f);

	if(UseDirectSolver.set) UpdateFilteredTexture(exactSolution);
	else UpdateFilteredTexture(multigridFilteringVariables[0].x);
	return 1;
}

template<class Real>
void TextureFilter<Real>::InitializeVisualization(){


	visualization.textureWidth = textureWidth;
	visualization.textureHeight = textureHeight;

	visualization.colorTextureBuffer = new unsigned char[textureHeight*textureWidth * 3];
	memset(visualization.colorTextureBuffer, 128, textureHeight * textureWidth * 3 * sizeof(unsigned char));



	int tCount = mesh.triangles.size();

	visualization.triangles.resize(tCount);
	visualization.vertices.resize(3 * tCount);
	visualization.colors.resize(3 * tCount, Point3D<double>(0.75, 0.75, 0.75));
	visualization.textureCoordinates.resize(3 * tCount);
	visualization.normals.resize(3 * tCount);


	for (int i = 0; i < tCount; i++) for (int k = 0; k < 3; k++) visualization.triangles[i][k] = 3 * i + k;

	for (int i = 0; i<tCount; i++){
		for (int j = 0; j < 3; j++){
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

	visualization.callBacks.push_back(Visualization::KeyboardCallBack(&visualization, 's', "export texture", "Output Texture", ExportTextureCallBack));
	visualization.callBacks.push_back(Visualization::KeyboardCallBack(&visualization, 'y', "stiffness weight", "Stiffnes Weight", StiffnessWeightCallBack));
	visualization.callBacks.push_back(Visualization::KeyboardCallBack(&visualization, 'g', "sharpness weight", "Sharpness Weight", SharpnessWeightCallBack));

	char* ext = GetFileExtension(InputTexture.value);
	if (!strcasecmp(ext, "png")){
		visualization.textureType = COLOR_TEXTURE;
	}
	else if (!strcasecmp(ext, "normap")){
		visualization.textureType = NORMAL_TEXTURE;
		visualization.normalProgram = new GLSLProgram("../Shaders/normal_texture_vertex.vs", "../Shaders/normal_texture_fragment.fs", "normalTexture");
		visualization.normalProgram->setup();
	}
	else {
		printf("Unexpected texture type! \n");
	}

	visualization.UpdateVertexBuffer();
	visualization.UpdateFaceBuffer();
	visualization.UpdateTextureBuffer(filteredTexture);

	visualization.maskBufferValues = new unsigned char[textureHeight*textureWidth * 3];
	memset(visualization.maskBufferValues, 128, textureHeight * textureWidth * 3 * sizeof(unsigned char));
	for (int i = 0; i < textureNodes.size(); i++) {
			int ci = textureNodes[i].ci;
			int cj = textureNodes[i].cj;
			int offset = 3 * (textureWidth*cj + ci);
			visualization.maskBufferValues[offset + 0] = (unsigned char)(0.8 * 255.0);
			visualization.maskBufferValues[offset + 1] = (unsigned char)(0.8 * 255.0);
			visualization.maskBufferValues[offset + 2] = (unsigned char)(0.8 * 255.0);
	}
	visualization.UpdateMaskTextureBuffer();
}

template<class Real>
int TextureFilter<Real>::Init(){
	levels = std::max<int>(Levels.value,1);
	stiffnessWeight = StiffnessWeight.value;
	sharpnessWeight = SharpnessWeight.value;
	updateVCycles = MultigridUpdateVcycles.value;
	cycleCount = updateVCycles;

	if (!ReadTexturedMesh(mesh, Mesh.value, InputTexture.value)){
		printf("Unable to read mesh data\n");
		return 0;
	}

	textureWidth = mesh.texture.width();
	textureHeight = mesh.texture.height();

	//Define centroid and scale for visualization
	Point3D<double> centroid(0.f, 0.f, 0.f);
	for (int i = 0; i < mesh.vertices.size(); i++) centroid += mesh.vertices[i];
	centroid /= double(mesh.vertices.size());
	double radius = 0.f;
	for (int i = 0; i < mesh.vertices.size(); i++) radius = std::max<double>(radius, Point3D<double>::Length(mesh.vertices[i] - centroid));
	for (int i = 0; i < mesh.vertices.size(); i++) mesh.vertices[i] = (mesh.vertices[i] - centroid) / radius;

	if (1) for (int i = 0; i < mesh.textureCoordinates.size(); i++)mesh.textureCoordinates[i][1] = 1.0 - mesh.textureCoordinates[i][1];

	if (RandomJitter.set) {
		printf("Jittering \n");
		double jitterScale = 1e-3 / double(std::max<int>(textureWidth, textureHeight));
		AddRandomJitter(mesh.textureCoordinates, jitterScale);
	}

	ComputePadding(padding, textureWidth, textureHeight, mesh.textureCoordinates);
	if (padding.nonTrivial) {
		PaddTextureCoordinates(padding, textureWidth, textureHeight, mesh.textureCoordinates);
		PaddImage(padding, mesh.texture);
		textureWidth += (padding.left + padding.right);
		textureHeight += (padding.bottom + padding.top);
	}

	printf("Resolution %d x %d \n", textureWidth, textureHeight);
	if (!InitializeSystem(textureWidth, textureHeight)){
		printf("Unable to initialize system\n");
		return 0;
	}

	//Assign position to exterior nodes using baricentric-exponential map
	FEM::RiemannianMesh< double > * Rmesh = new FEM::RiemannianMesh< double >(GetPointer(mesh.triangles), mesh.triangles.size());
	Rmesh->setMetricFromEmbedding(GetPointer(mesh.vertices));
	Rmesh->makeUnitArea();
	Pointer(FEM::CoordinateXForm< double >) xForms = Rmesh->getCoordinateXForms();

	for (int i = 0; i<textureNodes.size(); i++){
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
	//

	textureNodePositions.resize(textureNodes.size());
	for (int i = 0; i < textureNodePositions.size(); i++){
		Point2D<double> barincetricCoords = textureNodes[i].baricentricCoords;
		int tId = textureNodes[i].tId;
		Point3D<float> surfacePosition = mesh.vertices[mesh.triangles[tId][0]] * (1.0 - barincetricCoords[0] - barincetricCoords[1]) +
										 mesh.vertices[mesh.triangles[tId][1]] * barincetricCoords[0] +
										 mesh.vertices[mesh.triangles[tId][2]] * barincetricCoords[1];
		textureNodePositions[i] = surfacePosition;
	}

	unifromTexelSharpenningMask.resize(textureNodes.size(), 0.5);

	for (int c = 0; c < 3; c++) texelStiffness[c].resize(textureNodes.size());

	cellSharpenningMask.resize(cellIndices.size(), 1);
	uniformCellSharpenningMask.resize(cellIndices.size(), 0.5);

	cellCenterPositions.resize(cellIndices.size());
	for (int i = 0; i < cellIndices.size(); i++){
		cellCenterPositions[i] = (textureNodePositions[cellIndices[i][0]] +
								  textureNodePositions[cellIndices[i][1]] +
								  textureNodePositions[cellIndices[i][2]] +
								  textureNodePositions[cellIndices[i][3]]) / 4.0;
	}

	if (1) {
		int multiChartTexelCount = 0;
		Image<int> texelId;
		texelId.resize(textureWidth, textureHeight);
		for (int i = 0; i < texelId.size(); i++)texelId[i] = -1;
		for (int i = 0; i < textureNodes.size(); i++) {
			int ci = textureNodes[i].ci;
			int cj = textureNodes[i].cj;
			if (texelId(ci, cj) != -1) {
				if(0) printf("WARNING: Texel (%d %d) belong to multiple charts! \n", ci, cj);
				multiChartTexelCount++;
			}
			texelId(ci, cj) = i;
		}
		if (multiChartTexelCount)printf("WARNING: %d texels belong to multiple charts! \n", multiChartTexelCount);
	}
	return 1;
}

template<class Real>
int _main(int argc, char* argv[])
{
	if (!TextureFilter<Real>::Init()) return 0;

	if (!OutputTexture.set){
		glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
		TextureFilter<Real>::visualization.displayMode = DisplayMode.value;
		if (DisplayMode.value == ONE_REGION_DISPLAY) {
			TextureFilter<Real>::visualization.screenWidth = 800;
			TextureFilter<Real>::visualization.screenHeight = 800;
		}
		else if (DisplayMode.value == TWO_REGION_DISPLAY) {
			TextureFilter<Real>::visualization.screenWidth = 1600;
			TextureFilter<Real>::visualization.screenHeight = 800;
		}
		else if (DisplayMode.value == THREE_REGION_DISPLAY) {
			TextureFilter<Real>::visualization.screenWidth = 1200;
			TextureFilter<Real>::visualization.screenHeight = 800;
		}
		else if (DisplayMode.value == FOUR_REGION_DISPLAY) {
			TextureFilter<Real>::visualization.screenWidth = 1500;
			TextureFilter<Real>::visualization.screenHeight = 600;
		}
		TextureFilter<Real>::visualization.UpdateMainFrameSize();
		glutInitWindowSize(TextureFilter<Real>::visualization.screenWidth, TextureFilter<Real>::visualization.screenHeight);
		glutInit(&argc, argv);
		char windowName[1024];
		sprintf(windowName, "Texture Filtering");
		glutCreateWindow(windowName);
		if (glewInit() != GLEW_OK) fprintf(stderr, "[ERROR] glewInit failed\n"), exit(0);
		glutDisplayFunc(TextureFilter<Real>::Display);
		glutReshapeFunc(TextureFilter<Real>::Reshape);
		glutMouseFunc(TextureFilter<Real>::MouseFunc);
		glutMotionFunc(TextureFilter<Real>::MotionFunc);
		glutKeyboardFunc(TextureFilter<Real>::KeyboardFunc);
		if (!UseDirectSolver.set) glutIdleFunc(TextureFilter<Real>::Idle);
		if (CameraConfig.set) TextureFilter<Real>::visualization.ReadSceneConfigurationCallBack(&TextureFilter<Real>::visualization, CameraConfig.value);
		TextureFilter<Real>::InitializeVisualization();
		TextureFilter<Real>::visualization.showSlideBar = true;
		glutMainLoop(); 
	}
	else{
		TextureFilter<Real>::ExportTextureCallBack(&TextureFilter<Real>::visualization,OutputTexture.value);
	}

	return 0;
}

int main(int argc, char* argv[])
{
	cmdLineParse(argc - 1, argv + 1, params);
	if (!Mesh.set || !InputTexture.set){ ShowUsage(argv[0]); return EXIT_FAILURE; }
	omp_set_num_threads(Threads.value);
	printf("Number of threads %d. Number of processors %d. \n", omp_get_max_threads(), omp_get_num_procs());
	_main<REAL_TYPE>(argc, argv); 
	return 0;
}
