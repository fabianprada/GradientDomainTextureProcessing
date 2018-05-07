#ifndef CHART_DECOMPOSITION_INCLUDE
#define CHART_DECOMPOSITION_INCLUDE
#include <Src/SimpleMesh.h>
#include <queue>

int AddComponent(std::vector<int> & vertexComponent, int vIndex, int currentComponent, const std::vector<std::vector<int>> & neighbours){
	vertexComponent[vIndex] = currentComponent;
	std::queue<int> visitingQueue;
	visitingQueue.push(vIndex);

	while (!visitingQueue.empty()){
		int currentVertex = visitingQueue.front();
		visitingQueue.pop();
		const std::vector<int> & vertexNeighbours = neighbours[currentVertex];
		for (int i = 0; i < vertexNeighbours.size(); i++){
			if (vertexComponent[vertexNeighbours[i]] == -1) {
				vertexComponent[vertexNeighbours[i]] = currentComponent;
				visitingQueue.push(vertexNeighbours[i]);
			}
			else if (vertexComponent[vertexNeighbours[i]] == currentComponent){}
			else{
				printf("Unexpected Condition On A Connected Component. Expected %d. Obtained %d\n", currentComponent, vertexComponent[vertexNeighbours[i]]);
				return 0;
			}
		}
	}
	return 1;
}

int InitializeTriangleChartIndexing(const TexturedMesh & mesh, std::vector<int> & chartIndex, int & numCharts) {
	std::unordered_map<unsigned long long, int> edgeIndex;
	for (int i = 0; i < mesh.triangles.size(); i++) {
		for (int k = 0; k < 3; k++) {
			unsigned long long  edgeKey = SetMeshEdgeKey(mesh.triangles[i][k], mesh.triangles[i][(k + 1) % 3]);
			if (edgeIndex.find(edgeKey) == edgeIndex.end()) {
				edgeIndex[edgeKey] = 3 * i + k;
			}
			else {
				printf("Non manifold mesh!! \n");
				return 0;
			}
		}
	}

	std::vector<std::vector<int>> neighbours(mesh.triangles.size());
	for (int i = 0; i < mesh.triangles.size(); i++) {
		for (int k = 0; k < 3; k++) {
			unsigned long long edgeKey = SetMeshEdgeKey(mesh.triangles[i][(k + 1) % 3], mesh.triangles[i][k]);
			if (edgeIndex.find(edgeKey) != edgeIndex.end()) {
				int tIndex = edgeIndex[edgeKey] / 3;
				int kIndex = edgeIndex[edgeKey] % 3;
				if (mesh.textureCoordinates[3 * i + ((k + 1) % 3)][0] == mesh.textureCoordinates[3 * tIndex + kIndex][0] &&
					mesh.textureCoordinates[3 * i + ((k + 1) % 3)][1] == mesh.textureCoordinates[3 * tIndex + kIndex][1] &&
					mesh.textureCoordinates[3 * i + k][0] == mesh.textureCoordinates[3 * tIndex + ((kIndex + 1) % 3)][0] &&
					mesh.textureCoordinates[3 * i + k][1] == mesh.textureCoordinates[3 * tIndex + ((kIndex + 1) % 3)][1]
					) {
					neighbours[i].push_back(tIndex);
				}
			}
		}
	}
	chartIndex.clear();
	chartIndex.resize(mesh.triangles.size(), -1);
	int currentComponent = -1;
	for (int v = 0; v < mesh.triangles.size(); v++) {
		if (chartIndex[v] == -1) {
			currentComponent++;
			AddComponent(chartIndex, v, currentComponent, neighbours);
		}
	}
	numCharts = currentComponent + 1;
	return 1;
}


#endif //CHART_DECOMPOSITION_INCLUDE