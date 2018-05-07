#ifndef MAPLOOP_INCLUDED
#define MAPLOOP_INCLUDED

#include <unordered_map>
#include <unordered_set>


int LoopVertices(std::unordered_map<unsigned long long, unsigned long long> & forwardMap, std::vector<std::vector<unsigned long long>> & loopVertices){
	loopVertices.clear();
	std::unordered_set<unsigned long long> alreadyVisitedVertex;
	unsigned int loopCounter = 0;
	unsigned int edgeCounter = 0;
	for (auto iter = forwardMap.begin(); iter != forwardMap.end(); iter++){
		unsigned long long sourceVertex = (*iter).first;
		if (alreadyVisitedVertex.find(sourceVertex) == alreadyVisitedVertex.end()){
			std::vector< unsigned long long> currentLoop;
			unsigned long long currentVertex = sourceVertex;
			bool terminate = false;
			unsigned int startEdgeCounter = edgeCounter;
			do{
				if (alreadyVisitedVertex.find(currentVertex) == alreadyVisitedVertex.end()){
					alreadyVisitedVertex.insert(currentVertex);
					currentLoop.push_back(currentVertex);
					auto mappedVertex = forwardMap.find(currentVertex);
					if (mappedVertex != forwardMap.end()){
						unsigned long long nextVertex = (*mappedVertex).second;
						edgeCounter++;
						currentVertex = nextVertex;
					}
					else{
						printf("Vertex to dead end! \n");
						return 0;
					}
				}
				else{
					if (currentVertex != sourceVertex){
						printf("Non simple loop! \n");
						printf("Node %llu \n", currentVertex);
						return 0;
					}
					terminate = true;
				}
			} while (!terminate);
			loopCounter++;
			loopVertices.push_back(currentLoop);
		}
	}
	return 1;
}

int ListVerticesSimpleLoop(std::unordered_map<int, int> & forwardMap, std::vector<int> & vertexList) {
	vertexList.clear();
	std::unordered_set<int> alreadyVisitedVertex;
	int loopCounter = 0;

	int sourceVertex = (*forwardMap.begin()).first;
	int currentVertex = sourceVertex;
	bool terminate = false;
	do {
		if (alreadyVisitedVertex.find(currentVertex) == alreadyVisitedVertex.end()) {
			alreadyVisitedVertex.insert(currentVertex);
			vertexList.push_back(currentVertex);
			auto mappedVertex = forwardMap.find(currentVertex);
			if (mappedVertex != forwardMap.end()) {
				unsigned long long nextVertex = (*mappedVertex).second;
				currentVertex = nextVertex;
			}
			else {
				printf("Vertex to dead end! \n");
				return 0;
			}
		}
		else {
			if (currentVertex != sourceVertex) {
				printf("Non simple loop! \n");
				return 0;
			}
			terminate = true;
		}
	} while (!terminate);

	return 1;
}

#endif //MAPLOOP_INCLUDED