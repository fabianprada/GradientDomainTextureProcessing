#pragma once

int InitializeInteriorTexelToCellLines(std::vector<InteriorTexelToCellLine> & interiorTexeltoCellLine, const GridAtlas & gridAtlas) {
	const std::vector<RasterLine> & rasterLines = gridAtlas.rasterLines;
	const std::vector<GridNodeInfo> & nodeInfo = gridAtlas.nodeInfo;
	const std::vector<GridChart> & gridCharts = gridAtlas.gridCharts;
	interiorTexeltoCellLine.resize(rasterLines.size());
	for (int i = 0; i < rasterLines.size(); i++) {
		int interiorTexelStart = rasterLines[i].lineStartIndex;
		int ci = nodeInfo[interiorTexelStart].ci;
		int cj = nodeInfo[interiorTexelStart].cj;
		int chartId = nodeInfo[interiorTexelStart].chartId;

		interiorTexeltoCellLine[i].texelStartIndex = rasterLines[i].lineStartIndex;
		interiorTexeltoCellLine[i].texelEndIndex = rasterLines[i].lineEndIndex;
		interiorTexeltoCellLine[i].coeffOffset = rasterLines[i].coeffStartIndex;

		if (gridCharts[chartId].cellType(ci - 1, cj - 1) != 1) {
			printf("ERROR: Non interior cell! \n");
		}
		interiorTexeltoCellLine[i].previousCellStartIndex = gridCharts[chartId].localCellIndex(ci - 1, cj - 1) + gridCharts[chartId].globalIndexCellOffset;

		if (gridCharts[chartId].cellType(ci - 1, cj) != 1) {
			printf("ERROR: Non interior cell! \n");
		}
		interiorTexeltoCellLine[i].nextCellStartIndex = gridCharts[chartId].localCellIndex(ci - 1, cj) + gridCharts[chartId].globalIndexCellOffset;
	}
	return 1;
}
