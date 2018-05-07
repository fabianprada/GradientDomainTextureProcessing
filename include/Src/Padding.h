#pragma once

class Padding {
public:
	Padding() {
		left = bottom = right = top = 0;
		nonTrivial = false;
	}
	int left;
	int bottom;
	int right;
	int top;
	bool nonTrivial;
};

void ComputePadding(Padding & padding, const int & width, const int & height, std::vector<Point2D<double>>  & textureCoordinates) {
	Point2D<double> pixMinCorner(0.5 / double(width), 0.5 / double(height));
	Point2D<double> pixMaxCorner((double(width) - 0.5) / double(width), (double(height) - 0.5) / double(height));

	Point2D<double> texMinCorner = textureCoordinates[0];
	Point2D<double> texMaxCorner = textureCoordinates[0];
	for (int i = 0; i < textureCoordinates.size(); i++) {
		for (int c = 0; c < 2;c++) texMinCorner[c] = std::min<double>(texMinCorner[c], textureCoordinates[i][c]);
		for (int c = 0; c < 2;c++) texMaxCorner[c] = std::max<double>(texMaxCorner[c], textureCoordinates[i][c]);
	}

	padding.left = texMinCorner[0] < pixMinCorner[0] ? ceil((pixMinCorner[0] - texMinCorner[0]) * double(width)) : 0;
	padding.bottom = texMinCorner[1] < pixMinCorner[1] ? ceil((pixMinCorner[1] - texMinCorner[1]) * double(height)) : 0;

	padding.right = texMaxCorner[0] > pixMaxCorner[0] ? ceil((texMaxCorner[0] - pixMaxCorner[0]) * double(width)) : 0;
	padding.top = texMaxCorner[1] > pixMaxCorner[1] ? ceil((texMaxCorner[1] - pixMaxCorner[1]) * double(height)) : 0;

	//Make image dimensions are multiples of 8 (Hardware texture mapping seems to fail if not)
	int newWidth = width + padding.left + padding.right;
	int newHeight = height + padding.bottom + padding.top;

	int paddedWidth = 8 * (((newWidth - 1) / 8) + 1);
	int paddedHeight = 8 * (((newHeight - 1) / 8) + 1);
	padding.left += (paddedWidth - newWidth);
	padding.bottom += (paddedHeight - newHeight);

	if (padding.left || padding.right || padding.bottom || padding.top) padding.nonTrivial = true;
	if(padding.nonTrivial) printf("Padding applied : Left %d. Right %d. Bottom %d. Top %d.\n", padding.left, padding.right, padding.bottom, padding.top);
	else printf("No padding required! \n");

}

template<class DataType>
void PaddImage(Padding & padding, Image<DataType> & im){
	int newWidth = im.width() + padding.left + padding.right;
	int newHeight = im.height() + padding.bottom + padding.top;

	Image<Point3D<float>> newIm;
	newIm.resize(newWidth, newHeight);
	for (int i = 0; i < newWidth; i++) for (int j = 0; j < newHeight; j++) {
		int ni = std::min<int>(std::max<int>(0, i - padding.left), im.width() - 1);
		int nj = std::min<int>(std::max<int>(0, j - padding.bottom), im.height() - 1);
		newIm(i, j) = im(ni, nj);
	}
	im = newIm;
}

template<class DataType>
void UnpaddImage(Padding & padding, Image<DataType> & im){
	int outputWidth = im.width() - padding.left - padding.right;
	int outputHeight = im.height() - padding.bottom - padding.top;
	Image<Point3D<float>> newIm;
	newIm.resize(outputWidth, outputHeight);
	for (int i = 0; i < outputWidth; i++)for (int j = 0; j < outputHeight; j++)newIm(i, j) = im(padding.left + i, padding.bottom + j);
	im = newIm;
}


void PaddTextureCoordinates(Padding & padding, int width, int height, std::vector<Point2D<double>>  & textureCoordinates){
	int newWidth = width + padding.left + padding.right;
	int newHeight = height + padding.bottom + padding.top;

	for (int i = 0; i < textureCoordinates.size(); i++) {
		textureCoordinates[i][0] = (textureCoordinates[i][0] * double(width) + double(padding.left)) / double(newWidth);
		textureCoordinates[i][1] = (textureCoordinates[i][1] * double(height) + double(padding.bottom)) / double(newHeight);
	}
}

void UnpaddTextureCoordinates(Padding & padding, int width, int height, std::vector<Point2D<double>>  & textureCoordinates) {
	int newWidth = width + padding.left + padding.right;
	int newHeight = height + padding.bottom + padding.top;

	for (int i = 0; i < textureCoordinates.size(); i++) {
		textureCoordinates[i][0] = (textureCoordinates[i][0] * double(newWidth) - double(padding.left)) / double(width);
		textureCoordinates[i][1] = (textureCoordinates[i][1] * double(newHeight) - double(padding.bottom)) / double(height);
	}
}