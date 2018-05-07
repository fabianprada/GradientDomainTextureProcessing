#pragma once
#include <Src/TexturedMeshVisualization.h>

enum {
	COLOR_TEXTURE,
	NORMAL_TEXTURE,
	TEXTURE_TYPE_COUNT
};

class TextureFilteringVisualization : public TexturedMeshVisualization{
public:
	TextureFilteringVisualization();
	int textureType;
	double slideBarCursorPosition = 0.5;
	double slideBarCursorOldPosition = 0.5;
	bool isBrushActive;
	bool isSlideBarActive;
	int diskX, diskY;
	bool showDisk;
	bool showSlideBar;
	GLuint maskTextureBuffer = 0;
	unsigned char * maskBufferValues;
	unsigned char * colorTextureBuffer;
	int textureWidth;
	int textureHeight;
	void UpdateTextureBuffer(const Image<Point3D<float>> & image);
	void UpdateColorTextureBuffer();
	void UpdateMaskTextureBuffer();

	void display(void);
	void LoadGeometryData();
	void DrawSlideBar();
};

#include "TextureFilteringVisualization.inl"