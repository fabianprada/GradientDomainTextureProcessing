# GradientDomainTextureProcessing

-----------------------------------------------------------------
INSTALLATION

This code was compiled and tested on Windows using Visual Studio 2017, and on Linux using GCC 5.5. 

To compile on Windows copy /include and /lib from 4Windows.zip to the main directory, and build the applications from the provided solution file. Add the content of /dll to /x64/Release to run the applications.

To compile in Linux, install PNG, GLEW and GLUT libraries on your system, and build the binaries from the Makefile provided for each application. 
-----------------------------------------------------------------

SUPPORTED FILE FORMATS

Meshes. The input. ply mesh must store uv-coordiante per triangle corner . Meshes that store uv-coordinates per vertex and duplicate vertices at the seam are not supported.

Textures. For color textures the supported format is .png. For normal textures we use a binary file with the extension. normap. This stores a two-dimensional array where each normal correspond to a triplet of double precision numbers.

Tangent Vector Fields. Vector fields are stored in binary files with extension. vf. This stores a one-dimensional array where each tangent vector corresponds to a pair of double precision numbers.
-----------------------------------------------------------------
 
APPLICATIONS


(1) TextureFiltering: Smooth or sharpen texture detail using local and global edition tools. Run the application as:

TextureFiltering --mesh mesh.ply --inTexture texture.png/texture.normap


(2) Geodesics: Compute the geodesic distance to a source point on the surface. Run the application as:

Geodesics --mesh mesh.ply

(2) Line Integral Convolution: Visualize a tangential vector field by computing the line integral convolution over a random texture. Run the application as:

LineIntegralConvolution --mesh mesh.ply --major(default)/--minor/--vf vectorField.vf


Run each application with no parameters to print the command line options.
-----------------------------------------------------------------
 
DEMO

View a demo of this applications from: 

http://www.cs.jhu.edu/~fpradan1/code/




