#ifndef CONSTRAINED_TRIANGULATION_INCLUDED
#define CONSTRAINED_TRIANGULATION_INCLUDED

#include<Misha/Geometry.h>

#undef REAL
#undef dest
#define ANSI_DECLARATORS
#define TRILIBRARY
#define NO_TIMER

#include "triangle.c"

void TriangulatePolygon(const std::vector<Point2D<double>> & vertices, std::vector<TriangleIndex> & outputTriangles){

	struct triangulateio in, out;

	/* Define input points. */

	in.numberofpoints = vertices.size();

	in.pointlist = (REAL *)malloc(in.numberofpoints * 2 * sizeof(REAL));
	in.pointmarkerlist = (int *)malloc(in.numberofpoints * sizeof(int));
	for (int i = 0; i < in.numberofpoints; i++) {
		in.pointlist[2 * i] = vertices[i][0];
		in.pointlist[2 * i + 1] = vertices[i][1];
		in.pointmarkerlist[i] = 1; //Check boundary markers documentation
	}

	in.numberofsegments = vertices.size();
	in.segmentlist = (int *)malloc(in.numberofsegments * 2 * sizeof(int));
	in.segmentmarkerlist = (int *)malloc(in.numberofsegments * sizeof(int));
	for (int i = 0; i < in.numberofsegments; i++) {
		in.segmentlist[2 * i] = i + 1;
		in.segmentlist[2 * i + 1] = i < (in.numberofsegments - 1) ? (i + 2) : 1;
		in.segmentmarkerlist[i] = 1;
	}
	in.numberofholes = 0;
	in.numberofregions = 0;
	in.numberofpointattributes = 0;

	out.pointlist = (REAL *)NULL;
	out.trianglelist = (int *)NULL;
	out.segmentlist = (int *)NULL;
	out.pointmarkerlist = (int *)NULL;
	out.triangleattributelist = (REAL *)NULL;
	out.segmentmarkerlist = (int *)NULL;
	/* Refine the triangulation according to the attached */
	/*   triangle area constraints.                       */

	triangulate("pQ", &in, &out, (struct triangulateio *) NULL);

	outputTriangles.resize(out.numberoftriangles);
	for (int i = 0; i < out.numberoftriangles; i++) outputTriangles[i] = TriangleIndex(out.trianglelist[3 * i] - 1, out.trianglelist[3 * i + 1] - 1, out.trianglelist[3 * i + 2] - 1);

	free(in.pointlist);
	free(in.pointmarkerlist);
	free(in.segmentlist);
	free(in.segmentmarkerlist);


	free(out.pointlist);
	free(out.trianglelist);
	free(out.segmentlist);
	free(out.pointmarkerlist);
	free(out.triangleattributelist);
	free(out.segmentmarkerlist);
}


#undef REAL
#undef dest
#undef ANSI_DECLARATORS
#undef TRILIBRARY
#undef NO_TIMER

#endif //CONSTRAINED_TRIANGULATION_INCLUDED
