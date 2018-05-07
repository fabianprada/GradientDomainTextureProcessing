#ifndef HSV_INCLUDED
#define HSV_INCLUDED

#include <Misha/Geometry.h>
// h [0,1) ,s [0,1], v [0,1]
// Hue is mapped to [0,240) deg to avoid circularity
Point3D<float> HSV2RGB(const double h, const double s, const double v){
	const double c = s*v;
	// use next line for full [0,360) hue
	const double _h = h * 6.0;
	// use next line for partial [0,240) hue
	//const double _h = h * 4.0;
	const double _h_mod_2 = _h - floor(_h / 2.0)*2.0;
	double x = c*(1 - fabs(double(_h_mod_2 - 1)));
	double r, g, b;
	r = g = b = 0.0;
	if (_h <= 1){
		r = c;
		g = x;
	}
	else if (_h <= 2){
		r = x;
		g = c;
	}
	else if (_h <= 3){
		g = c;
		b = x;
	}
	else if (_h <= 4){
		g = x;
		b = c;
	}
	else if (_h <= 5){
		r = x;
		b = c;
	}
	else{
		r = c;
		b = x;
	}
	double m = v - c;
	return Point3D<float>(r + m, g + m, b + m);
}

#endif //HSV_INCLUDED
