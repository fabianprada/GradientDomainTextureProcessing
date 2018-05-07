#include <Eigen/Core>

#if USE_FLOAT
#define REAL_TYPE float
typedef Eigen::VectorXf EigenVector;
#include <vectorclass/vectorclass.h>

#define VECTOR_TYPE Vec4f
inline Vec4f ZeroData() {
	return Vec4f(0);
}
inline Vec4f SetData(float a, float b, float c) {
	return Vec4f(a, b, c, 0);
}
inline float DotData(const Vec4f & v1, const Vec4f & v2) {
	return horizontal_add(v1*v2);
}

#else 
#define REAL_TYPE double
typedef Eigen::VectorXd EigenVector;

#include <vectorclass/vectorclass.h>
#define VECTOR_TYPE Vec4d
inline Vec4d ZeroData() {
	return Vec4d(0);
}
inline Vec4d SetData(double a, double b, double c) {
	return Vec4d(a, b, c, 0);
}
inline double DotData(const Vec4d & v1, const Vec4d & v2) {
	return horizontal_add(v1*v2);
}

#define DOUBLE_VECTOR_TYPE Vec4d
inline Vec4d ZeroDataDouble() {
	return Vec4d(0);
}

inline double DotDataDouble(const Vec4d & v1, const Vec4d & v2) {
	return horizontal_add(v1*v2);
}
#endif
