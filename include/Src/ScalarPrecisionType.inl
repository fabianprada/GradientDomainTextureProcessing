#include <Eigen/Core>

#if USE_FLOAT
#define REAL_TYPE float
typedef Eigen::VectorXf EigenVector;
inline float ZeroData() {
	return 0.f;
}
inline float DotData(const float & f1, const float & f2) {
	return f1*f2;
}
#else 
#define REAL_TYPE double
typedef Eigen::VectorXd EigenVector;
inline double ZeroData() {
	return 0;
}
inline double DotData(const double & f1, const double & f2) {
	return f1*f2;
}
#endif