#ifndef MATH0X_VECTOR_H
#define MATH0X_VECTOR_H

#include <math0x/eigen.h>

namespace math0x { 
  template<class U = RR, int Rows = Eigen::Dynamic> 
  using vector = Eigen::Matrix<U, Rows, 1, 
			       Eigen::Matrix<U, Rows, 1>::Options,
			       Eigen::Matrix<U, Rows, 1>::MaxRowsAtCompileTime,
			       Eigen::Matrix<U, Rows, 1>::MaxColsAtCompileTime>;

	// shorthand
  typedef vector<RR, 2> vec2;
  typedef vector<RR, 3> vec3;
  typedef vector<RR, 4> vec4;
  typedef vector<RR, 6> vec6;

	typedef vector<RR> vec;
}
#endif
