#ifndef GROUP_VECTOR_H
#define GROUP_VECTOR_H

#include <group/eigen.h>

template<class U = real, int Rows = Eigen::Dynamic> 
using vector = Eigen::Matrix<U, Rows, 1, 
			     Eigen::Matrix<U, Rows, 1>::Options,
			     Eigen::Matrix<U, Rows, 1>::MaxRowsAtCompileTime,
			     Eigen::Matrix<U, Rows, 1>::MaxColsAtCompileTime>;

typedef vector<real, 2> vec2;
typedef vector<real, 3> vec3;

#endif
