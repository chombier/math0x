#ifndef MATH0X_MATRIX_H
#define MATH0X_MATRIX_H

#include <math0x/eigen.h>

namespace math0x { 

  template<class U = RR, int Rows = Eigen::Dynamic, int Cols = Eigen::Dynamic> 
  using matrix = Eigen::Matrix<U, Rows, Cols, 
			       Eigen::Matrix<U, Rows, Cols>::Options,
			       Eigen::Matrix<U, Rows, Cols>::MaxRowsAtCompileTime,
			       Eigen::Matrix<U, Rows, Cols>::MaxColsAtCompileTime>;

  typedef matrix<RR, 2, 2> mat22;
  typedef matrix<RR, 3, 3> mat33;
  typedef matrix<RR, 4, 4> mat44;
  typedef matrix<RR, 6, 6> mat66;

	typedef matrix<RR> mat;
}

#endif
