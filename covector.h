#ifndef MATH0X_COVECTOR_H
#define MATH0X_COVECTOR_H

#include <Eigen/Core>

namespace math0x { 

  template<class U = real, int Cols = Eigen::Dynamic> 
  using covector = Eigen::Matrix<U, 1, Cols, 
				 Eigen::Matrix<U, 1, Cols>::Options,
				 Eigen::Matrix<U, 1, Cols>::MaxRowsAtCompileTime,
				 Eigen::Matrix<U, 1, Cols>::MaxColsAtCompileTime>;

  typedef covector<real, 2> form2;
  typedef covector<real, 3> form3;

}
#endif
