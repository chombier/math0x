#ifndef MATH0X_QUATERNION_H
#define MATH0X_QUATERNION_H

#include <Eigen/Geometry>

namespace math0x { 
  template<class U = real> 
  using quaternion = Eigen::Quaternion<U, Eigen::Quaternion<U>::Coefficients::Options >;

  typedef quaternion<RR> HH;

}
#endif
