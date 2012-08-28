#ifndef GROUP_QUATERNION_H
#define GROUP_QUATERNION_H

#include <Eigen/Geometry>

template<class U = real> 
using quaternion = Eigen::Quaternion<U, Eigen::Quaternion<U>::Coefficients::Options >;

typedef quaternion<RR> HH;


#endif
