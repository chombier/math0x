#ifndef MATH0X_TYPES_H
#define MATH0X_TYPES_H

#include <math0x/traits.h>

namespace math0x {
  // geometry
  namespace euclid {
    template<class E> class space;
  }

  namespace lie {
    template<class G> class group;
  }

  // basic types
  typedef int integer;
  typedef unsigned int natural;
  typedef double real;

  // // vector/matrices
  // static constexpr int dynamic_size = -1;

  // template<class U = real, int M = dynamic_size> class vector;
  // template<class U = real, int M = dynamic_size> class linear_form;

  // template<class U = real, int M = dynamic_size, int N = dynamic_size> class matrix;

  // some functions
  namespace func {
  
    template<class F> struct default_push;
    template<class F> struct default_pull;

    template<class G> struct id;

    template<class E> struct minus;
  
    template<class Outer, class Inner> struct comp;

    template<class Domain, class Range, class Error> struct error;

    template<class ... > struct tuple;

    template<class G> struct apply;  

    template<class, int I = -1> struct get;
    template<class, int I = -1> struct part;
  }


  // matrix groups
  // template<class U = real> struct quaternion;


  // template<class U = real> struct dual_quaternion;

  // // handy typedefs
  // typedef vector<real> vec; 
  // typedef linear_form<real> form;
  // typedef matrix<real> mat;

  // typedef vector<real, 2> vec2;
  // typedef vector<real, 3> vec3;

  // typedef linear_form<real, 2> form2;
  // typedef linear_form<real, 3> form3;

  // typedef matrix<real, 2, 2> mat22;
  // typedef matrix<real, 3, 3> mat33;
  // typedef matrix<real, 3, 2> mat32;
  // typedef matrix<real, 2, 3> mat23;
  // typedef matrix<real, 4, 4> mat44;
  // typedef matrix<real, 4, 3> mat43;
  // typedef matrix<real, 6, 6> mat66;
  // typedef matrix<real, 3, 6> mat36;

  // mathematical notations
  typedef natural NN;
  typedef integer ZZ;
  typedef real RR;

  // special orthogonal
  template<int, class = RR> class SO;
	typedef SO<3> SO3;
	
  template<int, class = RR> struct so;
	
  // special euclidean
  template<int, class = RR> struct SE;
	typedef SE<3> SE3;
	
	template<int, class = RR> class se;

  // special linear
  template<int, class = RR> class SL;
  template<int, class = RR> class sl;
	
}

#endif
