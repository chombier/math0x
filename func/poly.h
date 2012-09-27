#ifndef GROUP_FUNC_POLY_H
#define GROUP_FUNC_POLY_H

#include <group/types.h>

namespace func {

  // a monomial of degree n: x -> x^n
  template<class U = real>
  struct poly {
    typedef poly base;
    
    NN degree;
    
    // fast exponentiation
    U operator()(const U& x) const {
      if( degree == 0 ) return 1.0;
      if( degree == 1 ) return x;
	
      U res_sqrt = poly{degree / 2}(x);
      
      if( degree % 2 == 0 ) {
	return res_sqrt * res_sqrt;
      } else {
	return x * res_sqrt * res_sqrt;
      }
	
    }

    struct push {
      U value;

      push(const poly& of, const U& at) 
	: value( of.degree * poly{of.degree - 1}(at) ) {

      }

      U operator()(const U& dx) const {
	return value * dx;
      }
      
    };

    typedef push pull;
    
  };

}


#endif
