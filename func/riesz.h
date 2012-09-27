#ifndef GROUP_FUNC_RIESZ_H
#define GROUP_FUNC_RIESZ_H

#include <group/euclid.h>

namespace func {

  // Riesz representation, given metric M: x -> (Mx)^T
  template<class E, class Metric = id<E> > 
  struct riesz {
    
    euclid::space<E> primal;
    euclid::space< euclid::dual<E> > dual;
    
    Metric metric;
    
    riesz(const euclid::space<E>& primal = {},
	  const Metric& metric = {} )
      : primal(primal),
	dual( *primal ),
	metric( metric ) {

    }
    
    euclid::dual<E> operator()(const E& x) const {

      const E& Mx = metric(x);
      
      euclid::dual<E> res = dual.zero();
      
      NN i = 0;
      dual.each(res, [&](euclid::field<E>& res_i) {
	  res_i = primal.coord(i, Mx);
	  ++i;
	});
      
      return res;
    }
    

    // TODO push/pull ?
  };


  // TODO add transpose function for canonical dual vector ?

}


#endif
