#ifndef MATH0X_FUNC_RIESZ_H
#define MATH0X_FUNC_RIESZ_H

#include <math0x/euclid.h>

namespace math0x { 
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

	      // hold a reference 
				const E& Mx = metric(x);
      
				euclid::dual<E> res = dual.zero();
      
				euclid::range< euclid::dual<E> > rres( res );
				euclid::range< E > rMx( Mx );
				
				for( ;!rMx.empty(); rMx.pop(), rres.pop() ) {
					const_cast< euclid::field<E>& >(rres.front()) = rMx.front();
				}
				
				return res;
      }
    

      // TODO push/pull ?
    };


  }

}
#endif
