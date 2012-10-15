#ifndef MATH0X_FUNC_MINUS_H
#define MATH0X_FUNC_MINUS_H

#include <math0x/euclid.h>

namespace math0x { 
	namespace func {

		template<class E>
		struct minus {
			typedef minus base;
    
			euclid::space<E> space;

			E operator()(const E& x) const {
				return space.minus(x);
			}

    
			struct push;
			struct pull;
    
		};
  

		template<class E>
		struct minus<E>::push : minus {
      
			push(const minus& of, const E& ) : push::base(of) { }
    
		};

		template<class E>
		struct minus<E>::pull : minus< euclid::dual<E> > {
    
			pull(const minus& of, const E& ) : pull::base{ *of.space } {  }
    
		};
  
		template<class E>
		minus<E> make_minus(const euclid::space<E>& space) { return {space}; }
		
	}

}

#endif
