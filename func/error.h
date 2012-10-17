#ifndef MATH0X_FUNC_ERROR_H
#define MATH0X_FUNC_ERROR_H

#include <math0x/types.h>
#include <math0x/lie.h>

namespace math0x { 
	namespace func {

		// throw an error
		template<class Domain, class Range, class What >
		struct error {
			typedef error base;
    
			What what;
    
			Range operator()(const Domain& ) const { 
				throw what;
			}
    
    
			struct push : error< lie::algebra<Domain>, lie::algebra<Range>, What > {
      
				push(const error& of, const Domain& ) : push::base{of.what} { }
      
			};


			struct pull : error< lie::coalgebra<Range>, lie::coalgebra<Domain>, What > {
				pull(const error& of, const Domain& ) : pull::base{of.what} { }
			};

		};


	};

}
#endif
