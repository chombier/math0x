#ifndef MATH0X_FUNC_DEFAULT_H
#define MATH0X_FUNC_DEFAULT_H

#include <math0x/func/error.h>
#include <stdexcept>

#include <math0x/lie.h>

namespace math0x {
	namespace func {

		// default pushforward/pullback
		template<class F>
		struct default_push : func::error< lie::algebra< func::domain<F> >,
		                                   lie::algebra< func::range<F> >,
		                                   std::logic_error > {
			default_push(const F&, const func::domain<F>& )
				: default_push::base{ std::logic_error( meta::name<F>() + " has no pushforward") }  {
					
			}      
      
		};

		template<class F>
		struct default_pull : func::error< lie::coalgebra< func::range<F> >,
		                                   lie::coalgebra< func::domain<F> >,
		                                   std::logic_error > {
			default_pull(const F&, const func::domain<F>& )
				: default_pull::base{ std::logic_error( meta::name<F>() + " has no pullback") } {
				
			}      
			
		};
	}
}

#endif
