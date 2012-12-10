#ifndef MATH0X_FUNC_COORD_H
#define MATH0X_FUNC_COORD_H

#include <math0x/vector.h>
#include <math0x/euclid.h>
#include <math0x/func.h>

#include <math0x/meta.h>
#include <math0x/coords.h>

namespace math0x {

	namespace func {

		// express an euclidean function in coordinates
		template<class F>
		struct coord {

			F of;
			euclid::space< func::domain<F> > dmn;
			euclid::space< func::range<F> > rng;
		
			typedef euclid::coords< func::domain<F> > domain;
			typedef euclid::coords< func::range<F> > range;
		
			range operator()(const domain& x) const {
				func::domain<F> xx = dmn.zero();
			
				dmn.set(xx, x);
			
				range res; res.resize( rng.dim() );
				rng.get(res, of(xx) );
			
				return res; 

			}


			// TODO push/pull ?
		
		};

		template<class F>
		coord< meta::decay<F> > make_coords(F&& of, 
		                                    const euclid::space< func::domain< meta::decay<F> > > & dmn = {},
		                                    const euclid::space< func::range< meta::decay<F> > > & rng = {} ) {
			return { std::forward<F>(of), dmn, rng};
		}
	
	}

}
#endif
