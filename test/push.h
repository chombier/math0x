#ifndef MATH0X_TEST_PUSH_H
#define MATH0X_TEST_PUSH_H

#include <math0x/func.h>

#include <math0x/test/random.h>

#include <math0x/func/num.h>
#include <math0x/func/num.h>
#include <math0x/func/norm2.h>
#include <math0x/func/minus.h>
#include <math0x/func/sum.h>
#include <math0x/func/push.h>

#include <math0x/error.h>
#include <math0x/func/operators.h>


#include <iostream>

namespace math0x {

	namespace test {

		// numerically compare a function pushforward with its numerical
		// derivative, using random base-point/tangent vector
		template<class F>
		RR push(const F& f, RR step,
		        const lie::group< func::domain<F> >& dmn = {},
		        const lie::group< func::range<F> >& rng = {} ) {
			
			// random base-point
			func::domain<F> at = random(dmn);

			// random tangent vector
			lie::algebra< func::domain<F> > v = random_tangent(dmn);

			using namespace func;			
			auto num = d( make_num(f, step) )(at);
			auto df = d(f)(at);

			auto diff = make_sum(rng.alg()) << 	make_tie( df, make_minus(rng.alg()) <<  num );

			auto norm2 = make_norm2( rng.alg() );
			


			try {
				RR den = std::sqrt( (norm2 << num)(v) );
				
				// return relative error
				RR num = std::sqrt( (norm2 << diff)(v) );

				if( den < math0x::epsilon() ) {
					std::cerr << "warning: div by zero lol" << std::endl;
					return num;
				}
				return num / den;				
			} 
			catch( const math0x::error& e ){
				std::cerr << "exception caught: " << e.what() << std::endl;
				return 42;
			}
			
		}
		
	}
}


#endif
