#ifndef MATH0X_TEST_PULL_H
#define MATH0X_TEST_PULL_H

#include <math0x/func.h>
#include <math0x/func/pull.h>

#include <math0x/test/random.h>

#include <math0x/func/num.h>
#include <math0x/func/num.h>
#include <math0x/func/norm2.h>
#include <math0x/func/minus.h>

namespace math0x {

	namespace test {

		// numerically compare a function pullback with its numerical
		// derivative, using random base-point/tangent vector
		template<class F>
		RR pull(const F& f, RR step,
		        const lie::group< func::domain<F> >& dmn = {},
		        const lie::group< func::range<F> >& rng = {} ) {
			
			// random base-point
			func::domain<F> at = random(dmn);

			// random tangent vector
			lie::coalgebra< func::range<F> > v = random_cotangent( rng );
			
			using namespace func;			
			auto num = dT( make_num(f, step) )(at);
			auto dTf = dT(f)(at);
			
			typedef decltype( num ) num_type;
			typedef decltype( dTf ) dTf_type;

			auto diff = make_sum(*dmn.alg()) << make_tie( dTf, (make_minus(*dmn.alg()) << num) );
			
			return std::sqrt( (make_norm2(*dmn.alg()) << diff)(v) );
		}
		
	}
}


#endif
