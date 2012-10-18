#ifndef MATH0X_TEST_FUNC_H
#define MATH0X_TEST_FUNC_H

#include <math0x/test/push.h>
#include <math0x/test/pull.h>

#include <math0x/debug.h>
#include <math0x/meta.h>

namespace math0x {

	namespace test {

		
		template<class F>
		void func(const F& f, 
		          const lie::group< func::domain<F> >& dmn = {},
		          const lie::group< func::range<F> >& rng = {},
		          RR epsilon = 1e-8, RR step = 1e-5) {
			
			debug("testing",meta::name<F>());
			
			// pushforward
			{
				RR error = push(f, step, dmn, rng);
				
				if( error > epsilon ) {
					debug("push error:", error);
				} else {
					debug("push ok");
				}
			}

			// pullback
			{
				RR error = pull(f, step, dmn, rng);
				
				if( error > epsilon ) {
					debug("pull error:", error);
				} else {
					debug("pull ok");
				}
			}
			
		}


	}

}


#endif
