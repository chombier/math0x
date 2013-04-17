#ifndef MATH0X_TEST_FUNC_H
#define MATH0X_TEST_FUNC_H

#include <math0x/test/push.h>
#include <math0x/test/pull.h>

#include <math0x/debug.h>
#include <math0x/meta.h>

#include <math0x/test/io.h>

namespace math0x {

	namespace test {

		template<class F>
		bool func(const F& f = {}, 
		          const lie::group< func::domain<F> >& dmn = {},
		          const lie::group< func::range<F> >& rng = {},
		          RR epsilon = 1e-7, RR step = 1e-4) {
			
			bool res = true;
			
			debug("testing", meta::name<F>());
			
			// pushforward
			{
				RR error = push(f, step, dmn, rng);
				
				if( error > epsilon ) {
					debug(fail(), "push error:", error);
					res = false;
				} else {
					debug(ok(), "push ok");
				}
			}

			// pullback
			{
				RR error = pull(f, step, dmn, rng);
				
				if( error > epsilon ) {
					debug(fail(), "pull error:", error);
					res = false;
				} else {
					debug(ok(), "pull ok");
				}
			}
			
			return res;
		}


	}

}


#endif
