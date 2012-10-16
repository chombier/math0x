#ifndef MATH0X_TEST_FUNC_H
#define MATH0X_TEST_FUNC_H

#include <math0x/debug.h>
#include <math0x/test/func.h>
#include <math0x/meta.h>

namespace math0x {

	namespace test {

		
		template<class F>
		void func(const F& f, 
		          const lie::group< func::domain<F> >& dmn = {},
		          const lie::group< func::range<F> >& rng = {},
		          RR epsilon = 1e-7, RR step = 1e-5) {
			
			debug("testing",meta::name<F>());
			
			{
				RR error = push(f, step, dmn, rng);
				
				if( error > epsilon ) {
					debug("push error:", error);
				} else {
					debug("push ok");
				}
			}
			
		}


	}

}


#endif
