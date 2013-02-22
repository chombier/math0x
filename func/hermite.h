#ifndef MATH0X_FUNC_HERMITE_H
#define MATH0X_FUNC_HERMITE_H

#include <math0x/func/id.h>
#include <math0x/func/val.h>
#include <math0x/func/poly.h>
#include <math0x/func/operators.h>
#include <math0x/macro.h>

namespace math0x {
	namespace func {
		
		template<class U = RR>
		class hermite {
			static constexpr id<U> t = {};
			typedef value<U, U> val;

		public:
			
			static auto h00() -> macro_returns( 2.0 * (t^3) - 3.0 * (t^2) + val(1.0) );

			static auto h10() -> macro_returns( (t^3) - 2.0 * (t^2) + t );
			
			static auto h01() -> macro_returns( -2.0 * (t^3) + 3.0 * (t^2) );

			static auto h11() -> macro_returns( (t^3) - (t^2) );
		
		};

	}
}


#endif
