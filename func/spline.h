#ifndef MATH0X_FUNC_SPLINE_H
#define MATH0X_FUNC_SPLINE_H

#include <math0x/func/hermite.h>

namespace math0x {
	namespace func {
	
		template<class U = RR>
		class spline {
			static constexpr id<U> t = {};
			typedef value<U, U> val;

			static U width(U start, U end) { 
				assert( start < end && "bounds must be increasing lol");
				return end - start;
			}
				
			static auto scaling(U start, U end) -> macro_returns( (t - val(start)) / width(start, end) );
				
			public:
			
			typedef func::hermite<U> basis;
		public:

			// hermite spline coefficients, see
			// http://en.wikipedia.org/wiki/Cubic_Hermite_spline
			// p(t) = s0(t) * p0 + s1(t) * m0 + s2(t) * p1 + s3(t) * m1
			// coefficients for points(pi) and tangents(mi)
			static auto hermite(U start = 0, U end = 1) -> 
				macro_returns( make_tie( basis::h00(),
				                         width(start, end) * basis::h10(),
				                         basis::h01(),
				                         width(start, end) * basis::h11()) << scaling(start, end) );
		};
			
	}
}


#endif
