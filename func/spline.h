#ifndef MATH0X_FUNC_SPLINE_H
#define MATH0X_FUNC_SPLINE_H

#include <math0x/func/hermite.h>
#include <math0x/func/get.h>
#include <math0x/func/subgroup.h>
#include <math0x/func/prod.h>

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
				
			static auto scaling(U start, U end) -> 
				macro_returns( (t - val(start)) / width(start, end) );
			
			typedef math0x::tuple::repeat<U, 4> coeffs;
			
			static constexpr get<coeffs, 0> get0 = {};
			static constexpr get<coeffs, 1> get1 = {};
			static constexpr get<coeffs, 2> get2 = {};
			static constexpr get<coeffs, 3> get3 = {};
			

			static auto hermite_to_cardinal(U alpha) -> 
			macro_returns( make_tie( -alpha * get1,
			                         get0 - alpha * get3,
			                         get2 + alpha * get1,
			                         alpha * get3 ) );
			
			typedef func::hermite<U> basis;

			// cumulative form conversion, see http://portal.acm.org/citation.cfm?id=218380.218486
			static auto cumulative() ->
				macro_returns( make_tie(get3,
				                        get3 + get2,
				                        get2 + get1,
				                        get1 + get0) );


			template<class G>
			static auto subgroups(const G& g0, 
			                      const G& g1, 
			                      const G& g2, 
			                      const G& g3,
			                      const lie::group<G>& group,
			                      const lie::log<G>& log) ->
				macro_returns( make_tuple( subgroup<G>( log( g0 ), group ),
				                           subgroup<G>( log( group.prod( group.inv(g0), g1) ), group ),
				                           subgroup<G>( log( group.prod( group.inv(g1), g2) ), group ),
				                           subgroup<G>( log( group.prod( group.inv(g2), g3) ), group ) )				              
				               );


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

			// cardinal spline with tension, see TODO (UNEQUALLY spaced nodes !)
			static auto cardinal(U start = 0, U end = 1, U tension = 0) -> 
				macro_returns( hermite_to_cardinal( (1.0 - tension) / (2 * width(start, end)) )
				               << hermite(start, end) );
			
			// feed this with coefficients above to get a general lie group
			// spline patch based on [Kim95] quaternion splines
			template<class G>
			static auto patch(const G& g0,
			                  const G& g1,
			                  const G& g2,
			                  const G& g3,
			                  const lie::group<G>& group = {} ) ->
				macro_returns( prod<G, 4>(group) << subgroups(g0, g1, g2, g3, group, group.log()) << cumulative() );
		};
		
		

	}
}


#endif
