#ifndef MATH0X_FUNC_SPLINE_H
#define MATH0X_FUNC_SPLINE_H

#include <math0x/func/hermite.h>
#include <math0x/func/get.h>
#include <math0x/func/subgroup.h>
#include <math0x/func/prod.h>

namespace math0x {
	namespace func {
	
		// some spline-related tools
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
			
			// width is between t-1 and t1 !
			static U alpha( U tension, U width ) {
				return (1 - tension) / width;
			}
			

			// rearrange hermite coefficients to obtain cardinal ones,
			// uniform. alpha is (1 - tension) / (2 * width)
			static auto hermite_to_cardinal(U alpha) -> 
			macro_returns( make_tie( -alpha * get1,
			                         get0 - alpha * get3,
			                         get2 + alpha * get1,
			                         alpha * get3 ) );

			//  rearrange hermite coefficients to obtain cardinal ones,
			// non_uniform. alpha_prev = (1 - tension) / (t1 - t-1)
			// alpha_next = (1 - tension) / (t2 - t0)
			static auto hermite_to_cardinal_nonuniform(U alpha_prev, U alpha_next) -> 
				macro_returns( make_tie( -alpha_prev * get1,
				                         get0 - alpha_next * get3,
				                         get2 + alpha_prev * get1,
				                         alpha_next * get3 ) );

	
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

					
			typedef func::hermite<U> basis;

		public:

			// hermite spline coefficients for interpolating between
			// start/end, see
			// http://en.wikipedia.org/wiki/Cubic_Hermite_spline p(t) =
			// s0(t) * p0 + s1(t) * m0 + s2(t) * p1 + s3(t) * m1
			// coefficients for points(pi) and tangents(mi)
			static auto hermite(U start = 0, U end = 1) -> 
				macro_returns( make_tie( basis::h00(),
				                         width(start, end) * basis::h10(),
				                         basis::h01(),
				                         width(start, end) * basis::h11() ) << scaling(start, end) );

			// cardinal spline with tension, uniform nodes
			static auto cardinal(U start = 0, U end = 1, U tension = 0) -> 
				macro_returns( hermite_to_cardinal( alpha(tension, 2 * width(start, end) ) )
				               << hermite(start, end) );

			static auto cardinal(U t_m1, U t_0, U t_p1, U t_p2, U tension = 0) -> 
				macro_returns( hermite_to_cardinal( alpha(tension, width(t_m1, t_p1) ),
				                                    alpha(tension, width(t_p2, t_0) ) )
				               << hermite(start, end) );

			
			// feed this with coefficients functions above to get a general
			// lie group spline patch based on [Kim95] quaternion splines
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
