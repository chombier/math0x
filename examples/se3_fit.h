
#include <math0x/types.h>

#include <math0x/vector.h>
#include <math0x/SE3.h>

#include <math0x/func/get.h>
#include <math0x/func/tie.h>
#include <math0x/func/comp.h>
#include <math0x/func/default.h>

#include <math0x/real.h>
#include <math0x/levmar.h>
#include <math0x/func/tuple.h>

#include <math0x/func/push.h>
#include <math0x/func/pull.h>

namespace math0x {
	struct se3_fit {
	
		typedef vector<vec3> samples_type;
		typedef std::tuple< vec3, vec3, vector<SE<3>> > result_type;
		

		struct delta {
			
			typedef vector<SE<3>> domain;
			typedef vector<SE<3>> range;
			
			range operator()(const domain& x) const {
				range res; res.resize( x.size() );
				
				each(res, [&](NN i) {
						if( !i ) res(i) = x(i);
						else {
							res(i) = x(i-1).inv() * x(i); 
						}
					});
				
				return res;
			}


		};


		template<int I>
		struct fit {

			typedef result_type domain;
			typedef samples_type range;

			range operator()(const domain& x) const {
				range res;
				res.resize( std::get<2>(x).size() );
				func::apply< SE<3> > apply;
				
				each(res, [&](NN i) {
						res(i) = apply( std::make_tuple(std::get<2>(x)(i),
						                                std::get<I>(x)) );
					});
				return res;
			}
			
		};



		result_type operator()(const samples_type& ai,
		                       const samples_type& bi) const {
			assert( ai.size() == bi.size() );
			
			NN n = ai.size();

			func::get<result_type, 0> a;
			func::get<result_type, 1> b;
			func::get<result_type, 2> g;
			
			// gi.a = ai
			
			lie::group< vector<SE<3> > > se3_n( n );
			
			auto full = func::make_tie(a, b, delta{} << g,
			                           fit<0>{}, fit<1>{} );
			typedef decltype(full) full_type;

			auto rhs = func::range< full_type >( vec3::Zero(), vec3::Zero(), 
			                                     se3_n.id(), 
			                                     ai, bi);

			
			levmar opt;

			opt.outer.bound = 10;
			opt.inner.bound = 10;
			
			result_type res;
			std::get<2>(res) = se3_n.id();
			
			opt.sparse(res, full, rhs);
						
			return res;
		}

	};
}
