
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

#include <math0x/func/num.h>

#include <math0x/func/prod.h>
#include <math0x/func/inv.h>

#include <math0x/each.h>

#include <math0x/test/func.h>

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


			struct push {
				domain x;
				push(const delta& , const domain& at) : x(at) {
					
				}

				lie::algebra<range> operator()(const lie::algebra<domain>& v) const {
			
					lie::algebra<range> res; res.resize(v.size());
					
					each(res, [&](NN i) {
							if( !i ) res(i) = v(i);
							else {
								auto fun = func::prod< SE<3> >{} << func::make_tuple( func::inv< SE<3> >{},
								                                                      func::id< SE<3> > {} );
								std::tuple<SE<3>, SE<3> > at(x(i-1), x(i));
								std::tuple<vec6, vec6> vv(v(i-1), v(i));
								res(i) = d(fun)(at)(vv);
							}
						});
					
					return res;
				}
				
			};

			struct pull {
				domain x;
				pull(const delta&, const domain& at) : x(at) {
					
				}

				lie::coalgebra<domain> operator()(const lie::coalgebra<range>& p) const {
					lie::coalgebra<domain> res; res.resize(p.size());
					
					res = euclid::space_of(res).zero();
					
					each(res, [&](NN i) {
							if( !i ) res(i) = p(i);
							else {
								auto fun = func::prod< SE<3> >{} << func::make_tuple( func::inv< SE<3> >{},
								                                                      func::id< SE<3> > {} );
								std::tuple<SE<3>, SE<3> > at(x(i-1), x(i));

								auto pp = dT(fun)(at)(p(i));
								res(i - 1) += std::get<0>(pp);
								res(i) += std::get<1>(pp);
							}
						});
					
					return res;
				}
			};

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


			
			struct push {
				domain x;
				
				push(const fit&, const domain& x) : x(x) { }

				lie::algebra<range> operator()(const lie::algebra<domain>& v) const {

					lie::algebra<range> res;
					res.resize( std::get<2>(v).size() );

					func::apply< SE<3> > apply;
					
					// derp when @each is used
					for(unsigned i = 0, n = res.size(); i < n; ++i) {
						auto at = std::make_tuple(std::get<2>(x)(i), std::get<I>(x));
						res(i) = d(apply)(at)( std::make_tuple( std::get<2>(v)(i), std::get<I>(v)));
					}
					
					return res;
				}
			};


			struct pull {
				domain x;
				pull(const fit&, const domain& x) : x(x) { }
				
				lie::coalgebra<domain> operator()(const lie::coalgebra<range>& p) const {
					lie::coalgebra<domain> res;
					std::get<2>(res).resize( p.size() );
					
					res = euclid::space_of(res).zero();
					
					func::apply< SE<3> > apply;
					
					for( unsigned i = 0, n = p.size(); i < n; ++i) {
						auto at = std::make_tuple(std::get<2>(x)(i), std::get<I>(x));

						auto pp = dT(apply)(at)(p(i));
							
						std::get<2>(res)(i) += std::get<0>(pp);
						std::get<I>(res) += std::get<1>(pp);
					}
					return res;
				}
			};
			
		};



		result_type operator()(const samples_type& ai,
		                       const samples_type& bi) const {
			assert( ai.size() == bi.size() );
			
			NN n = ai.size();

			func::get<result_type, 0> a;
			func::get<result_type, 1> b;
			func::get<result_type, 2> g;
			
			lie::group< vector<SE<3> > > se3_n( n );
			
			result_type res;
			std::get<2>(res) = se3_n.id();
		
			auto full = func::make_tie(a, b, delta{} << g,
			                           fit<0>{}, fit<1>{} );
			typedef decltype(full) full_type;

			if( !test::func(full, lie::group_of(res), lie::group_of(full(res)) ) ) {
				throw std::logic_error("derp !");
			}
			
			auto rhs = func::range< full_type >( vec3::Zero(), vec3::Zero(), 
			                                     se3_n.id(), 
			                                     ai, bi );
			levmar opt;

			opt.outer.bound = 50;
			opt.inner.bound = -1;
			opt.inner.epsilon = 1e-3;
			

			opt.outer.cb = [&](NN i, RR eps) {
				std::cout << "outer: " << i << " eps: " << eps << std::endl;
			};

			opt.inner.cb = [&](NN i, RR eps) {
				std::cout << "inner: " << i << " eps: " << eps << std::endl;
			};

			
			opt.sparse(res, full, rhs);
						
			return res;
		}

	};
}
