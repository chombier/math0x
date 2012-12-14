
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
#include <math0x/func/scal.h>

#include <math0x/each.h>

#include <math0x/test/func.h>

#include <math0x/debug.h>

namespace math0x {
	struct se3_fit {
	
		typedef vector<vec3> samples_type;
		typedef std::tuple< RR, RR, vector<SE<3>> > result_type;
		
		levmar opt;

		struct delta {
			
			typedef vector<SE<3>> domain;
			typedef vector<SE<3>> range;
			
			range operator()(const domain& x) const {
				range res; res.resize( x.size() );
				
				each(res, [&](NN i) {
						if(i) {
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
							if( !i ) res(i).setZero();
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
							if(i) {
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

				auto fun = func::apply<SE<3> >{} << func::make_tuple( func::id< SE<3> >{},
				                                                      func::line<vec3>{vec3::UnitY()} );
				each(res, [&](NN i) {
						res(i) = fun( std::make_tuple(std::get<2>(x)(i),
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

					auto fun = func::apply<SE<3> >{} << func::make_tuple( func::id< SE<3> >{},
					                                                      func::line<vec3>{vec3::UnitY()} );
		
					// derp when @each is used
					for(unsigned i = 0, n = res.size(); i < n; ++i) {
						auto at = std::make_tuple(std::get<2>(x)(i), std::get<I>(x));
						res(i) = d(fun)(at)( std::make_tuple( std::get<2>(v)(i), std::get<I>(v)));
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
					
					auto fun = func::apply<SE<3> >{} << func::make_tuple( func::id< SE<3> >{},
					                                                      func::line<vec3>{ vec3::UnitY() } );
		
					for( unsigned i = 0, n = p.size(); i < n; ++i) {
						auto at = std::make_tuple(std::get<2>(x)(i), std::get<I>(x));

						auto pp = dT(fun)(at)(p(i));
							
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
		
			// first guess for translation
			for(NN i = 0; i < n; ++i) {
				std::get<2>(res)(i).translation = 0.5 * (ai(i) + bi(i));
			}


			real small = 1;
			auto full = func::make_tie( func::scal<RR>(small)<< a,  
			                            func::scal<RR>(small)<< b,  
			                            // func::get< vector<SE<3> > >{0} << g, 
			                            // delta{} << g,
			                            fit<0>{}, fit<1>{} );
			typedef decltype(full) full_type;

			// test::func(full, lie::group_of(res), lie::group_of(full(res)));
			
			auto rhs = func::range< full_type >( 0.0, 0.0,
			                                     // se3_n.id(), 
			                                     ai, bi );
		
	
			debug("search space dim:", lie::group_of(res).alg().dim() );

			opt.sparse(res, full, rhs);
						
			return res;
		}

	};
}
