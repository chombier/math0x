
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
#include <math0x/func/val.h>
// #include <math0x/func/ops.h>

#include <math0x/each.h>

#include <math0x/test/func.h>

#include <math0x/debug.h>
#include <math0x/macro.h>

namespace math0x {


	namespace func {

		template<>
		struct apply< std::tuple<RR3, SO3> > {

			typedef RR3 vec_type;
			typedef std::tuple<vec_type, SO3> rigid_type;

			typedef std::tuple< rigid_type, vec_type > domain;
			typedef vec_type range;
			
			static get<domain, 0> rigid() { return {};}
			static get<domain, 1> point() { return {}; }
			
			static get<rigid_type, 0> translation() { return {};} 
			static get<rigid_type, 1> rotation() { return {}; }

			range operator()(const domain& x) const {
				const rigid_type& g = std::get<0>(x);

				return translation()(g) + rotation()(g)( std::get<1>(x) );
			}


			struct push {
				
				domain at;
				push(const apply&, const domain& at) : at(at)  { }


				range operator()(const lie::algebra<domain>& v) const {

					auto fun = sum<RR3>{} << make_tie( translation() << rigid(),
					                                   apply<SO3>{} << make_tie( rotation() << rigid(),
					                                                             point() ));
					return d(fun)(at)(v);
					
				}
			};

			struct pull {
				
				domain at;
				pull(const apply&, const domain& at) : at(at)  { }
				
				
				lie::coalgebra<domain> operator()(const lie::coalgebra<range>& f) const {
					
					auto fun = sum<RR3>{} << make_tie( translation() << rigid(),
					                                   apply<SO3>{} << make_tie( rotation() << rigid(),
					                                                             point() ));
					return dT(fun)(at)(f);
				}
				
			};

		};
			


	}



	struct se3_fit {
	
		typedef vector<vec3> samples_type;

		typedef vec3 RR3;

		// typedef SE3 rigid_type;
		typedef std::tuple<RR3, SO3> rigid_type;
		typedef std::tuple< RR, vector< rigid_type > > result_type;
		
		levmar opt;

		struct delta {
			
			typedef vector< rigid_type> domain;
			typedef vector< rigid_type> range;
			

			static auto diff() -> 
				macro_returns( func::prod< rigid_type >{} << func::make_tuple( func::inv< rigid_type >{}, func::id< rigid_type > {} ) );
			
			mutable range res;
			
			const range& operator()(const domain& x) const {
				// range res;
				res.resize( x.size() );
				
				auto fun = diff();
						
				omp_each(res, [&](NN i) {
						if(i) {
							res(i) = fun( std::make_tuple(x(i-1), x(i))); 
						} else {
							res(i) = lie::group<rigid_type>{}.id();
						}
					});
				
				return res;
			}


			struct push {
				domain x;
				
				push(const delta& , const domain& at) : x(at) { }
				push(const delta& , domain&& at) : x( std::move(at) ) { }
				
				mutable lie::algebra<range> res;
				const lie::algebra<range>& operator()(const lie::algebra<domain>& v) const {
			
					// lie::algebra<range> res; 
					res.resize( v.size() );
					
					auto fun = diff();
					
					auto dfun = d(fun);
					
					omp_each(res, [&](NN i) {
							if( !i ) res(i) = lie::group<rigid_type>{}.alg().zero();
							else {
								std::tuple<rigid_type, rigid_type > at(x(i-1), x(i));
								std::tuple< lie::algebra<rigid_type>, lie::algebra<rigid_type> > vv(v(i-1), v(i));
								res(i) = dfun(at)(vv);
							}
						});
					
					return res;
				}
				
			};

			struct pull {
				domain x;

				pull(const delta&, const domain& at) : x(at) { }
				pull(const delta&, domain&& at) : x( std::move(at) ) { }

				mutable lie::coalgebra<domain> res;
				const lie::coalgebra<domain>& operator()(const lie::coalgebra<range>& p) const {
					// lie::coalgebra<domain> res;
					res.resize(p.size());
					
					res = euclid::space_of(res).zero();

					auto fun = diff();
							
					auto dTfun = dT( fun );
					auto coalg = *lie::group<rigid_type>{}.alg();
					
					
					// each(res, [&](NN i) {
					// 		if(!i) return;
							
					// 		std::tuple<rigid_type, rigid_type > at(x(i-1), x(i));
							
					// 		auto pp = dTfun(at)(p(i));
					// 		res(i - 1) = coalg.sum(res(i - 1), std::get<0>(pp));
					// 		res(i) = coalg.sum( res(i), std::get<1>(pp));

					// 	});

					omp_each(res, [&](NN i) {
							if(i) {
								std::tuple<rigid_type, rigid_type > at(x(i-1), x(i));
								auto pp = dTfun(at)(p(i));
								res(i) = coalg.sum( res(i), std::get<1>(pp));
							}

							if( i < (res.size() - 1) ) {
								std::tuple<rigid_type, rigid_type > at(x(i), x(i + 1));
								auto pp = dTfun(at)(p(i + 1));
								
								res(i) = coalg.sum(res(i), std::get<0>(pp));
							}
						});


					
					return res;
				}
			};

		};
 

		template<bool flip>
		struct fit {
			
			typedef result_type domain;
			typedef samples_type range;

			static RR sign() { return flip ? -1 : 1; }
			
			mutable range res;
			const range& operator()(const domain& x) const {
				// range res;
				res.resize( std::get<1>(x).size() );

				auto fun = func::apply< rigid_type >{} << func::make_tuple( func::id< rigid_type >{},
				                                                            func::line<vec3>{ sign() * vec3::UnitY()} );
				omp_each(res, [&](NN i) {
						res(i) = fun( std::make_tuple(std::get<1>(x)(i),
						                              std::get<0>(x)) );
					});
				return res;
			}

			struct push {
				domain x;
				
				push(const fit&, const domain& x) : x(x) { }
				push(const fit&, domain&& x) : x( std::move(x) ) { }

				mutable lie::algebra<range> res;

				const lie::algebra<range>& operator()(const lie::algebra<domain>& v) const {

					// lie::algebra<range> res;
					res.resize( std::get<1>(v).size() );
					
					auto fun = func::apply< rigid_type >{} << func::make_tuple( func::id< rigid_type >{},
					                                                            func::line<vec3>{ sign() * vec3::UnitY()} );
					
					// compilation derp when @each is used
					unsigned n = res.size();

#pragma omp parallel for
					for(unsigned i = 0; i < n; ++i) {
						auto at = std::make_tuple(std::get<1>(x)(i), std::get<0>(x));
						res(i) = d(fun)(at)( std::make_tuple( std::get<1>(v)(i), std::get<0>(v)));
					}
					
					return res;
				}
			};
			
			
			struct pull {
				domain x;
				
				pull(const fit&, const domain& x) : x(x) { }
				pull(const fit&, domain&& x) : x( std::move(x) ) { }
				
				lie::coalgebra<domain> operator()(const lie::coalgebra<range>& p) const {
					lie::coalgebra<domain> res;
					std::get<1>(res).resize( p.size() );
					
					res = euclid::space_of(res).zero();
					
					auto fun = func::apply<rigid_type>{} << func::make_tuple( func::id< rigid_type >{},
					                                                          func::line<vec3>{ sign() * vec3::UnitY() } );

					lie::group< rigid_type > g;
					auto g_star = *g.alg();
					
					for( unsigned i = 0, n = p.size(); i < n; ++i) {
						auto at = std::make_tuple(std::get<1>(x)(i), std::get<0>(x));
						
						auto pp = dT(fun)(at)(p(i));
						
						std::get<1>(res)(i) = g_star.sum( std::get<1>(res)(i), std::get<0>(pp));
						std::get<0>(res) += std::get<1>(pp);
					}
				
					return res;
				}
			};
			
		};


		// main call
		result_type operator()(const samples_type& ai,
		                       const samples_type& bi) const {
			assert( ai.size() == bi.size() );
			
			NN n = ai.size();

			func::get<result_type, 0> a;
			func::get<result_type, 1> g;
			
			lie::group< vector< rigid_type > > se3_n( n );
			
			result_type res;
			std::get<1>(res) = se3_n.id();
			
			// first guess for translation
			for(NN i = 0; i < n; ++i) {
				// std::get<2>(res)(i).translation = 0.5 * (ai(i) + bi(i));
				std::get<0>(std::get<1>(res)(i)) = 0.5 * (ai(i) + bi(i));
			}

			std::get<0>(res) = 0.5 * (ai(0) - bi(0)).norm();
			
			struct {
				real fit;
				real smooth;
			} weight;

			weight.fit = 1;

			euclid::space< vector<vec3> > RR3_n(n);
			
			auto full = func::make_tie( // func::scal<RR>(small) << a,
			                            // func::scal<RR>(small)<< b,  
			                            // func::get< vector<SE<3> > >{0} << g, 
			                            delta{} << g,
			                            fit<0>{}, 
			                            fit<1>{} );
			
			typedef decltype(full) full_type;

			// test::func(full, lie::group_of(res), lie::group_of(full(res)));
			
			auto rhs = func::range< full_type >( // 0.0, // 0.0,
			                                     se3_n.id(), 
			                                     ai, 
			                                     bi );
			
			
			debug("search space dim:", lie::group_of(res).alg().dim() );

			opt.sparse(res, full, rhs);
			// opt.dense(res, full, rhs);
						
			debug("length:", 2 * std::get<0>(res),  "should be approx.", (ai(0) - bi(0)).norm());
			
			return res;
		}

	};
} 
