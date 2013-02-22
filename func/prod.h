#ifndef MATH0X_FUNC_PROD_H
#define MATH0X_FUNC_PROD_H

#include <math0x/lie.h>
#include <math0x/tuple/repeat.h>

namespace math0x { 
	namespace func {


		// lie group product
		template<class G, int N = 2>
		struct prod {
			typedef prod base;
    
			lie::group<G> group;
    
			prod(const lie::group<G>& group = {} ) : group(group) { }
			
			static_assert( N >= 2, "N must be greater than 2");

			typedef math0x::tuple::repeat<G, N> domain;
			typedef math0x::tuple::range< N > each;
			
			struct call {
				G& res;
				const lie::group<G>& group;
				const domain& x;
	  
				template< int I >
				void operator()() const {
					if( I > 0 ) res = group.prod(res, std::get<I>(x));
				}
			};

			G operator()(const domain& g) const {
				G res = std::get<0>(g);
				each::apply( call{res, group, g} );
				return res;
			}
    

			struct fill_base {
				domain& at;							// id, x^-1(0), x^-1(1).x^-1(0), ...
				
				const lie::group<G>& group;
				const domain& x;
				
				template< int I >
				void operator()() const {
					if( I == 0 ) std::get<N - 1>(at) = group.id();
					else if (I == 1) std::get<N - 2>(at) = group.inv(std::get<N-1>(x));
					else {
						static constexpr int II = !I ? N - 1 : N - I;
						
						std::get<N - 1 - I>(at) = group.prod( std::get<II>(at),
						                                      group.inv( std::get<II>(x) ) );
					}
				}
				
			};
			
			// TODO skip this when group is trivial 
			static domain make_base(const domain& x,
			                        const lie::group<G>& group) {
				domain res;
				each::apply( fill_base{res, group, x} );
				
				return res;
			}
			
    
			struct push {

				lie::Ad< domain > Ad;
				prod< lie::algebra<G>, N > sum;

				push(const prod& of, const domain& at)
					: Ad( make_base(at, of.group ) ),
					  sum( of.group.alg().zero() ) // TODO optimize this ?!
				{
					
				}

				lie::algebra<G> operator()(const lie::algebra<domain>& v) const {
					return sum( Ad(v) );
				}
				
			};



			struct pull : meta::unpack_args<func::tie, typename lie::AdT<domain>::args_type>  {

					
				pull(const prod& of, const domain& at) 
					: pull::base( lie::AdT<domain>( make_base(at, of.group) ).args ) {
					
				}
				
			};
    

		};


		// convenience 
		template<class G>
		prod<G> make_prod(const lie::group<G>& group) { return {group}; }
		

	}

}
#endif
