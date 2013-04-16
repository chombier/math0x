#ifndef MATH0X_FUNC_PAIR_H
#define MATH0X_FUNC_PAIR_H

#include <math0x/euclid.h>
#include <math0x/tuple.h>
#include <math0x/func/form.h>

namespace math0x {
	namespace func {

		// TODO: might be convenient to abstract all bilinear forms in
		// bilinear<F> ?
		
		// natural pairing: E* x E -> RR
		template<class E>
		struct pair {
			
			typedef std::tuple< euclid::dual<E>, E > domain;
			typedef euclid::field<E> range;

			euclid::space< E> space;
			
		public:

			pair(const euclid::space<E>& space = {}) 
				: space(space) { }
			
			range operator()(const domain& x) const {
				return space.pair( std::get<0>(x), std::get<1>(x) );
			}
			
			// easy-peasy
			struct push : form< domain >{
				
				push(const pair&, const domain& at) 
					: push::base( std::make_tuple( std::get<1>(at), std::get<0>(at))) {
					
				}
				
			};
			
			struct pull : line< euclid::dual<domain> > {
				
				pull(const pair&, const domain& at)
					: pull::base( std::make_tuple( std::get<1>(at), std::get<0>(at)) ) {

				}

			};
			
			
		};

		template<class E>
		pair<E> make_pair(const euclid::space<E>& space) { return {space}; }
		
	}
}


#endif
