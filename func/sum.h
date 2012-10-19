#ifndef MATH0X_FUNC_SUM_H
#define MATH0X_FUNC_SUM_H

#include <math0x/euclid.h>

#include <math0x/func/tie.h>
#include <math0x/func/id.h>

namespace math0x { 
	namespace func {

		template<class E>
		struct sum {
			typedef sum base;
    
			euclid::space<E> space;
			sum( const euclid::space<E>& space = {}) : space(space) { }
			
			typedef std::tuple<E, E> domain;

			E operator()(const domain& x) const {
				return space.sum( std::get<0>(x), std::get<1>(x) );
			}
    
			struct push;
			struct pull;
    
		};



		template<class E>
		struct sum<E>::push : sum {
			push(const sum& of, const domain& ) : push::base(of) { }
		};



		template<class E>
		struct sum<E>::pull : tie< id< euclid::dual<E> >,
		                           id< euclid::dual<E> > > {

			pull(const sum&, const domain& )  { }

		}; 

		
		template<class E>
		sum<E> make_sum(const euclid::space<E>& space) { return {space}; }


	}
}
#endif
