#ifndef MATH0X_FUNC_DOT_H
#define MATH0X_FUNC_DOT_H

#include <math0x/euclid.h>
#include <math0x/tuple.h>

#include <math0x/func/id.h>

#include <math0x/func/riesz.h>
#include <math0x/func/form.h>

namespace math0x { 
	namespace func {

		// dot product, given metric M: (x, y) -> x^T M y
		template<class E, class Metric = id<E> > 
		struct dot {
			typedef dot base;
    
			riesz<E, Metric> impl;
    
			dot(const euclid::space<E>& space = {},
			    const Metric& metric = {}) 
				: impl(space, metric) { }
    
			typedef std::tuple<E, E> domain;
			euclid::field<E> operator()(const domain& x) const {
				return form<E>( impl( std::get<0>(x)) )( std::get<1>(x));
			};
    
			struct push {
      
				riesz<E, Metric> impl;
				func::push< form<E> > lhs, rhs;
      
				push(const dot& of, const domain& at)
					: impl(of.impl),
					  lhs( impl( std::get<1>(at) ) ),
					  rhs( impl( std::get<0>(at) ) ) {
	
				}
      
				euclid::field<E> operator()(const domain& dx) const {
					return lhs( std::get<0>(dx) ) + rhs( std::get<1>(dx) );
				}
      
			};
    
			struct pull {

				riesz<E, Metric> impl;
				line< euclid::dual<domain> > res;
				
				pull(const dot& of, const domain& at)
					: impl(of.impl),
					  res( std::make_tuple(impl( std::get<1>(at) ),
																 impl( std::get<0>(at) ) ) ) {
								 
				}
				
				euclid::dual<domain> operator()(const euclid::field<E>& f) const {
					return res(f);
				}
			
			};
  
		};

	}

}
#endif
