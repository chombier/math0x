#ifndef MATH0X_FUNC_NORM2_H
#define MATH0X_FUNC_NORM2_H

#include <math0x/euclid.h>

#include <math0x/func/id.h>

#include <math0x/func/riesz.h>
#include <math0x/func/form.h>


namespace math0x { 
	namespace func {

		// squared norm, given a metric M: x -> x^T M x
		template<class E, class Metric = id<E> >
		struct norm2 {
			typedef norm2 base;
    
			riesz<E, Metric> impl;
    
			norm2(const euclid::space<E>& space = {},
			      const Metric& metric = {} ) : impl(space, metric) { }
    
			euclid::field<E> operator()(const E& x) const {
				return form<E>( impl(x) )(x);
			}
    
			struct push : form<E> {
      
				push(const norm2& of, const E& at)
					: push::base( of.impl.dual.scal( 2.0, of.impl(at) ) ) {

				}
      
			};

			struct pull : line< euclid::dual<E> > {
				pull(const norm2& of, const E& at) 
					: pull::base( of.impl.dual.scal( 2.0, of.impl(at) ) ) {
	
				}
			};
    
		};

		// convenience
		template<class E>
		norm2<E> make_norm2(const euclid::space<E>& space) { return {space}; }

	}
}
#endif
