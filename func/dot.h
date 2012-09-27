#ifndef GROUP_FUNC_DOT_H
#define GROUP_FUNC_DOT_H

#include <group/euclid.h>
#include <group/tuple.h>

#include <group/func/id.h>

#include <group/func/riesz.h>
#include <group/func/form.h>

namespace func {

  // dot product, given metric M: (x, y) -> x^T M y
  template<class E, class Metric = id<E> > 
  struct dot {
    typedef dot base;
    
    riesz<E, Metric> impl;
    
    dot(const euclid::space<E>& space = {},
	const Metric& metric = {}) 
      : impl(space, metric) { }
    
    
    euclid::field<E> operator()(const std::tuple<E, E>& x) const {
      return form<E>( impl( std::get<0>(x)) )( std::get<1>(x));
    };
    
    
    struct push {
      
      riesz<E, Metric> impl;
      form<E> lhs, rhs;
      
      push(const dot& of, const E& at)
	: impl(of.impl),
	  lhs( impl( std::get<1>(at) ) ),
	  rhs( impl( std::get<0>(at) ) ) {
	
      }
      
      euclid::field<E> operator()(const std::tuple<E, E>& dx) const {
	return lhs( std::get<0>(dx) ) + rhs( std::get<1>(dx) );
      }
      
    };
    
    // TODO pull
  
  };

}


#endif
