#ifndef GROUP_FUNC_NORM2_H
#define GROUP_FUNC_NORM2_H

#include <group/euclid.h>

#include <group/func/id.h>

#include <group/func/riesz.h>
#include <group/func/form.h>

namespace func {

  // squared norm, given a metric M: x -> x^T M x
  template<class E, class Metric = id<E> >
  struct norm2 {
    typedef norm2 self;
    
    riesz<E, Metric> impl;
    
    euclid::field<E> operator()(const E& x) const {
      return form<E>( impl(x) )(x);
    }
    
    struct push : form<E> {
      
      push(const norm2& of, const E& at)
	: push::self( of.impl.dual.scal( 2.0, of.impl(at) ) ) {

      }
      
    };

    struct pull : line< euclid::dual<E> > {
      pull(const norm2& of, const E& at) 
	: pull::self( of.impl.dual.scal( 2.0, of.impl(at) ) ) {
	
      }
    };
    
  };
}

#endif
