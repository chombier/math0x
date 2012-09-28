#ifndef MATH0X_FUNC_FORM_H
#define MATH0X_FUNC_FORM_H

#include <math0x/euclid.h>
#include <math0x/func/line.h>

namespace func {

  template<class E> struct line;

  // linear form over a vector space: E -> field
  template<class E>
  struct form {
    typedef form base;
    
    euclid::dual<E> value;
    
    euclid::space< euclid::dual<E> > dual;
    euclid::space<E> primal;
    
    form(const euclid::dual<E>& value)
      : value( value ),
	dual( value ),
	primal( *dual ) {
      
    }

    form(euclid::dual<E>&& value)
      : value( std::move(value) ),
	dual( value ),
	primal( *dual ) {
      
    }
    
    euclid::field<E> operator()(const E& x) const {
      euclid::field<E> res = 0;
      
      NN i = 0;
      dual.each(value, [&](const euclid::field<E>& vi) {
	  res += vi * primal.coord(i, x);
	  ++i;
	});
      
      return res;      
    }
    
    
    struct push;
    struct pull;
    
    
  };

  template<class E>
  struct form<E>::push : form {
    push(const form& of, const E& ) : form(of) { }
  };


  template<class E>
  struct form<E>::pull : line< euclid::dual<E> > {
      
    pull( const form& of, const E& )
      : pull::base(of.value) { }
      
  };

}


#endif
