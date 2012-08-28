#ifndef GROUP_FUNC_FORM_H
#define GROUP_FUNC_FORM_H

#include <group/euclid.h>

// linear form over a vector space E 
namespace func {
  
  template<class E>
  struct form {
    typedef form self;
    
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
    
    
    struct push : form {
      
      push(const form& of, const E& ) : form(of) { }
      
    };


    struct pull : line< euclid::dual<E> > {
      
      pull( const form& of, const E& )
	: pull::self(of.value) { }
      
    };
    
    
  };

}


#endif
