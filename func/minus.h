#ifndef GROUP_FUNC_MINUS_H
#define GROUP_FUNC_MINUS_H

#include <group/euclid.h>

namespace func {

  template<class E>
  struct minus {
    typedef minus base;
    
    euclid::space<E> space;

    E operator()(const E& x) const {
      return space.minus(x);
    }

    
    struct push;
    struct pull;
    
  };
  

  template<class E>
  struct minus<E>::push : minus {
      
    push(const minus& of, const E& ) : push::base(of) { }
    
  };

  template<class E>
  struct minus<E>::pull : minus< euclid::dual<E> > {
    
    pull(const minus& of, const E& ) : pull::base{ *of.space } {  }
    
  };
  
}


#endif
