#ifndef GROUP_FUNC_SCAL_H
#define GROUP_FUNC_SCAL_H

#include <group/euclid.h>

namespace func {

  template<class E>
  struct scal {
    typedef scal self;
    
    euclid::field<E> lambda;
    euclid::space<E> space;
    
    scal( const euclid::field<E>& lambda,
	  const euclid::space<E>& space = euclid::space<E>() ) 
      : lambda( lambda ),
	space( space ) 
    {
      
    }

    
    E operator()(const E& x) const { return space.scal(lambda, x); }
    
    struct push : scal { 
      push(const scal& of, const E& ) : push::self(of) { }
    };

    struct pull : scal< euclid::dual<E> > { 
      pull(const scal& of, const E& ) : pull::self(of.lambda, *of.space) { }
    };
    
  };


}


#endif
