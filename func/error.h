#ifndef GROUP_FUNC_ERROR_H
#define GROUP_FUNC_ERROR_H

#include <group/types.h>
#include <group/lie.h>

namespace func {

  // throw an error
  template<class Domain, class Range, class What >
  struct error {
    typedef error self;
    
    What what;
    
    Range operator()(const Domain& ) const { 
      throw what;
    }
    
    
    struct push : error< lie::algebra<Domain>, lie::algebra<Range>, What > {
      
      push(const error& of, const Domain& ) : push::self{of.what} { }
      
    };


    struct pull : error< lie::coalgebra<Range>, lie::coalgebra<Domain>, What > {
      pull(const error& of, const Domain& ) : pull::self{of.what} { }
    };

  };


};


#endif
