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
    
    
    struct push : error< lie::alg<Domain>, lie::alg<Range>, What > {
      
      push(const error& of, const Domain& ) : push::self{of.what} { }
      
    };


    struct pull : error< lie::coalg<Range>, lie::coalg<Domain>, What > {
      pull(const error& of, const Domain& ) : pull::self{of.what} { }
    };

  };


};


#endif
