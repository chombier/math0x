#ifndef GROUP_FUNC_LINE_H
#define GROUP_FUNC_LINE_H

#include <group/euclid.h>

// vector line
namespace func {
  
  template<class E>
  struct line {
    typedef line self;
    
    E dir;
    euclid::space<E> space;
    
    line(const E& dir) 
    : dir(dir),
      space(dir) {
      
    }
    
    E operator()(const euclid::field<E>& alpha) const {
      return space.scal(alpha, dir);
    }
    
    struct push : line {
      
      push(const line& of, const euclid::field<E>& ) : line(of) { }
      
    };
    
    struct pull : form< euclid::dual<E> > {
      
      pull(const line& of, const euclid::field<E>& ) 
	: pull::self(of.dir) { }
      
    };
    
  };

}


#endif
