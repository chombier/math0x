#ifndef GROUP_FUNC_LINE_H
#define GROUP_FUNC_LINE_H

#include <group/euclid.h>
#include <math0x/func/form.h>


namespace func {
  
  template<class E> struct form;

  // vector line: field -> E 
  template<class E>
  struct line {
    typedef line base;
    
    E dir;
    euclid::space<E> space;
    
    line(const E& dir) 
    : dir(dir),
      space(dir) {
      
    }
    
    E operator()(const euclid::field<E>& alpha) const {
      return space.scal(alpha, dir);
    }
    
    struct push;
    struct pull;
  };



  template<class E>
  struct line<E>::push : line<E> {
    push(const line& of, const euclid::field<E>& ) : line(of) { }
  };



  template<class E>
  struct line<E>::pull: form< euclid::dual<E> > {
      
    pull(const line& of, const euclid::field<E>& ) 
      : pull::base(of.dir) { }
      
  };
  
}


#endif
