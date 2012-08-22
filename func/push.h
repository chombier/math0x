#ifndef GROUP_FUNC_PUSH_H
#define GROUP_FUNC_PUSH_H

#include <group/func.h>
#include <group/meta.h>

namespace func {
  
  // differential of a function: d(f)(x) returns the derivative of f
  // in the direction x
  template<class F>
  struct pushforward {
    F of;
      
    push<F> operator()(const domain<F>& x) const {
      return {of, x};
    }
    
  };
  
  template<class F>
  pushforward< meta::decay<F> > d(F&& f) { return {std::forward<F>(f)}; }
}


#endif
