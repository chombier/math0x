#ifndef GROUP_FUNC_PULL_H
#define GROUP_FUNC_PULL_H

#include <group/func.h>
#include <group/meta.h>

namespace func {

  template<class F>
  struct pullback {
    F of;
      
    pull<F> operator()(const domain<F>& x) const {
      return {of, x};
    }
      
  };
    
  template<class F>
  pullback< meta::decay<F> > dT(F&& f) { return {std::forward<F>(f)}; }
  
}


#endif
