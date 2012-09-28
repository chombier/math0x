#ifndef MATH0X_FUNC_PULL_H
#define MATH0X_FUNC_PULL_H

#include <math0x/func.h>
#include <math0x/meta.h>

namespace math0x { 
  namespace func {

    template<class F>
    struct pullback {
      F of;
      
      pull<F> operator()(const domain<F>& x) const {
	return pull<F>(of, x);
      }
      
    };
    
    template<class F>
    pullback< meta::decay<F> > dT(F&& f) { return {std::forward<F>(f)}; }
  
  }

}
#endif
