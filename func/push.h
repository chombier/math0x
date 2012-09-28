#ifndef MATH0X_FUNC_PUSH_H
#define MATH0X_FUNC_PUSH_H

#include <math0x/func.h>
#include <math0x/meta.h>

namespace math0x { 
  namespace func {
  
    // differential of a function: d(f)(x) returns the derivative of f
    // in the direction x
    template<class F>
    struct pushforward {
      F of;
      
      push<F> operator()(const domain<F>& x) const {
	return push<F>(of, x);
      }
    
    };
  
    template<class F>
    pushforward< meta::decay<F> > d(F&& f) { return {std::forward<F>(f)}; }
  }

}
#endif
