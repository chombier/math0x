#ifndef MATH0X_META_H
#define MATH0X_META_H

#include <utility>

// meta-programming
namespace meta {
  template<class F>
  using decay = typename std::decay<F>::type;

  template<class F>
  using remove_pointer = typename std::remove_pointer<F>::type;
  
  template<class ... Args>
  inline void noop(Args&& ... ) { }

}

#endif
