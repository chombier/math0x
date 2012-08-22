#ifndef GROUP_META_H
#define GROUP_META_H

#include <utility>

// meta-programming
namespace meta {
  template<class F>
  using decay = typename std::decay<F>::type;

  template<class F>
  using remove_pointer = typename std::remove_pointer<F>::type;
  
  template<class ... Args>
  void noop(Args&& ... ) { }
}

#endif
