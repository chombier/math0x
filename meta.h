#ifndef MATH0X_META_H
#define MATH0X_META_H

#include <utility>

// meta-programming
namespace meta {

	// some short-hands because we don't like writing typename all the
	// time
  template<class F>
  using decay = typename std::decay<F>::type;

  template<class F>
  using remove_pointer = typename std::remove_pointer<F>::type;
  
	// does nothing :-) use it to remove unused variables warnings
  template<class ... Args>
  inline void noop(Args&& ... ) { }
	
}

#endif
