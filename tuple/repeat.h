#ifndef GROUP_TUPLE_REPEAT_H
#define GROUP_TUPLE_REPEAT_H

#include <group/tuple/push.h>

namespace tuple {
  
  namespace impl {
    template<class A, int M>
    struct repeat {
      typedef push::back< typename repeat<A, M-1>::type, A > type;
    };
    
    template<class A>
    struct repeat<A, 0> {
      typedef std::tuple<> type;
    };
    
  }
  
  template<class A, int M>
  using repeat = typename impl::repeat<A, M>::type;
  
}


#endif
