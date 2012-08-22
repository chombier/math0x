#ifndef GROUP_TUPLE_RANGE_H
#define GROUP_TUPLE_RANGE_H

#include <group/tuple/index.h>

namespace tuple {

  namespace impl {
    template<int I>
    struct range {
      typedef typename range<I-1>::type::next_type type;
    };


    template<>
    struct range<0> {
      typedef index<> type;
    };
  }

  template<class ... Args>
  using range = typename impl::range< sizeof...(Args) >::type;
  
  template<class ...Args>
  range<Args...> make_range(const std::tuple<Args...>& ) { return {}; }
  
}


#endif
