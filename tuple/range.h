#ifndef MATH0X_TUPLE_RANGE_H
#define MATH0X_TUPLE_RANGE_H

#include <math0x/tuple/index.h>

// TODO this should be an external lib of its own

namespace math0x { 
  
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

}
#endif
