#ifndef GROUP_TUPLE_ELEMENT_H
#define GROUP_TUPLE_ELEMENT_H

#include <tuple>

namespace tuple {

  template<int M, class Tuple>
  using element = typename std::tuple_element<M, Tuple>::type;
  
}


#endif
