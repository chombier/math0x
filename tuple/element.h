#ifndef MATH0X_TUPLE_ELEMENT_H
#define MATH0X_TUPLE_ELEMENT_H

#include <tuple>

namespace math0x { 
  namespace tuple {

    template<int M, class Tuple>
    using element = typename std::tuple_element<M, Tuple>::type;
  
  }

}
#endif
