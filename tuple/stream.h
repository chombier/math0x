#ifndef GROUP_TUPLE_STREAM_H
#define GROUP_TUPLE_STREAM_H

#include <group/tuple/range.h>
#include <ostream>

namespace tuple {

  namespace impl {

    template<class ... Args>
    struct stream {
      
      std::ostream& out;
      const std::tuple<Args...>& data;
    
      template<int I>
      void operator()() const {
	out << std::get<I>(data);
	
	if( I < (sizeof...(Args) - 1) ) out << ", ";
      }
      
    };
  }
}


template<class ... Args>
std::ostream& operator<<(std::ostream& out,
			 const std::tuple<Args...>& args) {
  tuple::make_range(args).apply( tuple::impl::stream<Args...>{out, args} ); 
  return out;
}

#endif
