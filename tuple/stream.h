#ifndef MATH0X_TUPLE_STREAM_H
#define MATH0X_TUPLE_STREAM_H

#include <math0x/tuple/range.h>
#include <ostream>

template<class ... Args>
std::ostream& operator<<(std::ostream& out,
			 const std::tuple<Args...>& args);
namespace math0x { 
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
}

template<class ... Args>
std::ostream& operator<<(std::ostream& out,
			 const std::tuple<Args...>& args) {
  math0x::tuple::make_range(args).apply( math0x::tuple::impl::stream<Args...>{out, args} ); 
  return out;
}

#endif
