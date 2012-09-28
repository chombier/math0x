#ifndef MATH0X_TUPLE_FIND_H
#define MATH0X_TUPLE_FIND_H

#include <math0x/tuple/repeat.h>

// #include <group/log.h>

namespace math0x { 
  namespace tuple {

    namespace impl {

      template<int Start, int End, class U, class ... Args, class F>
      void find(const std::tuple<U, Args...>& ordered,
		U value,
		const F& f) {
	constexpr int Size = 1 + sizeof...(Args);

	// log("start:", Start,
	// 	  "end:", End,
	// 	  "size:", Size);
      
	static_assert( End <= Size, "derp !");
	static_assert( std::is_same< std::tuple<U, Args...>, tuple::repeat<U, Size> >::value,
		       "needs all same types");
      
	if( (Start + 1) == End) f.template operator()<Start>();
	else {
	  constexpr int Mid =  (Start + End) / 2;
	  
	  U mid = std::get< Mid >(ordered);
	
	  if( mid <= value ) find<Mid, End>(ordered, value, f);
	  else find<Start, Mid>(ordered, value, f);
	}
      
      }
    }
  
    template<class U, class ... Args, class F>
    void find( const std::tuple<U, Args...>& ordered,
	       U value,
	       const F& f) {
      return impl::find<0, sizeof...(Args) + 1>(ordered, value, f);
    }
  
  }

}
#endif
