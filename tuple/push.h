#ifndef MATH0XS_TUPLE_PUSH_H
#define MATH0XS_TUPLE_PUSH_H

#include <tuple>

namespace math0x { 
  namespace tuple {

    namespace push {

      namespace impl {
      
	template<class, class> struct back;
	template<class, class> struct front;

	template<class ... Args, class A>
	struct back< std::tuple<Args...>, A > {
	  typedef std::tuple<Args..., A> type;
	};

	template<class ... Args, class A>
	struct front< std::tuple<Args...>, A > {
	  typedef std::tuple<A, Args...> type;
	};
      
      }


      template<class Tuple, class A>
      using back = typename impl::back<Tuple, A>::type;

      template<class Tuple, class A>
      using front = typename impl::front<Tuple, A>::type;
    
    }

  }

}
#endif
