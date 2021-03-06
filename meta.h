#ifndef MATH0X_META_H
#define MATH0X_META_H

#include <utility>

#include <string>

#ifdef __GNUC__
#include <cxxabi.h>
#endif

namespace math0x {
	// meta-programming
	namespace meta {

		// some short-hands because we don't like writing typename all the
		// time
		template<class F>
		using decay = typename std::decay<F>::type;

		template<class F>
		using remove_pointer = typename std::remove_pointer<F>::type;
  
		template<bool b, class T = void>
		using enable_if = typename std::enable_if<b, T>::type;
		
		
		namespace impl { 
			template<class, class T>
			struct second { 
				typedef T type;
			};
		}

		template<class T1, class T2> using second = typename impl::second<T1, T2>::type;
		
		
		// does nothing :-) use it to remove unused variables warnings
		template<class ... Args>
		inline void noop(Args&& ... ) { }
		
		

		// overload priority control
		template<unsigned> struct priority;
		
		// lowest priority
		template<> struct priority<0> { };

		// higher priority = more specialized types
		template<unsigned I> struct priority : priority<I - 1> { };
		

		// address
		template<class F>
		decay<F>* addr(F&& f) { return &f; }
		

		namespace impl {

			template< template<class ...> class temp, 
			          class Args > 
			struct unpack_args; 
		   
			template< template<class ...> class temp,
			          class ... Args > 
			struct unpack_args<temp, std::tuple<Args...> > {
				typedef unpack_args base;
				typedef temp<Args...> type;
			};		
			
		}

		// use this to unpack variadic template args from a std::tuple type
		template< template<class...> class T,
		          class Args > 
		using unpack_args = typename impl::unpack_args<T, Args>::type;
	
		// type name
		template<class F>
		std::string name() {
#ifdef __GNUC__
			int     status;
			char   *realname;
			
			realname = abi::__cxa_demangle(typeid(F).name(), 0, 0, &status);
			std::string res( realname );
			std::free(realname);
			
			return res;
#else
			return typeid(F).name();
#endif
			
		}
		
	}
}
#endif
