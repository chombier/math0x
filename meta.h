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
