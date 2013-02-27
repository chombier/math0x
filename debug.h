#ifndef MATH0X_DEBUG_H
#define MATH0X_DEBUG_H

#warning debug header included lol
#include <iostream>

namespace {
	// TODO bloat ?
	inline void debug() { std::cout << std::endl; }
	
	template<class A, class ... Args>
	void debug(const A& a, const Args&... args) {
		std::cout << a << " ";
		debug( args... );
	}
}

#endif
