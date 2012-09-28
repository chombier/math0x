#ifndef MATH0X_LOG_H
#define MATH0X_LOG_H

#include <iostream>

void debug() { std::cout << std::endl; }

template<class A, class ... Args>
void debug(const A& a, const Args&... args) {
  std::cout << a << " ";
  debug( args... );
}

#endif
