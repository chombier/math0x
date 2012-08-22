#ifndef GROUP_LOG_H
#define GROUP_LOG_H

#include <iostream>

void log() { std::cout << std::endl; }

template<class A, class ... Args>
void log(const A& a, const Args&... args) {
  std::cout << a << " ";
  log( args... );
}

#endif
