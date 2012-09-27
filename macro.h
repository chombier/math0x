#ifndef GROUP_MACRO_H
#define GROUP_MACRO_H

#define macro_auto(...)	       \
  decltype( __VA_ARGS__ ) {    \
  return __VA_ARGS__;	       \
  }			       \
  

#endif
