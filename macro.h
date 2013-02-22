#ifndef MATH0X_MACRO_H
#define MATH0X_MACRO_H

#define macro_returns(...)	  \
	decltype( __VA_ARGS__ ) { \
		return __VA_ARGS__; \
  }			       \
	

#endif
