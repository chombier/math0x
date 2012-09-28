#ifndef GROUP_COMMA_H
#define GROUP_COMMA_H

namespace math0x {

  // comma trick: ( f(), on_void(5.0) ).value is equal to f() when f
  // returns a double, or 5.0 when f returns void
  template<class T>
  struct comma {
    T value;
  };

  template<class T>
  comma<T> on_void(const T& value) {
    return { value };
  }
  
  template<class T>
  comma<T> operator,(const T& value, const comma<T>& ) {
    return {value};
  }
  
}


#endif
