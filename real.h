#ifndef GROUP_REAL_H
#define GROUP_REAL_H

#include <group/euclid.h>
#include <cassert>

namespace euclid {

  template<>
  struct traits<RR> {
    typedef RR field;
    typedef RR dual;

    traits(const RR& ) { }
    traits() { }
    
    NN dim() const { return 1; }

    RR zero() const { return 0; }

    field& coord(NN i, RR& x) const { 
      assert( i == 0 );
      return x; 
    }
    const field& coord(NN i, const RR& x) const { 
      assert( i == 0 );
      return x; 
    }
    
  };


}


#endif
