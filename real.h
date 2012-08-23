#ifndef GROUP_REAL_H
#define GROUP_REAL_H

#include <group/types.h>
#include <group/func/id.h>

#include <cassert>

// real numbers, as defined in types.h

// euclidean structure
namespace euclid {

  template<>
  struct traits<RR> {
    
    typedef RR field;
    typedef RR dual;

    traits(const RR& ) { }
    traits() { }
    
    NN dim() const { return 1; }

    RR zero() const { return 0; }

    space<dual> operator*() const { return {}; } 

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

// (flat) lie group structure
namespace lie {
  
  template<> 
  struct traits<RR> {
    
    typedef RR algebra;
    
    struct adjoint : func::id< RR > {
      adjoint(const RR& ) { }
    };

    struct coadjoint : func::id< RR > {
      coadjoint(const RR& ) { }
    };
    
    struct exponential : func::id<RR> {
      exponential(const group<RR>& ) { }
    };

    struct logarithm : func::id<RR> {
      logarithm(const group<RR>& ) { }
    };
    
    RR id() const { return 0; }
    RR inv(const RR& x) const { return -x; }
    RR prod(const RR& x, const RR& y) const { return x + y; }
    
    traits() { }
    traits(const RR& ) { }
    
  };

}



#endif
