#ifndef GROUP_LIE_H
#define GROUP_LIE_H

#include <group/euclid.h>

namespace lie {
  
  template<class G>
  class group {
    typedef traits<G> impl_type;
    impl_type impl;
    
    void constraints() {
      // TODO more

      G x = id();
      group lie( x );
      
      // TODO check exp/log domain/range

      noop(x, lie);
    }
    
  public:
    
    ~group() {
      noop( &group::constraints );
    }

    // identity
    G id() const { return impl.id(); }

    // inverse
    G inv(const G& x) const { return impl.inv(x); }

    // product
    G prod(const G& x, const G& y) const { 
      return impl.prod(x, y); 
    }

    // TODO rvalue overloads for efficiency ?
    
    // algebra euclidean structure
    euclid::space< algebra > alg() const { 
      return impl.alg();
    }


    adjoint<G> ad(const G& g) const { return {g}; }
    coadjoint<G> adT(const G& g) const { return {g}; }

    exponential<G> exp() const { return {*this}; }
    logarithm<G> log() const { return {*this}; }
    
    template<class ... Args>
    group(Args&& ... args) : impl(std::forward<Args>(args)...) {  }
    
    group(const group& ) = default;
    group(group&& ) = default;
    
  };


}


#endif
