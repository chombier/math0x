#ifndef GROUP_LIE_H
#define GROUP_LIE_H

#include <group/euclid.h>
#include <group/meta.h>

namespace lie {
  
  template<class G>
  class group {
    // friend class traits<G>;
    typedef traits<G> impl_type;

  public:
    impl_type impl;
  private:
    void constraints() {
      // TODO more

      G x = id();
      group lie( x );
      
      // static_assert( std::is_same< G, func::range< lie::exp<G> > >::value, "exp range error" );
      // static_assert( std::is_same< lie::alg<G>, func::domain< lie::exp<G> > >::value, "exp domain error" );

      // static_assert( std::is_same< G, func::domain< lie::log<G> > >::value, "log range error" );
      // static_assert( std::is_same< lie::alg<G>, func::range< lie::log<G> > >::value, "exp domain error" );

      // static_assert( std::is_same< lie::alg<G>, func::range< lie::ad<G> > >::value, "ad range error" );
      // static_assert( std::is_same< lie::alg<G>, func::domain< lie::ad<G> > >::value, "ad domain error" );

      // static_assert( std::is_same< lie::coalg<G>, func::range< lie::adT<G> > >::value, "adT range error" );
      // static_assert( std::is_same< lie::coalg<G>, func::domain< lie::adT<G> > >::value, "adT domain error" );
      
      meta::noop(x, lie);
    }
    
  public:
    
    ~group() {
      meta::noop( &group::constraints );
    }

    // identity
    G id() const { 
      return impl.id(); 
    }

    // inverse
    G inv(const G& x) const { 
      return impl.inv(x); 
    }
    
    // product
    G prod(const G& x, const G& y) const { 
      return impl.prod(x, y); 
    }

    // TODO rvalue overloads for efficiency ?
    
    // algebra euclidean structure
    euclid::space< lie::algebra<G> > alg() const { 
      return impl.alg();
    }
    
    // adjoint
    lie::Ad<G> Ad(const G& g) const { return {g}; }

    // coadjoint
    lie::AdT<G> AdT(const G& g) const { return {g}; }
    
    // exponential
    lie::exp<G> exp() const { return {*this}; }

    // logarithm
    lie::log<G> log() const { return {*this}; }
    
    // lie bracket
    algebra<G> bracket(const algebra<G>& x, const algebra<G>& y) { 
      return impl.bracket(x, y);
    }

    // lie cobracket
    coalgebra<G> cobracket(const coalgebra<G>& x, const coalgebra<G>& y) { 
      return impl.cobracket(x, y);
    }
    
    // forward ctor
    template<class ... Args>
    group(Args&& ... args) : impl(std::forward<Args>(args)...) {  }
    
    // 
    group(const group& ) = default;
    group(group&& ) = default;
    
  };


}


#endif
