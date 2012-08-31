#ifndef GROUP_SO_H
#define GROUP_SO_H

#include <group/types.h>
#include <group/lie.h>

#include <group/quaternion.h>

#include <group/vector.h>
#include <group/covector.h>

#include <group/error.h>

// 3-dimensional rotation group
template<class U>
class SO<3, U> {
  typedef quaternion<U> quat_type;
  quat_type quat;
public:
  typedef SO self;
  
  SO( const quat_type& quat = quat_type::Identity() ) 
    : quat( quat ) {
    assert( std::abs( quat.norm() -  1 ) < 1e-14 );
  }

  SO operator*(const SO& other) const {
    // TODO normalize ?
    return {quat * other.quat};
  }
  
  SO inv() const { 
    return { quat.conjugate() };
  }
  
  typedef vector<U, 3> vec3_type;

  vec3_type operator()(const vec3_type& x) const { 
    return quat * x;
  }
  
  // TODO push/pull
  
};


namespace lie {

  template<class U>
  struct traits< SO<3, U> > {

    typedef vector<U, 3> algebra;
    
    typedef SO<3, U> G;

    traits() { }
    traits(const G& ) { }

    G id() const { return {}; }
    G inv(const G& x) const { return x.inv(); }
    G prod(const G& x, const G& y) const { return x * y; }
    
    euclid::space<algebra> alg() const { return {}; }
    
    struct ad :  G {
      ad(const G& at) : G(at) { }
    };
    
    // TODO better ?
    struct adT { 
      G at;
      
      adT( const G& g ) : at( g.inv() ) { }
      
      lie::coalg<G> operator()(const lie::coalg<G>& f) const {
      	return at(f.transpose()).transpose();
      }
      
    };


    struct exp {
      typedef exp self;
      
      exp(const group<G>& = group<G>() ) { }

      G operator()(const lie::alg<G>& x) const { 
	// throw error("not implemented");
      }

      

    };


    struct log {
      typedef log self;
      
      log(const group<G>& = group<G>() ) { }

      lie::alg<G> operator()(const G& x) const { 
      	// throw error("not implemented");
      }
      
      // struct push { push(const log&, const G&) { } };
      // struct pull { pull( const log&, const G& ) { } }; 
      
    };
    
    
  };


}

#endif
