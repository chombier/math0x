#ifndef MATH0X_SO3_H
#define MATH0X_SO3_H

#include <math0x/types.h>
#include <math0x/lie.h>

#include <math0x/quaternion.h>

#include <math0x/vector.h>
#include <math0x/covector.h>

#include <math0x/error.h>

namespace math0x { 

  // 3-dimensional rotation group
  template<class U>
  class SO<3, U> {
    typedef quaternion<U> quat_type;
    quat_type quat;
  public:
    typedef SO base;
  
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
  
    struct push : SO {
      
      push(const SO& of, const vec3_type& ) 
	: push::base(of) {

      }

    };

    struct pull {
      
      SO of_inv;
      
      pull(const SO& of, const vec3_type& ) 
	: of_inv(of.inv()) {
	
      }
      
      lie::coalgebra<vec3_type> operator()(const lie::coalgebra<vec3_type>& f) const {
	return of_inv(f.transpose()).transpose();
      }
      
    };
  

  };


  namespace lie {

    template<class U>
    struct traits< SO<3, U> > {

      // TODO make so<3, U> and wrap skew-symmetric matrices
      typedef vector<U, 3> algebra;
    
      typedef SO<3, U> G;

      traits() { }
      traits(const G& ) { }

      G id() const { return {}; }
      G inv(const G& x) const { return x.inv(); }
      G prod(const G& x, const G& y) const { return x * y; }
    
      euclid::space<algebra> alg() const { return {}; }
    
      
      algebra bracket(const algebra& x, const algebra& y) const {
	return x.cross(y);
      }


      struct Ad :  G {
	Ad(const G& at) : G(at) { }
      };
    
      // TODO better ?
      struct AdT { 
	typedef AdT base;
      
	G at;
      
	AdT( const G& g ) : at( g.inv() ) { }
      
	lie::coalgebra<G> operator()(const lie::coalgebra<G>& f) const {
	  return at(f.transpose()).transpose();
	}
      
      };


      struct exp {
	typedef exp base;
      
	exp(const group<G>& = group<G>() ) { }

	G operator()(const lie::algebra<G>& ) const { 
	  throw error("not implemented");
	}

      

      };


      struct log {
	typedef log base;
      
	log(const group<G>& = group<G>() ) { }
      
	lie::algebra<G> operator()(const G& ) const { 
	  throw error("not implemented");
	}
      
	// struct push { push(const log&, const G&) { } };
	// struct pull { pull( const log&, const G& ) { } }; 
      
      };
    
    
    };


  }


  namespace func {

    template<class U>
    struct apply< SO<3, U> > {
      
      typedef SO<3, U> G;

      typedef std::tuple<G, func::domain<G> > domain;
      typedef func::range<G> range;
      
      range operator()(const domain& x) const {
	return std::get<0>(x)(std::get<1>(x));
      }
      

      struct push {

	domain at;
	lie::group<G> group;
	
	push( const apply&, const domain& at) : at(at) { }

	lie::algebra<range> operator()(const lie::algebra<domain>& dx) const {
	  
	  return std::get<0>(at)( group.bracket(std::get<0>(dx), std::get<1>(at))
				  + std::get<1>(dx) );
	}
	
      };
      

    };

  }


}



#endif
