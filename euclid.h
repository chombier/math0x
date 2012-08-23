#ifndef GROUP_EUCLID_H
#define GROUP_EUCLID_H

#include <group/types.h>
#include <group/meta.h>

namespace euclid {
  
  template<class E>
  class space {
    typedef traits<E> impl_type;
    traits<E> impl;

    void constraints() {
      // TODO more
      
      E x = zero();
      space euclid(x);

      meta::noop(x, euclid);
    }    
    
  public:

    ~space() {
      meta::noop( &space::constraints );
    }
      
    template<class...Args>
    space(Args&& ... args) 
      : impl(std::forward<Args>(args)...) { }
    
    space(const space& ) = default;
    space(space&& ) = default;
    
    // dual geometry
    space< dual<E> > operator*() const { 
      return *impl;
    }

    // coordinate accessors 
    field<E>& coord(NN i, E& e) const { 
      return impl.coord(i, e); 
    }
    
    const field<E>& coord(NN i, const E& e) const { 
      return impl.coord(i, e); 
    }
    
    // space dimension
    NN dim() const { return impl.dim(); }

    // zero element
    E zero() const { return impl.zero(); }
    
    // coordinate iterators
    template<class F>
    void each(const E& x, const F& f) const {
      for(natural i = 0, n = dim(); i < n; ++i) {
	f( coord(i, x) );
      }
    };

    template<class F>
    void each(E& x, const F& f) const {
      for(natural i = 0, n = dim(); i < n; ++i) {
	f( coord(i, x) );
      }
    };


    // vector space operations
    E minus(E&& x) const {
      each(x, [&](field<E>& xi) {
	  xi = -xi;
	});
      
      return std::move(x);
    } 
    

    E minus(const E& res) const {
      return minus( E(res) );
    } 

    E scal(field<E> lambda, E&& x) const {
      each(x, [&](field<E>& xi) {
	  xi = lambda * xi;
	});
      return std::move(x);
    }
    
    E scal(field<E> lambda, const E& res) const {
      return scal(lambda, E(res) );
    }
    

    E sum(E&& x, const E& y) const {
      NN i = 0;

      each(x, [&](field<E>& xi) {
	  xi += this->coord(i, y);
	  ++i;
	});
      
      return std::move(x);
    }

    E sum(const E& x, E&& y) const {
      return sum(std::move(y), x);
    }

    E sum(E&& x, E&& y) const {
      return sum(std::move(x), y);
    }

    E sum(const E& x, const E& y) const {
      return sum(E(x), y);
    }
    
  };
  
}


#endif
