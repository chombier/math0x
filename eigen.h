#ifndef GROUP_EIGEN_H
#define GROUP_EIGEN_H

#include <group/euclid.h>
#include <group/lie.h>
#include <group/error.h>

#include <Eigen/Core>

// adapts eigen vectors for use with this library.

namespace impl {
  template<int K>
  NN static_size() {
    static_assert( K != Eigen::Dynamic, "should not be called for dynamic size" );
    return K;
  }
  
  template<class Xpr, class F>
  Eigen::CwiseUnaryOp<F, Xpr > unary(const Xpr& xpr, const F& f) {
    return {xpr, f};
  }

  template<class LHS, class RHS, class F>
  Eigen::CwiseBinaryOp<F, LHS, RHS > binary(const LHS& lhs,
					    const RHS& rhs,
					    const F& f) {
    return {lhs, rhs, f};
  }
  
  
}
namespace euclid {

  template<class U, int M, int N>
  struct traits< Eigen::Matrix<U, M, N> > {

    static_assert( (M == 1) || (N == 1),
		   "sorry, unimplemented: at least one dimension must be 1");
    
    NN n;
    space<U> sub;
    
    traits(NN n = impl::static_size<M>() * impl::static_size<N>(),
	   const space<U>& sub = space<U>())
      : n(n),
	sub(sub) {
      assert( dim() );
    }
    
    typedef euclid::field<U> field;
    
    typedef Eigen::Matrix<U, M, N> E;
    
    traits( const E& x ) 
    : n( x.size() ),
      sub( x(0)) {
      assert( dim() );
    }
    
    // dual space: swaps dimensions
    typedef Eigen::Matrix<U, N, M> dual;

    NN dim() const { return n * sub.dim(); }
    
    E zero() const {
      E res;
      res.resize( n );

      for( NN i = 0; i < n; ++i) {
	res(i) = sub.zero();
      }
      
      return res;
    }
    
    field& coord(NN i, E& x) const {
      return sub.coord(i % sub.dim(), x( i / sub.dim() ) );
    }
    
    const field& coord(NN i, const E& x) const {
      return sub.coord(i % sub.dim(), x( i / sub.dim() ) );
    }
    
    space< dual > operator*() const {
      return { n, *sub };
    }
    
  };

  
}


namespace lie {
  
  // TODO termwise operations ?
  template<class U, int M, int N>
  struct traits< Eigen::Matrix<U, M, N> > {

    static_assert( (M == 1) || (N == 1),
		   "sorry, unimplemented: at least one dimension must be 1");
    
    NN n;
    group<U> sub;

    traits(NN n = impl::static_size<M>() * impl::static_size<N>(),
	   const group<U>& sub = group<U>())
      : n(n),
	sub(sub) {
      assert( n );
    }

    typedef Eigen::Matrix< lie::alg<U>, M, N > algebra;

    typedef Eigen::Matrix<U, M, N> G;
    
    traits( const G& x ) 
    : n( x.size() ),
      sub( x(0)) {
      assert( n );
    }
    
    struct op {

      struct inv {
	const group<U>& sub;

	U operator()(const U& x) const {
	  return sub.inv(x);
	}
	
      };

      struct prod {
	const group<U>& sub;
	
	U operator()(const U& x, const U& y) const {
	  return sub.prod(x, y);
	}
	
      };

    };


    G id() const {
      G res; res.resize(n);
      for(NN i = 0, n = res.size(); i < n; ++i) res(i) = sub.id(); 
      return res;
    }
    
    G inv( const G& x ) const { 
      return impl::unary(x, typename op::inv{sub} );
    }

    G prod( const G& x, const G& y) const { 
      return impl::binary(x, y, typename op::prod{sub} );
    }

    euclid::space< lie::alg<G> > alg() const { return { n, sub.alg() }; }

    struct ad {
      ad( const G& ) { };
      
      lie::alg<G> operator()(const lie::alg<G>& ) const {
	throw error("not implemented");
      }
      
    };
    
    struct adT {  
      adT( const G& ) { };
      
      lie::coalg<G> operator()(const lie::coalg<G>& ) const {
	throw error("not implemented");
      }
    };

    struct exp {
      exp( const group<G>& ) { }

      G operator()(const lie::alg<G>& x) const { 
	throw error("not implemented");
      }
      
    };

    struct log { 
      log( const group<G>& ) { }

      lie::alg<G> operator()(const G& x) const { 
	throw error("not implemented");
      }
      
    };
    
  };

}
#endif
