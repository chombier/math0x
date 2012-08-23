#ifndef GROUP_EIGEN_H
#define GROUP_EIGEN_H

#include <group/euclid.h>

#include <Eigen/Core>

// adapts eigen vectors for use with this library.

namespace impl {
  template<int K>
  NN static_size() {
    static_assert( K != Eigen::Dynamic, "should not be called for dynamic size" );
    return K;
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

    typedef Eigen::Matrix< lie::algebra<U>, M, N > algebra;

    typedef Eigen::Matrix<U, M, N> G;
    
    traits( const G& x ) 
    : n( x.size() ),
      sub( x(0)) {
      assert( n );
    }
    
    
    G id() const {
      // TODO lol 
    }

    G inv( const G& x ) const { 
      // TODO lol
    }

    G prod( const G& x, const G& y) const { 
      // TODO lol
    }

  };

}
#endif
