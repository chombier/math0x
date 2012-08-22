#ifndef GROUP_EIGEN_H
#define GROUP_EIGEN_H

#include <group/euclid.h>

#include <Eigen/Core>

// adapts eigen vectors for use with this library.
namespace euclid {

  template<class U, int M, int N>
  struct traits< Eigen::Matrix<U, M, N> > {

    static_assert( (M == 1) || (N == 1),
		   "sorry, unimplemented: at least one dimension must be 1");
    
    NN n;
    space<U> sub;
    
    template<int K>
    static NN static_size() {
      static_assert( K != Eigen::Dynamic, "should not be called for dynamic size" );
      return K;
    }

    traits(NN n = static_size<M>() * static_size<N>(),
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
      res = E::Zero( dim() );
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


};


#endif
