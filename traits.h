#ifndef MATH0X_TRAITS_H
#define MATH0X_TRAITS_H

namespace math0x {
  namespace euclid {

    template<class E> struct traits;

    // field of the vector space
    template<class E>
    using field = typename traits<E>::field;

    // vector dual space
    template<class E>
    using dual = typename traits<E>::dual;

	  template<class E>
	  using range = typename traits<E>::range;
	  
  }

  namespace lie {
  
    template<class G> struct traits;
  
    // Lie algebra
    template<class G>
    using algebra = typename traits<G>::algebra;
  
    // Lie coalgebra
    template<class G>
    using coalgebra = euclid::dual< algebra<G> >;
  
    // exponential: alg<G> -> G
    template<class G>
    using exp = typename traits<G>::exp;

    // logarithm: G -> alg<G>
    template<class G>
    using log = typename traits<G>::log;
  
    // group adjoint: alg<G> -> alg<G>
    template<class G>
    using Ad = typename traits<G>::Ad;

    // group adjoint transpose: coalg<G> -> coalg<G>
    template<class G>
    using AdT = typename traits<G>::AdT;

	  // // algebra adjoint: alg<G> -> alg<G>
	  // template<class G>
	  // using ad = typename traits<G>::ad;

	  // // algebra adjoint transpose: coalg<G> -> coalg<G>
	  // template<class G>
	  // using adT = typename traits<G>::adT;
	  
  }

  namespace func {
  
    template<class F, class = void> struct traits;

    // domain over which the function is defined
    template<class F>
    using domain = typename traits<F>::domain;
  
    // range of the function
    template<class F>
    using range = typename traits<F>::range;
  
    // pushforward (derivative) type
    template<class F>
    using push = typename traits<F>::push;
  
    // pullback (derivative transpose) type
    template<class F>
    using pull = typename traits<F>::pull;

  }
}

#endif
