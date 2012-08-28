#ifndef GROUP_TRAITS_H
#define GROUP_TRAITS_H

namespace euclid {

  template<class E> struct traits;

  template<class E>
  using field = typename traits<E>::field;

  template<class E>
  using dual = typename traits<E>::dual;
  
}

namespace lie {
  
  template<class G> struct traits;
  
  // Lie algebra
  template<class G>
  using alg = typename traits<G>::algebra;
  
  // Lie coalgebra
  template<class G>
  using coalg = euclid::dual< alg<G> >;
  
  // exponential: alg<G> -> G
  template<class G>
  using exp = typename traits<G>::exp;

  // logarithm: G -> alg<G>
  template<class G>
  using log = typename traits<G>::log;
  
  // adjoint: alg<G> -> alg<G>
  template<class G>
  using ad = typename traits<G>::ad;

  // coadjoint: coalg<G> -> coalg<G>
  template<class G>
  using adT = typename traits<G>::adT;
  
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


#endif
