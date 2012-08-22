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
 
  template<class G>
  using algebra = typename traits<G>::algebra;
  
  template<class G>
  using coalgebra = euclid::dual< algebra<G> >;
  
  template<class G>
  using exponential = typename traits<G>::exponential;

  template<class G>
  using logarithm = typename traits<G>::logarithm;
  
  template<class G>
  using adjoint = typename traits<G>::adjoint;

  template<class G>
  using coadjoint = typename traits<G>::coadjoint;

}

namespace func {
  template<class F, class = void> struct traits;
}


#endif
