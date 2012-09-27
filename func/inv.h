#ifndef GROUP_FUNC_INV_H
#define GROUP_FUNC_INV_H

#include <group/lie.h>

#include <group/func/minus.h>

namespace func {

  // lie group inverse
  template<class G>
  struct inverse {
    typedef inverse base;
    
    lie::group<G> group;
    inverse(const lie::group<G>& group = lie::group<G>() ) : group(group) { }


    G operator()(const G& g) const { return group.inv(g); } 

    struct push : comp< minus< lie::algebra<G> >, lie::Ad<G> > {

      push(const inverse& of,
	   const G& at) : push::base( minus< lie::algebra<G> >{of.group.alg()},
				      lie::Ad<G>{at} ) {

      }

    };

    struct pull : comp< minus< lie::coalgebra<G> >, lie::AdT<G> > {
       
      pull(const inverse& of,
	   const G& at) : pull::base( minus< lie::coalgebra<G> >{of.group.coalg()},
				      lie::AdT<G>{at} ) {
	 
      }
       
    };

  };

}


#endif
