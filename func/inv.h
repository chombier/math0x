#ifndef GROUP_FUNC_INV_H
#define GROUP_FUNC_INV_H

#include <group/lie.h>

#include <group/func/minus.h>

namespace func {

  // group inverse
  template<class G>
  struct inverse {
    typedef inverse self;
    
    lie::group<G> group;
    inverse(const lie::group<G>& group = lie::group<G>() ) : group(group) { }


    G operator()(const G& g) const { return group.inv(g); } 

    struct push : comp< minus< lie::alg<G> >, lie::ad<G> > {

      push(const inverse& of,
	   const G& at) : push::self( minus< lie::alg<G> >{of.group.alg()},
				      lie::ad<G>{at} ) {

      }

    };

    struct pull : comp< minus< lie::coalg<G> >, lie::adT<G> > {
       
      pull(const inverse& of,
	   const G& at) : pull::self( minus< lie::coalg<G> >{*of.group.alg()},
				      lie::adT<G>{at} ) {
	 
      }
       
    };

  };

}


#endif
