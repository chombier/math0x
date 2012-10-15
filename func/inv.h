#ifndef MATH0X_FUNC_INV_H
#define MATH0X_FUNC_INV_H

#include <math0x/lie.h>

#include <math0x/func/minus.h>

namespace math0x { 
	namespace func {

		// lie group inverse
		template<class G>
		struct inv {
			typedef inv base;
    
			lie::group<G> group;
			inv(const lie::group<G>& group = lie::group<G>() ) : group(group) { }


			G operator()(const G& g) const { return group.inv(g); } 

			struct push : comp< minus< lie::algebra<G> >, lie::Ad<G> > {

				push(const inv& of,
				     const G& at) : push::base( minus< lie::algebra<G> >{of.group.alg()},
				                                lie::Ad<G>{at} ) {

				}

			};

			struct pull : comp< minus< lie::coalgebra<G> >, lie::AdT<G> > {
       
				pull(const inv& of,
				     const G& at) : pull::base( minus< lie::coalgebra<G> >{of.group.coalg()},
				                                lie::AdT<G>{at} ) {
	 
				}
       
			};

		};


		template<class G>
		inv<G> make_inv(const lie::group<G>& group) { return {group}; }

	}

}
#endif
