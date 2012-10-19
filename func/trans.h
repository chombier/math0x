#ifndef MATH0X_FUNC_TRANS_H
#define MATH0X_FUNC_TRANS_H

#include <math0x/lie.h>
#include <math0x/func/id.h>

namespace math0x { 
	namespace func {
  
		// lie group left/right translations
		template<class G>
		struct left {
			typedef left base;

			G h;
			lie::group<G> group;
    
			left(const G& h) : h(h), group(h) { }

			G operator()(const G& g) const { return group.prod(h, g); }
    
			struct push : id< lie::algebra<G> > { 
				push(const left&, const G&) { } 
			};


			struct pull : id< lie::coalgebra<G> > { 
				pull(const left&, const G&) { } 
			};
    
		};

		template<class G>
		left< G > L(const G& g) { return {g}; }


		template<class G>
		struct right {
			typedef right base;

			G h;
			lie::group<G> group;

			right(const G& h) : h(h), group(h) { }
    
			G operator()(const G& g) const {  return group.prod(g, h); }
    
			struct push : lie::Ad<G> { 
				push(const right& of, const G&) : push::Ad(of.group.inv(of.h)) { } 
			};
    
			struct pull : lie::AdT<G> {
				pull(const right& of, const G&) : pull::AdT(of.group.inv(of.h)) { } 
			};
    

		};

		template<class G>
		right< G > R(const G& g) { return {g}; }
  

	}

}
#endif
