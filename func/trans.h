#ifndef GROUP_FUNC_TRANS_H
#define GROUP_FUNC_TRANS_H

#include <group/lie.h>
#include <group/func/id.h>

namespace func {
  
  // lie group left/right translations
  

  template<class G>
  struct left {
    typedef left self;

    G h;
    lie::group<G> group;
    
    left(const G& h) : h(h), group(h) { }

    G operator()(const G& g) const { return group.prod(h, g); }
    
    struct push : id< lie::alg<G> > { 
      push(const left&, const G&) { } 
    };


    struct pull : id< lie::coalg<G> > { 
      pull(const left&, const G&) { } 
    };
    
  };

  template<class G>
  left< G > L(const G& g) { return {g}; }


  template<class G>
  struct right {
    typedef right self;

    G h;
    lie::group<G> group;

    right(const G& h) : h(h), group(h) { }
    
    G operator()(const G& g) const {  return group.prod(g, h); }
    
    struct push : lie::ad<G> { 
      push(const right& of, const G&) : pull::self(of.group.inv(of.h)) { } 
    };
    
    struct pull : lie::adT<G> {
      pull(const right& of, const G&) : pull::self(of.group.inv(of.h)) { } 
    };
    

  };

  template<class G>
  right< G > R(const G& g) { return {g}; }
  

}


#endif
