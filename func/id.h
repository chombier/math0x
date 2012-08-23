#ifndef GROUP_FUNC_ID_H
#define GROUP_FUNC_ID_H

#include <group/types.h>

namespace func {

  // identity over a group G
  template<class G>
  struct id {
    typedef id self;
    
    // TODO this should be safe since rvalue references bind below ?
    const G& operator()(const G& x) const { return x; }
    
    G operator()(G&& x) const { return std::move(x); }
    
    struct push : id< lie::algebra<G> > { 
      push( const id&, const G& )  { }
    };

    struct pull : id< lie::coalgebra<G> > { 
      pull( const id&, const G& )  { }
    };
    
  };


}


#endif
