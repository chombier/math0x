#ifndef GROUP_FUNC_SUBGROUP_H
#define GROUP_FUNC_SUBGROUP_H

#include <group/lie.h>

#include <group/func/line.h>
#include <group/func/form.h>

namespace func {

  // one parameter subgroup: t -> exp( t * dir )

  // we don't simply use exp << line(dir) as the push/pull are simpler
  template<class G>
  struct subgroup {
    typedef subgroup base;
    
    typedef euclid::field< lie::algebra< G > > field;
    
    lie::algebra< G > dir;
    lie::group<G> group;
    euclid::space< lie::algebra<G> > alg;
    lie::exp< G > exp;
    
    subgroup(const G& dir) 
    : dir(dir), 
      group(dir),
      alg(group.alg()),
      exp( group.exp() )
    { }
    
    
    G operator()(const field& t ) const {
      return exp( alg.scal(t, dir) );
    }
    
    // mmmhh not sure if this will produce good second derivatives...
    struct push : line< lie::algebra<G> > {
      

      push(const subgroup& of, const field& )
	: push::base( of.dir ) { 

      }

    };

    struct pull : form< lie::coalgebra<G> > {
      
      pull(const subgroup& of, const field& )
	: pull::base( of.dir ) { 
	
      }

    };

    


  };



}


#endif
