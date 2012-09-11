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
    typedef subgroup self;
    
    typedef euclid::field< lie::alg< G > > field;
    
    lie::alg< G > dir;
    lie::group<G> group;
    euclid::space< lie::alg<G> > alg;
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
    struct push : line< lie::alg<G> > {
      

      push(const subgroup& of, const field& )
	: push::self( of.dir ) { 

      }

    };

    struct pull : form< lie::coalg<G> > {
      
      pull(const subgroup& of, const field& )
	: pull::self( of.dir ) { 
	
      }

    };

    


  };



}


#endif
