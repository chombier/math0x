#ifndef GROUP_FUNC_PROD_H
#define GROUP_FUNC_PROD_H

#include <group/lie.h>

namespace func {

  // lie group product
  template<class G>
  struct product {
    lie::group<G> group;
    
    product(const lie::group<G>& group = {} ) : group(group) { }
    
    typedef std::tuple<G, G> domain;
    
    G operator()(const domain& g) const {
      return group.prod( std::get<0>(g), 
			 std::get<1>(g) );
    }
    
    
    struct push {

      euclid::space< lie::algebra<G> > alg;
      lie::Ad<G> Ad;

      push(const product& of, const domain& at)
	: alg(of.group.alg()), 
	  Ad(of.group.ad( of.group.inv(std::get<1>(at) ) ) )
      {
	
      }

      lie::algebra<G> operator()(const lie::algebra<domain>& v) const {
	return alg.sum( Ad(std::get<0>(v)), 
			std::get<1>(v) );
      }
      
    };


    struct pull {

      lie::AdT<G> AdT;

      pull(const product& of, const domain& at) 
	: AdT(of.group.ad( of.group.inv(std::get<1>(at) ) ) )  {
	
      }

      lie::coalgebra<domain> operator()(const lie::coalgebra<G>& p) const {
	return std::make_tuple(AdT(p), p);
      }
      
    };

    

  };

}


#endif
