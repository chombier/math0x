#ifndef MATH0X_FUNC_PROD_H
#define MATH0X_FUNC_PROD_H

#include <math0x/lie.h>

namespace math0x { 
	namespace func {

		// lie group product
		template<class G>
		struct prod {
			typedef prod base;
    
			lie::group<G> group;
    
			prod(const lie::group<G>& group = {} ) : group(group) { }
    
			typedef std::tuple<G, G> domain;
    
			G operator()(const domain& g) const {
				return group.prod( std::get<0>(g), 
				                   std::get<1>(g) );
			}
    
    
			struct push {

				euclid::space< lie::algebra<G> > alg;
				lie::Ad<G> Ad;

				push(const prod& of, const domain& at)
					: alg(of.group.alg()), 
					  Ad(of.group.Ad( of.group.inv(std::get<1>(at) ) ) )
				{
	
				}

				lie::algebra<G> operator()(const lie::algebra<domain>& v) const {
					return alg.sum( Ad(std::get<0>(v)), 
					                std::get<1>(v) );
				}
      
			};


			struct pull {

				lie::AdT<G> AdT;

				pull(const prod& of, const domain& at) 
					: AdT(of.group.AdT( of.group.inv(std::get<1>(at) ) ) )  {
	
				}

				lie::coalgebra<domain> operator()(const lie::coalgebra<G>& p) const {
					return std::make_tuple(AdT(p), p);
				}
      
			};
    

		};


		// convenience 
		template<class G>
		prod<G> make_prod(const lie::group<G>& group) { return {group}; }
		

	}

}
#endif
