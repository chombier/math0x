#ifndef MATH0X_SO3_H
#define MATH0X_SO3_H

#include <math0x/types.h>
#include <math0x/lie.h>

#include <math0x/quaternion.h>

#include <math0x/vector.h>
#include <math0x/covector.h>
#include <math0x/matrix.h>

#include <math0x/error.h>

#include <boost/math/special_functions/sinc.hpp>

namespace math0x { 

	// 3-dimensional rotation group
	template<class U>
	class SO<3, U> {
		typedef math0x::quaternion<U> quaternion_type;
		quaternion_type quaternion;
	public:
		typedef SO base;

		const quaternion_type& quat() const { return quaternion; }

		SO( const quaternion_type& quaternion = quaternion_type::Identity() ) 
			: quaternion( quaternion ) {
			// TODO precision
			assert( std::abs( quat().norm() -  1 ) < 1e-14 );
		}
		
		SO operator*(const SO& other) const {
			// TODO normalize ?
			return {quat() * other.quat()};
		}
  
		SO inv() const { 
			return { quat().conjugate() };
		}
  
		typedef vector<U, 3> vec_type;

		vec_type operator()(const vec_type& x) const { 
			return quat() * x;
		}

		typedef math0x::matrix<U, 3, 3> matrix_type;
		matrix_type mat() const {
			return quat().toRotationMatrix();
		}
		
		struct push : SO {
      
			push(const SO& of, const vec_type& ) 
				: push::base(of) {

			}

		};

		struct pull {
      
			SO of_inv;
      
			pull(const SO& of, const vec_type& ) 
				: of_inv(of.inv()) {
	
			}
      
			lie::coalgebra<vec_type> operator()(const lie::coalgebra<vec_type>& f) const {
				return of_inv(f.transpose()).transpose();
			}
      
		};
  

	};


	namespace lie {

		template<class U>
		struct traits< SO<3, U> > {

			// TODO make so<3, U> and wrap skew-symmetric matrices
			typedef vector<U, 3> algebra;
    
			typedef SO<3, U> G;

			traits() { }
			traits(const G& ) { }

			G id() const { return {}; }
			G inv(const G& x) const { return x.inv(); }
			G prod(const G& x, const G& y) const { return x * y; }
    
			euclid::space<algebra> alg() const { return {}; }
    
      
			algebra bracket(const algebra& x, const algebra& y) const {
				return x.cross(y);
			}


			struct Ad :  G {
				Ad(const G& at) : G(at) { }
			};
    
			// TODO better ?
			struct AdT { 
				typedef AdT base;
      
				G at;
      
				AdT( const G& g ) : at( g.inv() ) { }
      
				lie::coalgebra<G> operator()(const lie::coalgebra<G>& f) const {
					return at(f.transpose()).transpose();
				}
      
			};


			struct exp {
				typedef exp base;
      
				exp(const group<G>& = group<G>() ) { }

				G operator()(const lie::algebra<G>& w) const { 
					const U theta = 0.5 * w.norm();
					const U sinc = 0.5 * boost::math::sinc_pi(theta);
					
					return quaternion<U>( std::cos(theta), 
					                      sinc * w.x(),
					                      sinc * w.y(),
					                      sinc * w.z() );
				}
				 
				// TODO push/pull
				
			};


			struct log {
				typedef log base;
      
				log(const group<G>& = group<G>() ) { }
      
				lie::algebra<G> operator()(const G& g) const { 
					quaternion<U> q = g.quat();
					if( q.w() < 0 ) q.coeffs() = - q.coeffs();
					
					// sanity clamp
					const U w = std::min(1.0, q.w());
					assert( w >= 0 );

					const U theta = std::acos( w );
					const U sinc = 0.5 * boost::math::sinc_pi(theta);
					
					return q.vec() / sinc;
				}
      
				// struct push { push(const log&, const G&) { } };
				// struct pull { pull( const log&, const G& ) { } }; 
      
			};
    
    
		};


	}


	namespace func {

		template<class U>
		struct apply< SO<3, U> > {
      
			typedef SO<3, U> G;

			typedef std::tuple<G, func::domain<G> > domain;
			typedef func::range<G> range;
      
			range operator()(const domain& x) const {
				return std::get<0>(x)(std::get<1>(x));
			}
      

			struct push {

				domain at;
				lie::group<G> group;
	
				push( const apply&, const domain& at) : at(at) { }

				lie::algebra<range> operator()(const lie::algebra<domain>& dx) const {
	  
					return std::get<0>(at)( group.bracket(std::get<0>(dx), std::get<1>(at))
					                        + std::get<1>(dx) );
				}
	
			};

			// TODO pull
      

		};

	}


}



#endif
