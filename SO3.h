#ifndef MATH0X_SO3_H
#define MATH0X_SO3_H

#include <math0x/types.h>
#include <math0x/lie.h>

#include <math0x/quaternion.h>

#include <math0x/vector.h>
#include <math0x/covector.h>
#include <math0x/matrix.h>
#include <math0x/so3.h>

#include <math0x/error.h>
#include <math0x/epsilon.h>

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
			
			assert( std::abs( quat().norm() -  1 ) < epsilon<U>() );
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

			push(const SO& of) : push::base(of) { }
			
		};

		struct pull {
      
			SO of_inv;
      
			pull(const SO& of, const vec_type& ) 
				: of_inv(of.inv()) {
	
			}

			pull(const SO& of) : of_inv(of.inv()) { }
      
			lie::coalgebra<vec_type> operator()(const lie::coalgebra<vec_type>& f) const {
				return of_inv(f.transpose()).transpose();
			}
			
		};
  

	};


	namespace lie {

		template<class U>
		struct traits< SO<3, U> > {

			// TODO make so<3, U> and wrap skew-symmetric matrices ?
			typedef vector<U, 3> algebra;
			typedef covector<U, 3> coalgebra;
    
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
    
			class AdT { 
				func::pull<G> atT;
			public:
				
				AdT( const G& g ) : atT( g ) { }
				
				lie::coalgebra<G> operator()(const lie::coalgebra<G>& f) const {
					return atT(f);
				}
				
			};


			struct ad {
				
				algebra x;
				
				ad(const algebra& x, const group<G>& = {} ) : x(x) { }

				algebra operator()(const algebra& y) const {
					return x.cross(y);
				}
				
			};


			struct exp {
				typedef exp base;
      
				exp(const group<G>& = {} ) { }

				G operator()(const lie::algebra<G>& w) const { 
					U theta = 0.5 * w.norm();
					U sinc = 0.5 * boost::math::sinc_pi(theta);
					
					return quaternion<U>( std::cos(theta), 
					                      sinc * w.x(),
					                      sinc * w.y(),
					                      sinc * w.z() );
				}
				 
				
				class push {
					algebra x;
					U theta2;
					U theta;

					bool origin;
					
					G qT;
				public:
					const G& value_inv() const { return qT; }
					

					push(const exp& of, const algebra& at) 
						: x(at),
						  theta2(x.squaredNorm()),
						  theta( std::sqrt(theta2) ),
						  origin( theta < epsilon<U>() ),
						  qT( origin ? G() : of(-at) )
					{ }
					
					algebra operator()(const algebra& h) const {
						
						// easy peasy
						if( origin ) return h;
      
						// we decompose h = u + ad(x).v, with u in ker(ad(x))
						algebra u = x * ( x.dot(h) / theta2 );

						// w = ad(x).v
						algebra w = h - u;
						
						// obtain v
						algebra v = w.cross(x) / theta2;
						
						// general result is u + (I - Ad(exp(-x))).v
						return u + v - qT(v);
					}
	
				};
      

				class pull {
					algebra x;
					
					U theta2;
					U theta;
					bool origin;
					
					G q;
				public:
					
					pull(const exp& of, const algebra& at) 
						: x( at),
						  theta2( x.squaredNorm() ),
						  theta( std::sqrt( theta2 ) ),
						  origin( theta < epsilon<U>() ),
						  q( origin ? G() : of(at) )
					{ }
					
					coalgebra operator()(const coalgebra& f) const {
						// easy peasy
						if( origin ) return f;
						
						// convenience
						algebra fT = f.transpose();
						
						// f pulled on u 
						algebra fTu = fT;

						// f pulled on v
						algebra fTv = fT - q(fT);
						
						// f pulled on w
						algebra fTw = x.cross(fTv) / theta2;
						
						// f pulled on dx
						algebra fTdx = fTw + x * (x.dot( fTu - fTw ) / theta2); 
						
						// and back into shape
						return fTdx.transpose();
					}
					
				};

			};
		


			struct log {
				typedef log base;
      
				log(const group<G>& = {} ) { }
      
				lie::algebra<G> operator()(const G& g) const { 
					quaternion<U> q = g.quat();
					if( q.w() < 0 ) q.coeffs() = - q.coeffs();
					
					// sanity clamp
					U w = std::min(1.0, q.w());
					assert( w >= 0 );

					U theta = std::acos( w );
					U sinc = 0.5 * boost::math::sinc_pi(theta);
					
					return q.vec() / sinc;
				}
			
				class push {
					
					G at;
					algebra x;
					U theta2, theta;
					bool origin;
					U scale;
					
				public:
					const algebra& value() const { return x; }
					
					push( const log& of, const G& at) 
						: at(at),
						  x(of(at)),
						  theta2( x.squaredNorm() ),
						  theta( std::sqrt( theta2 ) ),
						  origin( theta < epsilon<U>() ),
						  scale( origin ? 0 : 1 / ( 2 * (1 - std::cos( theta ) ) ) )
					{
						
					}
					
					algebra operator()(const algebra& h) const {
						
						// easy peasy
						if( origin ) return h;
						
						// we decompose h = u + (I - R^T).v, with u in ker(I - R^T) ~ x
						const algebra u = x * ( x.dot(h) / theta2 );
						
						// w = (I - R^T).v
						const algebra w = h - u;
						
						// obtain v: invert complex number 1 - e^{-i theta}
						const algebra v = scale * (w - at(w));
						
						return u + x.cross(v);
					}


				};


				class pull {
					
					G atT;
					algebra x;
					U theta2, theta;
					bool origin;
					U scale;
					
				public:
					pull( const log& of, const G& at) 
						: atT(at.inv()),
						  x(of(at)),
						  theta2( x.squaredNorm() ),
						  theta( std::sqrt( theta2 ) ),
						  origin( theta < epsilon<U>() ),
						  scale( origin ? 0 : 1 / ( 2 * (1 - std::cos( theta ) ) ) )
					{
						
					}
					
					coalgebra operator()(const coalgebra& f) const {
						
						// easy peasy
						if( origin ) return f;
      
						// convenience
						algebra fT = f.transpose();
						
						// f pulled on u 
						algebra fTu = fT;
						
						// f pulled on v
						algebra fTv = -x.cross(fT);
						
						// f pulled on w
						algebra fTw = scale * (fTv - atT(fTv));

						// f pulled on dg
						algebra fTdg = fTw + x * (x.dot( fTu - fTw ) / theta2); 
						
						return fTdg.transpose();

					}


				};

			
			};
    
    
		};


	}


	namespace func {

		template<class U>
		struct apply< SO<3, U> > {
			
			typedef SO<3, U> G;
			
			typedef std::tuple<G, func::domain<G> > domain;
			typedef func::range<G> range;
			
			apply( const lie::group<G>& = {} ) { }

			range operator()(const domain& x) const {
				return std::get<0>(x)(std::get<1>(x));
			}
      

			struct push {

				domain at;
				
				push( const apply&, const domain& at) : at(at) { }
				
				lie::algebra<range> operator()(const lie::algebra<domain>& v) const {
					return std::get<0>(at)( std::get<0>(v).cross(std::get<1>(at)) + std::get<1>(v) );
				}
				
			};
			
			struct pull {
				
				domain at;
				
				pull( const apply&, const domain& at) : at(at) { }
				
				lie::coalgebra<domain> operator()(const lie::coalgebra<range>& p) const {
					lie::algebra<range> fT = p.transpose();

					lie::algebra<range> RTfT = std::get<0>(at).inv()(fT);
					
					return lie::coalgebra<domain>{std::get<1>(at).cross(RTfT).transpose(), RTfT.transpose()};
				}
	
			};


		};

	}


}



#endif
