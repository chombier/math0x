#ifndef MATH0X_MINRES_H
#define MATH0X_MINRES_H

#include <cmath>

#include <math0x/vector.h>
#include <math0x/iter.h>

namespace math0x {

	// minimum residual for symmetric systems (see [Choi06])
	template<int M = -1, class U = RR>
	class minres {
    
		typedef U real;
		typedef vector<U, M> vec;
    
		static inline real abs(const real& x) { return x < 0 ? -x : x; }
		static inline real sign(const real& x) { return x < 0 ? -1.0 : 1.0; }
    
		static void sym_ortho(const real& a, const real& b,
		                      real& c, real& s, real& r) {
			// TODO asserts ! nans !
			if( !b  ) {
				s = 0;
				r = abs( a );
      
				c = a ? sign(a) : 1.0;
      
			} else if ( !a ) {
				c = 0;
				s = sign(b);
				r = abs(b);
			} else {

				const real aabs = abs(a);
				const real babs = abs(b);
      
				if( babs > aabs ) {
					const real tau = a / b;
	
					s = sign( b ) / std::sqrt( 1 + tau * tau );
					c = s * tau;
					r = b / s;
				} else if( aabs > babs ) { 
					const real tau = b / a;
	    
					c = sign( a ) / std::sqrt( 1 + tau * tau );
					s = c * tau;
					r = a / c;
				}  else {
	  
					const real tau = b / a;
	  
					c = sign( a ) / std::sqrt( 1 + tau * tau );
					s = c * tau;
					r = a / c;
				}
	
			}
		}

	public:

		struct lanczos {
			real& alpha;
			real& beta;
			vec& v;

			template<class Matrix>
			static void step(const Matrix& A, 
			                 const vec& v,
			                 const vec& v_prev,
			                 real beta,
			                 real sigma,
			                 lanczos res)  {
				// we use res.v as work vector
				vec& p = res.v;
	
				p = A(v) - sigma * v;
	
				res.alpha = v.dot( p );
      
				p -= res.alpha * v + beta * v_prev;
      
				res.beta = res.v.norm();
      
				if( res.beta ) res.v /= res.beta;
			}

		};



		struct data_type {
			real beta;

			vec v_prev, v, v_next;
      
			vec d_next, d, d_prev;

			vec r;										// residual
			real phi;									// residual norm
			real tau;
      
			real delta_1;
			real gamma_min;						// minimum non-zero eigenvalue

			real norm;								// matrix norm
			real cond;								// condition number

			real c, s;

			real eps;

			natural k;								// iteration 
    
			// r = b - Ax
			void init(const vec& rr, real threshold = 0) {
				natural n = rr.rows();
      
				auto zero = vec::Zero(n);
	
				r = rr;
	
				beta = r.norm();

				if( beta <= threshold ) {
					// we are done already lol
					phi = 0;
					return;
				}
      
				v_prev = zero;
				v = r.normalized();

				// TODO nans !
				// assert( !nan(v) ); 
 
				v_next = zero;

				d = zero;
				d_prev = zero;

				phi = beta;
				tau = beta;
	
				delta_1 = 0;
				gamma_min = 0;

				norm = 0;
				cond = 0;

				c = -1;
				s = 0;

				eps = 0;
				k = 1;
			}

			template<class Matrix>
			void step(vec& x, const Matrix& A, real sigma = 0) {
	
				if( !phi ) return;
	
				real alpha;
				real beta_prev = beta;
				lanczos res{alpha, beta, v_next};
				lanczos::step( A, v, v_prev, beta, sigma, res );
	
				real delta_2 = c * delta_1  +  s * alpha;
				real gamma_1 = s * delta_1  -  c * alpha;

				real eps_next = s * beta;
				real delta_1_next = -c * beta;
	  
				real gamma_2;
				sym_ortho(gamma_1, beta, c, s, gamma_2);
	
				tau = c * phi;
				phi = s * phi;
	
				r = (s * s) * r - (phi * c) * v_next;
  
				norm = (k == 1) ? std::sqrt( alpha * alpha  +  beta * beta ) 
					: std::max(norm, std::sqrt( alpha * alpha  +  beta * beta  +  beta_prev * beta_prev));
	  
				if( gamma_2 ) {
	  
					d_next = (v  -  delta_2 * d  -  eps * d_prev ) / gamma_2;
					x += tau * d_next;
	  
					gamma_min = (k == 1) ? gamma_2 
						: std::min(gamma_min, gamma_2);
	  
					assert( gamma_min );
	  
					cond = norm / gamma_min;
	  
					// swap pointers:
					// d <- d_next 
					d_prev.swap(d);
					d.swap(d_next);
	  
					// v <- v_next
					v_prev.swap(v);
					v.swap(v_next);
	  
					eps = eps_next;
					delta_1 = delta_1_next;
	  
				} else {
					// core::log( DERP );
				}
	   
				++k;
			}


		};


		math0x::iter iter;
		U sigma;
    
		minres(RR sigma = 0) : sigma(sigma) { }

		// solves (A - sigma * I) x = b, for symmetric A
		template<class Matrix>
		math0x::iter solve(vec& x, const Matrix& A, const vec& b) const {
			const NN n = b.size();
      
			if( !x.size() ) x = vec::Zero(n);
      
			data_type data;
			data.init(b, iter.epsilon);
      
			return iter( [&] {
					data.step(x, A, sigma);
					return data.phi;
				});
		};
	};

}


#endif
