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
    
    static void sym_ortho(const real& a, const real& b,
			  real& c, real& s, real& r) {

    }

  public:

    struct lanczos {
      real& alpha;
      real& beta;
      vector<real, M>& v;

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

      vec r;			// residual
      real phi;			// residual norm
      real tau;
      
      real delta_1;
      real gamma_min; 		// minimum non-zero eigenvalue

      real norm; 		// matrix norm
      real cond;		// condition number

      real c, s;

      real eps;

      natural k;		// iteration 
    
      // r = b - Ax
      void init(const vec& rr, real threshold = 0) {
	natural n = rr.rows();
      
	auto zero = vec::Zero(n);
	
	r = rr;
	assert(!nan(r));
      
	beta = r.norm();
	assert(!nan(beta));
      
	if( beta <= threshold ) {
	  // we are done already lol
	  phi = 0;
	  return;
	}
      
	v_prev = zero;
	v = r.normalized();

	assert( !nan(v) );
 
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


    iter it;
    U sigma;
    
    minres(RR sigma = 0) : sigma(sigma) { }

    // solves (A - sigma * I) x = b, for symmetric A
    template<class Matrix>
    iter solve(vec& x, const Matrix& A, const vec& b) const {
      const NN n = b.size();
      
      if( x.empty() ) x = vec::Zero(n);
      
      data_type data;
      data.init(b, it.epsilon);
      
      return it( [&] {
	  data.step(x, A, sigma);
	  return data.phi;
	});
    };
  };

}


#endif
