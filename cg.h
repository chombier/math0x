#ifndef MATH0X_CG_H
#define MATH0X_CG_H

#include <math0x/iter.h>
#include <math0x/vector.h>

namespace math0x {


	template<int M = -1, class U = RR>
	struct cg {
		
		typedef U real;
		typedef vector<U, M> vec;
   
		math0x::iter iter;
		
		template<class Matrix>
		void solve(vec& x, const Matrix& A, const vec& b) const {
			assert( x.size() == b.size() );
			
			// residual / gradient
			vec r = b - A(x);

			// conjugate gradient
			vec p = r;
			
			vec Ap;
			
			real ro = r.squaredNorm();

			iter([&] {
					Ap = A(p);
					
					real alpha = ro / p.dot( Ap );
					
					x += alpha * p;
					r -= alpha * Ap;
					
					real ro_prev = ro;
					
					ro = r.squaredNorm();

					real beta = ro / ro_prev;

					p = r + beta * p;
					
					return std::sqrt(ro);
				});
		}

	};



}


#endif
