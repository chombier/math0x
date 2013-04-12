#ifndef MATH0X_CR_H
#define MATH0X_CR_H

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

			vec Ar = A(r);

			real ro = r.dot(s);

			// conjugate gradient
			vec p = r;
			vec Ap = s;
			
			iter([&] {
					
					real alpha = ro / Ap.squaredNorm();
					
					x += alpha * p;
					r -= alpha * Ap;
					
					Ar = A(r);
					
					real ro_prev = ro;
					
					ro = r.dot(s);

					real beta = ro / ro_prev;
					
					p = r + beta * p;
					Ap = s + beta * Ap;
					
					return std::sqrt(r.squaredNorm());
				});
		}

	};

}

#endif

