#ifndef MATH0X_ALGO_BLUR_H
#define MATH0X_ALGO_BLUR_H

#include <math0x/vector.h>
#include <math0x/lie.h>
#include <math0x/each.h>

namespace math0x {
	namespace algo {
		
		// a time-domain bluring filter for general Lie group data,
		// generalized from "General Construction of Time-Domain Filters
		// for Orientation Data" by Lee & Shin, 2001
		
		template<class G>
		vector<G> blur(const vector<G>& data) {
			NN n = data.size();
			assert( n > 1 );
			
				vector<G> res; res.resize( n );

			auto group = lie::group_of(data(0));
			auto alg = group.alg();
			
			auto log = group.log();
			auto exp = group.exp();
			
			// linearized deltas
			vector< lie::algebra<G> > w;
			w.resize( n - 1 );
			
			each(w, [&](NN i) {
					w(i) = log(group.prod(group.inv(data(i)), data(i + 1) ) );
				});
			
			
			const lie::algebra<G> zero = alg.zero();
			
			each(data, [&](ZZ i) {
					
					lie::algebra<G> w_im2 = i >= 2 ? w(i - 2) : zero;
					lie::algebra<G> w_im1 = i >= 1 ? w(i - 1) : zero;
					lie::algebra<G> w_i =  w(i);
					lie::algebra<G> w_ip1 = i < int(n - 1) ? w(i + 1) : zero;
					
					// lisp style, i has it
					lie::algebra<G> sum = alg.scal
						( 1.0 / 16,
						  alg.sum
						  ( alg.scal(-1, w_im2),
						    alg.sum
						    ( alg.scal(-5, w_im1),
						      alg.sum
						      ( alg.scal(5, w_i),
						        alg.scal( 1, w_ip1) 
						        ) 
						      ) 
						    ) 
						  );
					
					res(i) = group.prod(data(i), exp(sum));
					
				});
			
			return res;
		}

	}
}


#endif
