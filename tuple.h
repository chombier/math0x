#ifndef GROUP_TUPLE_H
#define GROUP_TUPLE_H

#include <group/euclid.h>

#include <group/tuple/range.h>
#include <group/tuple/element.h>
#include <group/tuple/find.h>

#include <cassert>

namespace euclid {
  
  template<class ... Args>
  struct traits< std::tuple<Args...> > {
    
    // iterator type
    typedef tuple::range<Args...> each;
    
    typedef std::tuple< space<Args>... > impl_type;
    impl_type impl;
    
    typedef tuple::repeat<NN, sizeof...(Args) > offset_type;
    offset_type offset;
    
    NN dimension;

    typedef std::tuple<Args...> E;

    // choose first as field TODO check consistency
    typedef euclid::field< tuple::element<0, E> > field;
    typedef std::tuple< euclid::dual<Args>... > dual;

    // some helper functors for implementing operations
    struct helper {

      // term-wise euclidean structure
      struct impl {
	
	const E& value;

	template<int I>
	euclid::space< tuple::element<I, E> > operator()() const {
	  return { std::get<I>(value) };
	}
	
      };

      // term-wise dual-space
      struct dual {
  	const impl_type& impl;

  	template<int I>
  	auto operator()() const -> decltype( *std::get<I>(impl) ) {
  	  return *std::get<I>(impl);
  	}
      };

      // computes successive offsets
      struct offset {
  	NN& last;

  	offset_type& offset;
  	const impl_type& impl;
	
  	template<int I>
  	void operator()() const {
  	  std::get<I>(offset) = last;
  	  last += std::get<I>(impl).dim();
  	}

      };

      // coordinate access
      struct coord {
  	const E& value;
  	const impl_type& impl;
  	const offset_type& offset;
	
  	const field*& res;
  	NN i;
	
  	template<int I>
  	void operator()() const {
  	  const auto& impl_i = std::get<I>(impl);
  	  NN index = i - std::get<I>(offset);
  	  res = & impl_i.coord( index, std::get<I>(value) );
  	}
	
      };

      // term-wise zero element
      struct zero {
	
	const impl_type& impl;
	
	template<int I>
	auto operator()() const -> decltype( std::get<I>(impl).zero() ) {
	  return std::get<I>(impl).zero();
	}
	
      };
      
    };
    
    // builds offset tuple and returns total dimension
    NN make_offset() {
      NN last = 0;
      each::apply( typename helper::offset{last, offset, impl} );
      return last;
    }
    
    // dual geometry
    space< dual > operator*() const {
      return { each::map( typename helper::dual{impl} ) };
    }
    
    // coordinate access
    const field& coord(NN i, const E& x) const {
      assert( i < dim() );

      const field* res = 0;
      tuple::find(offset, i, typename helper::coord{x, impl, offset, res, i} );
      
      return *res;
    }

    field& coord(NN i, E& x) const {
      assert( i < dim() );

      const field* res = 0;
      tuple::find(offset, i, typename helper::coord{x, impl, offset, res, i} );

      // uglyness i has it      
      return *const_cast<field*>(res);
    }
    
    
    E zero() const {
      return each::map( typename helper::zero{impl} );
    }

    NN dim() const { return dimension; }

    // ctors
    traits( const impl_type& impl = impl_type() ) 
      : impl(impl) {
      dimension = make_offset();
    }

    traits( const space<Args>&... args ) 
      : impl(args...) {
      dimension = make_offset();
    }

    traits( const E& value )
    : impl( each::map( typename helper::impl{value} ) ) {
      dimension = make_offset();
    }
    
  };

}


// TODO finish lie traits

namespace lie {

  template<class ... Args>
  struct traits< std::tuple<Args...> > {
    
    // iterator type
    typedef tuple::range<Args...> each;
    
    typedef std::tuple< group<Args>... > impl_type;
    impl_type impl;
  
    typedef std::tuple< lie::algebra<Args>... > algebra;
    
    
    struct adjoint {
      // TODO 
    };

    struct coadjoint {
      // TODO 
    };

    struct exponential  {
      // TODO 
    };

    struct logarithm  {
      // TODO 
    };


  };
}


#endif
