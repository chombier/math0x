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

namespace lie {

  template<class ... Args>
  struct traits< std::tuple<Args...> > {
    
    typedef std::tuple<Args...> G;

    // iterator type
    typedef tuple::range<Args...> each;
    
    typedef std::tuple< group<Args>... > args_type;
    args_type args;
  
    typedef std::tuple< lie::alg<Args>... > algebra;

    template<int I>
    using type = ::tuple::element<I, G>;


    struct get_alg {
      const args_type& args;

      template<int I>
      euclid::space< lie::alg< type<I> > > operator()() const {
	return std::get<I>(args).alg();
      }

    };

    euclid::space< algebra > alg() const { 
      return { each::map( get_alg{args} ) };
    }

    struct get_id {
      const args_type& args;
      
      template<int I>
      type<I> operator()() const {
	return std::get<I>(args).id();
      }
      
    };

    struct get_inv {
      const args_type& args;
      const G& x;
      
      template<int I>
      type<I> operator()() const {
	return std::get<I>(args).inv( std::get<I>(x) );
      }
      
    };

    struct get_prod {
      const args_type& args;

      const G& lhs;
      const G& rhs;
      
      template<int I>
      type<I> operator()() const {
	return std::get<I>(args).prod( std::get<I>(lhs),
				       std::get<I>(rhs));
      }
      
    };


    G id() const {
      return each::map( get_id{args} );
    }

    G inv(const G& x) const {
      return each::map( get_inv{args, x} );
    }

    G prod(const G& lhs, const G& rhs) const {
      return each::map( get_prod{args, lhs, rhs} );
    }

    struct get_group {
      const G& g;

      template<int I>
      lie::group< type<I> > operator()() const {
	return { std::get<I>(g) };
      }
      
    };

    traits( const G& g ) : args( each::map( get_group{g} ) ) { }
    traits( const args_type& args = args_type() ) : args(args) { }


    struct ad : func::tuple< lie::ad<Args>... > {
      
      struct get {
	const G& at;
	
	template<int I>
	lie::ad< type<I> > operator()() const {
	  return { std::get<I>(at) };
	}
      };
      
      ad( const G& at) : ad::self{ each::map( get{at} ) } { }
      
    };

    struct adT : func::tuple< lie::adT<Args>... > {
      
      struct get {
    	const G& at;
	
    	template<int I>
    	lie::adT< type<I> > operator()() const {
    	  return { std::get<I>(at) };
    	}
      };
      
      adT( const G& at) : adT::self{ each::map( get{at} ) } { }
      
    };
    
    struct exp : func::tuple< lie::exp<Args>... >  {
      struct get {
    	const group<G>& g;

    	template<int I>
    	lie::exp< type<I> > operator()() const {
    	  return std::get<I>(g.impl.args).exp();
    	}
      };

      exp( const group<G>& g) : exp::self{ each::map( get{g} ) }  { }
      
      
    };

    struct log : func::tuple< lie::log<Args>... >  {
      struct get {
    	const group<G>& g;

    	template<int I>
    	lie::log< type<I> > operator()() const {
    	  return std::get<I>(g.impl.args).log();
    	}
      };
      
      log( const group<G>& g) : log::self{ each::map( get{g} ) }  { }
     
      
      
    };


  };
}


#endif
