#ifndef _psi_src_lib_libutil_ref_h_
#define _psi_src_lib_libutil_ref_h_

namespace psi {
  
  // For most object counted by Refs we can use the following simple object policy:
  class StandardObjectPolicy {
  public:
  	template<typename T> void dispose(T* object) {
  		delete object;
  	}
  };

  // Clearly, the above policy will not work for arrays allocated with operator new[].
  // A replacement policy for this case is trival, fortunately:
  class StandardArrayPolicy {
  public:
  	template<typename T> void dispose(T* array) {
  		delete[] array;
  	}
  };

  class SimpleReferenceCount {
  private:
  	size_t* _counter;
  public:
  	SimpleReferenceCount() {
  		_counter = NULL;
  	}
	
  public:
  	template<typename T> void init(T*) {
  		_counter = ::new size_t;
  		*_counter = 1;
  	}
	
  	template<typename T> void dispose(T*) {
  		::delete _counter;
  	}
	
  	template<typename T> void increment(T*) {
  		++*_counter;
  	}
	
  	template<typename T> void decrement(T*) {
  		--*_counter;
  	}
	
  	template<typename T> bool is_zero(T*) {
  		return *_counter == 0;
  	}
  };
  // Reference counting pointer:
  template<typename T, typename CounterPolicy = SimpleReferenceCount, typename ObjectPolicy = StandardObjectPolicy>
  class Ref : private CounterPolicy, private ObjectPolicy {
    protected:
    // shortcuts
    typedef CounterPolicy CP;
    typedef ObjectPolicy  OP;

    T* _object_pointed_to;		// Object referred to (or NULL if none)

    public:
    // default constructor
    Ref() {
      _object_pointed_to = NULL;
    }

    // Test of copy constructor
    Ref(T* p) {
      _object_pointed_to = NULL;
      if (p) {
        init(p);
      }
    }

    // explicit Ref(T* p) {
    //     init(p);
    // }

    Ref(Ref<T, CP, OP> const& cp) : CP((CP const&) cp), OP((OP const&)cp) {
      _object_pointed_to = NULL;
      attach(cp);
    }

    ~Ref() {
      detach();
    }

    Ref<T,CP,OP>& operator= (T* p) {
      assert(p != _object_pointed_to);
      detach();
      init(p);
      return *this;
    }

    Ref<T,CP,OP>& operator= (Ref<T,CP,OP> const& cp) {
      if (_object_pointed_to != cp._object_pointed_to) {
        detach();
        CP::operator=((CP const&)cp);
        OP::operator=((OP const&)cp);
        attach(cp);
      }
      return *this;
    }

    T* operator-> () const {
      return _object_pointed_to;
    }

    T& operator* () const {
      return *_object_pointed_to;
    }

    operator T*() const { return _object_pointed_to; }

    T* pointer() const {
    #ifdef DEBUG
      assert(_object_pointed_to != NULL);
    #endif
      return _object_pointed_to;
    }

    T& operator[] (int i) const {
      return _object_pointed_to[i];
    }
    private:
    void init(T* p) {
      if (p != NULL) {
        CounterPolicy::init(p);
      }
      _object_pointed_to = p;
    }

    void attach(Ref<T,CP,OP> const& cp) {
      _object_pointed_to = cp._object_pointed_to;
      if (cp._object_pointed_to != NULL) {
        CounterPolicy::increment(cp._object_pointed_to);
      }
    }

    void detach() {
      if (_object_pointed_to != NULL) {
        CounterPolicy::decrement(_object_pointed_to);
        if (CounterPolicy::is_zero(_object_pointed_to)) {
          CounterPolicy::dispose(_object_pointed_to);
          ObjectPolicy::dispose(_object_pointed_to);
        }
      }
    }
  };

}

#endif
