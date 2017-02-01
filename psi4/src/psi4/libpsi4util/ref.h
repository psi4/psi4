/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#ifndef _psi_src_lib_libmints_ref_h_
#define _psi_src_lib_libmints_ref_h_

/*!
    \file libmints/ref.h
    \ingroup MINTS
*/

#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <assert.h>

// Need to wrap the following into proprocessor wrappers
#ifdef HAVE_PTHREAD
    #include <pthread.h>
    #define __LOCK(l) pthread_mutex_lock(&l);
    #define __UNLOCK(l) pthread_mutex_unlock(&l);
#else
    #define __LOCK(l) 
    #define __UNLOCK(l) 
#endif

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
#ifdef HAVE_PTHREAD
    pthread_mutex_t _lock;
#else
#define _lock
#endif
public:
	SimpleReferenceCount() {
		_counter = NULL;
	}
	
public:
	template<typename T>
	void init(T*) {
		_counter = ::new size_t;
		*_counter = 1;
#ifdef HAVE_PTHREAD
		// Initialize pthread mutex
        pthread_mutex_init(&_lock, 0);
#endif
	}
	
	template<typename T> 
	void dispose(T*) {
        __LOCK(_lock);
		::delete _counter;
        _counter = 0;
        __UNLOCK(_lock);
#ifdef HAVE_PTHREAD
        pthread_mutex_destroy(&_lock);
#endif
	}
	
	template<typename T>
	void increment(T*) {
        __LOCK(_lock);
		++*_counter;
        __UNLOCK(_lock);
	}
	
	template<typename T> void decrement(T*) {
        __LOCK(_lock);
        if (_counter != 0)
		    --*_counter;
        __UNLOCK(_lock);
	}
	
	template<typename T> bool is_zero(T*) {
		return _counter?*_counter == 0:true;
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
    bool _managed;             // While false the counter is not incremented/decremented.
public:
    // default constructor
    Ref() : _managed(true) {
        _object_pointed_to = NULL;
    }

    // Test of copy constructor
    Ref(T* p) : _managed(true) {
        _object_pointed_to = NULL;
        if (p) {
            init(p);
        }
    }

    // explicit Ref(T* p) {
    //     init(p);
    // }

    Ref(Ref<T, CP, OP> const& cp) : CP(static_cast<CP const&>(cp)), OP(static_cast<OP const&>(cp)), _managed(cp._managed) {
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
            CP::operator=(static_cast<CP const&>(cp));
            OP::operator=(static_cast<OP const&>(cp));
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

    void set_managed(bool managed) { _managed = managed; }
    bool get_managed() const { return _managed; }
    
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
        _managed = cp._managed;
        if (cp._object_pointed_to != NULL && _managed) {
            CounterPolicy::increment(cp._object_pointed_to);
        }
    }

    void detach() {
        if (_object_pointed_to != NULL && _managed) {
            CounterPolicy::decrement(_object_pointed_to);
            if (CounterPolicy::is_zero(_object_pointed_to)) {
                CounterPolicy::dispose(_object_pointed_to);
                ObjectPolicy::dispose(_object_pointed_to);
            }
        }
    }
};

}

#endif /*Ref_H_*/