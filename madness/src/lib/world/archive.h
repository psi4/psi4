/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680


  $Id: archive.h 2173 2011-02-23 21:40:46Z justus.c79@gmail.com $
*/

#ifndef MADNESS_WORLD_ARCHIVE_H__INCLUDED
#define MADNESS_WORLD_ARCHIVE_H__INCLUDED

/*!
  \file world/archive.h
  \brief Interface templates for the archives (serialization)
  \addtogroup serialization

  The programmer should not need to include world/archive.h directly.  Instead,
  include the header file for the actual archive (binary file, text/xml
  file, vector in memory, ...) that you want to use.

  \par Background

  The interface and implementation are deliberately modelled, albeit
  loosely, upon the Boost serialization class (thanks boost!). The major
  differences are that this archive class does \em not break cycles
  and does \em not automatically store unique copies of data
  referenced by multiple objects.  Also, classes are responsbible
  for managing their own version information.  At the lowest level,
  the interface to an archive also differs to facilitate
  vectorization and high-bandwidth data transfer.  The
  implementation employs templates that are almost entirely inlined.
  This should enable low-overhead use of archives in applications
  such as interprocess communication.

  \par How to use an archive?

  An archive is a uni-directional stream of typed data to/from disk,
  memory, or another process.  Whether the stream is for input or for
  output, you can use the \c & operator to transfer data to/from the
  stream.  If you really want, you can also use the \c << and \c >>
  for output or input, respectively, but there is no reason to do so.
  The \c & operator chains just like \c << for \c cout or \c >> for \c
  cin.  You may discover in \c archive.h other interfaces but you
  should \em not use them --- use the \& operator!  The lower level
  interfaces will probably not, or only inconsistently incorpoate
  type information and may even appear to work when they are not.

  Unless type checking has not been implemented by an archive for
  reasons of efficiency (e.g., message passing) a C-string exception
  will be thrown on a type-mismatch when deserializing.  End-of-file,
  out-of-memory and other others also generate string exceptions.

  Fundamental types (see below), STL complex, vector, strings, pairs
  and maps, and tensors (int, long, float, double,
  float_complex, double_complex) all just work without you doing
  anything, as do fixed dimension arrays of the same (STL allocators
  are not presently accomodated). E.g.,
  \code
  bool finished=false;
  int info[3] = {1,33,2};
  map<int,double> fred;
  map[0]=55.0; map[1]=99.0;
  BinaryFstreamOutputArchive ar('restart.dat');
  ar & map & info & finished;
  \endcode
  Deserializing is identical, except that you need to use an input archive,
  c.f.,
  \code
  bool finished;
  int info[3];
  map<int,double> fred;
  BinaryFstreamInputArchive ar('restart.dat');
  ar & map & info & finished;
  \endcode

  Variable dimension and dynamically allocated arrays do not have
  their dimension encoded in their type.  The best way to
  (de)serialize them is to wrap them in an \c archive_array as follows.
  \code
  int a[n]; // n is not known at compile time
  double *p = new double[n];
  ar & wrap(a,n) & wrap(p,n);
  \endcode
  The \c wrap() function template is a factory function to simplify
  instantiation of a correctly typed \c archive_array template.
  Note that when deserializing you must have first allocated the
  array --- the above code can be used for both serializing and
  deserializing.  If you want the memory to be automatically allocated
  consider using either an STL vector or a madness tensor.

  To transfer the actual value of a pointer to a stream (is this really
  what you want?) then store an archive_ptr wrapping it.  The factor
  function wrap_ptr() assists in doing this, e.g., here
  for a function pointer
  \code
  int foo();
  ar & wrap_ptr(foo);
  \endcode

  \par User-defined types

  User-defined types require a little more effort.  Three
  cases are distinguished.
  - symmetric load and store
     - intrusive
     - non-intrusive
  - non-symmetric load and store

  We will examine each in turn, but we first need to discuss a little
  about the implementation.

  When transfering an object \c obj to/from an archive \c ar with \c ar&obj,
  you are invoking the templated function
  \code
  template <class Archive, class T>
  inline const Archive& operator&(const Archive& ar, T& obj);
  \endcode
  that then invokes other templated functions to redirect to input or
  output streams as appropriate, manage type checking, etc..  We would
  now like to overload the behaviour of these functions in
  order accomodate your fancy object.  However, function templates
  cannot be partially specialized.  Following the technique
  recommended <a href=http://www.gotw.ca/publications/mill17.htm>here</a>
  (look for moral#2), each of the templated functions directly calls
  a member of a templated class.  Classes, unlike functions, can be
  partially specialized so it is easy to control and predict what is
  happening.  Thus, in order to change the behaviour of all archives
  for an object you just have to provide a partial specialization of
  the appropriate class(es).  Do \em not overload any of the function
  templates.

  <em>Symmetric instrusive method</em>

  Many classes can use the same code for serializing and
  deserializing.  If such a class can be modified, the cleanest way
  of enabling serialization is to add a templated method as follows.
  \code
  class A {
      float a;
  public:
      A(float a = 0.0) : a(a) {}

      template <class Archive>
      inline void serialize(const Archive& ar) {ar & a;}
  };
  \endcode

  <em>Symmetric non-intrusive method</em>

  If a class with symmetric serialization cannot be modified, then
  you can define an external class template with the following signature in
  the \c madness::archive namespace (where \c Obj is the name of
  your type).
  \code
  namespace madness {
      namespace archive {
          template <class Archive>
          struct ArchiveSerializeImpl<Archive,Obj> {
              static inline void serialize(const Archive& ar, Obj& obj);
          };
      }
  }
  \endcode

  For example,
  \code
  class B {
  public:
      bool b;
      B(bool b = false) : b(b) {};
  };

  namespace madness {
      namespace archive {
	  template <class Archive>
	  struct ArchiveSerializeImpl<Archive,B> {
	      static inline void serialize(const Archive& ar, B& b) {ar & b.b;};
	  };
      }
  }
  \endcode

  <em>Non-symmetric non-intrusive</em>

  For classes that do not have symmetric (de)serialization you must
  define separate partial templates for the functions \c load and \c store
  with these signatures and again in the \c madness::archive namespace.
  \code
  namespace madness {
      namespace archive {
	  template <class Archive>
	  struct ArchiveLoadImpl<Archive,Obj> {
	      static inline void load(const Archive& ar, Obj& obj);
	  };

	  template <class Archive>
	  struct ArchiveStoreImpl<Archive,Obj> {
	      static inline void store(const Archive& ar, Obj& obj);
	  };
      }
  }
  \endcode

  First a simple, but artificial example.
  \code
  class C {
  public:
      long c;
      C(long c = 0) : c(c) {};
  };

  namespace madness {
      namespace archive {
          template <class Archive>
	  struct ArchiveLoadImpl<Archive,C> {
	      static inline void load(const Archive& ar, C& c) {ar & c.c;}
          };

	  template <class Archive>
	  struct ArchiveStoreImpl<Archive,C> {
	      static inline void store(const Archive& ar, const C& c) {ar & c.c;}
	  };
      }
  }
  \endcode

  Now a more complicated example that genuinely requires asymmetric
  load and store.  First, a class definition for a simple linked list.
  \code
  class linked_list {
      int value;
      linked_list *next;
  public:
      linked_list(int value = 0) : value(value), next(0) {};

      void append(int value) {
          if (next) next->append(value);
          else next = new linked_list(value);
      };

      void set_value(int val) {value = val;};

      int get_value() const {return value;};

      linked_list* get_next() const {return next;};
  };
  \endcode
  And this is how you (de)serialize it.
  \code
  namespace madness {
      namespace archive {
	  template <class Archive>
	  struct ArchiveStoreImpl<Archive,linked_list> {
	      static void store(const Archive& ar, const linked_list& c) {
		  ar & c.get_value() & bool(c.get_next());
		  if (c.get_next()) ar & *c.get_next();
	      };
	  };

	  template <class Archive>
	  struct ArchiveLoadImpl<Archive,linked_list> {
	      static void load(const Archive& ar, linked_list& c) {
		  int value;  bool flag;
		  ar & value & flag;
		  c.set_value(value);
		  if (flag) {
		      c.append(0);
		      ar & *c.get_next();
		  }
	      };
	  };
      }
  }
  \endcode

  Given the above implementation of a linked list, you can
  (de)serialize an entire list using a single statement.
  \code
  linked_list list(0);
  for (int i=1; i<=10; ++i) list.append(i);
  BinaryFstreamOutputArchive ar('list.dat');
  ar & list;
  \endcode

  \par Non-default constructor

  There are various options for objects that do not have a default
  constructor.  The most appealing and totally non-intrusive
  approach is to define load/store functions for a pointer to the
  object.  Then in the load method you can deserialize all of the
  information necessary to invoke the constructor and return a
  pointer to a new object.

  Things that you know are contiguously stored in memory
  and are painful to serialize with full type safety can be
  serialized by wrapping opaquely as byte streams using the
  \c wrap_opaque() interface.  However, this should be regarded
  as a last resort.

  \par Type checking and registering your own types

  To enable type checking for user-defined types you must register
  them with the system.  There are 64 empty slots for user types
  beginning at cookie=128.  Type checked archives (currently all
  except the MPI archive) store a cookie (byte with value 0-255) with
  each datum.  Unknown (user-defined) types all end up with the same
  cookie indicating unkown --- i.e., no type checking unless you
  register.

  Two steps are required to register your own types (e.g., here
  for the types \c Foo and \c Bar )
  -# In a header file after including world/archive.h, associate
  your types and pointers to them with cookie values
  \code
  namespace madness {
      namespace archive {
	  ARCHIVE_REGISTER_TYPE_AND_PTR(Foo,128);
	  ARCHIVE_REGISTER_TYPE_AND_PTR(Bar,129);
      }
  }
  \endcode
  -# In a single source file containing your initialization routine
  define a macro to force instantiation of relevant templates
  \code
  #define ARCHIVE_REGISTER_TYPE_INSTANTIATE_HERE
  \endcode
  and then in the initalization routine register the
  name of your types as follows
  \code
    ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(Foo);
    ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(Bar);
  \endcode
  Have a look at
  <a href=http://code.google.com/p/m-a-d-n-e-s-s/source/browse/local/trunk/src/lib/world/testar.cc>this test code</a>.
  to see things in action.

  \par Types of archive

  Presently provided are
  - world/textfsar.h --- (text \c std::fstream ) a file in text (XML)
  - world/binfsar.h --- (binary \c std::fstream ) a file in binary
  - world/vecar.h --- binary in memory using an \c std::vector<unsigned_char>
  - world/bufar.h --- binary in memory buffer (this is rather heavily specialized for internal
    use so applications should use a vector instead)
  - world/mpiar.h --- binary stream for point-to-point communication
    using MPI (non-typesafe for efficiency).
  - world/parar.h --- parallel archive to binary file with multiple readers/writers.
    Mostly to support efficient transfer of large WorldContainer (world/worlddc.h)
    and MADNESS Function (mra/mra.h) objects, though any serialiable object
    can employ it.

  The buffer and \c vector archives are bitwise identical to the
  binary file archive.

  \par Implementing a new archive

  Minimally, an archive must derive from either BaseInputArchive or
  BaseOutputArchive and define for arrays of fundamental types
  either a \c load or \c store method, as appropriate.  Additional
  methods can be provided to manipulate the target stream.  Here is
  a simple, but functional, implementation of a binary file archive.
  \code
  #include <fstream>
  #include <world/archive.h>
  using namespace std;
  class OutputArchive : public BaseOutputArchive {
    mutable ofstream os;
  public:
    OutputArchive(const char* filename)
      : os(filename, ios_base::binary | ios_base::out | ios_base::trunc) {};

    template <class T>
    void store(const T* t, long n) const {
      os.write((const char *) t, n*sizeof(T));
    }
  };

  class InputArchive : public BaseInputArchive {
    mutable ifstream is;
  public:
    InputArchive(const char* filename)
      : is(filename, ios_base::binary | ios_base::in) {};

    template <class T>
    void load(T* t, long n) const {
      is.read((char *) t, n*sizeof(T));
    }
  };
  \endcode
*/

#include <complex>
#include <iostream>
#include <cstdio>
#include <vector>
#include <map>
#include <world/typestuff.h>
//#include <world/worldprofile.h>
#include <world/enable_if.h>
#include <world/worldexc.h>

#define ARCHIVE_COOKIE "archive"
#define ARCHIVE_MAJOR_VERSION 0
#define ARCHIVE_MINOR_VERSION 1

//#define MAD_ARCHIVE_DEBUG_ENABLE

#ifdef MAD_ARCHIVE_DEBUG_ENABLE
#define MAD_ARCHIVE_DEBUG(s) s
//using std::endl;
#else
#define MAD_ARCHIVE_DEBUG(s)
#endif

namespace madness {


    // Forward declarations
    template <typename T> class Tensor;

    namespace archive {

        // Forward declarations
        template <class>
        class archive_array;
        template <class T>
        inline archive_array<T> wrap(const T*, unsigned int);
        template <class T>
        inline archive_array<unsigned char> wrap_opaque(const T*, unsigned int);
        template <class T>
        inline archive_array<unsigned char> wrap_opaque(const T&);

        // There are 64 empty slots for user types.  Free space for
        // registering user types begins at cookie=128.

#ifdef MAD_ARCHIVE_TYPE_NAMES_CC
        const char *archive_type_names[256];
#else
        extern const char *archive_type_names[256];
#endif
        void archive_initialize_type_names();

        // Used to enable type checking inside archives
        template <typename T>
        struct archive_typeinfo {
            static const unsigned char cookie = 255; ///< 255 indicates unknown type
        };

        // Returns the name of the type, or unknown if not registered.
        template <typename T>
        const char* get_type_name() {
            return archive_type_names[archive_typeinfo<T>::cookie];
        }


#if defined(ARCHIVE_REGISTER_TYPE_INSTANTIATE_HERE) && defined(ARCHIVE_REGISTER_TYPE_IBMBUG)
#define ARCHIVE_REGISTER_TYPE_XLC_EXTRA(T) \
        ; const unsigned char archive_typeinfo< T >::cookie
#else
#define ARCHIVE_REGISTER_TYPE_XLC_EXTRA(T)
#endif

        /// \def ARCHIVE_REGISTER_TYPE(T, cooky)
        /// \brief Used to associate type with cookie value inside archive
        ///
        /// \def  ARCHIVE_REGISTER_TYPE_AND_PTR(T, cooky)
        /// \brief Used to associate type and ptr to type with cookie value inside archive
        ///
        /// \def ARCHIVE_REGISTER_TYPE_NAME(T)
        /// \brief Used to associate names with types
        ///
        /// \def ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(T)
        /// \brief Used to associate names with types and pointers

#define ARCHIVE_REGISTER_TYPE(T, cooky) \
	template <> struct archive_typeinfo< T > { \
		static const unsigned char cookie = cooky; \
	} \
        ARCHIVE_REGISTER_TYPE_XLC_EXTRA(T)

#define ARCHIVE_REGISTER_TYPE_AND_PTR(T, cooky) \
        ARCHIVE_REGISTER_TYPE(T, cooky); \
        ARCHIVE_REGISTER_TYPE(T*, cooky+64)

#define ATN ::madness::archive::archive_type_names
#define ATI ::madness::archive::archive_typeinfo
#define ARCHIVE_REGISTER_TYPE_NAME(T) \
     if (strcmp(ATN[ATI< T >::cookie],"invalid")) {\
        std::cout << "archive_register_type_name: slot/cookie already in use! "<< #T << " " << ATN[ATI< T >::cookie]<< std::endl; \
        MADNESS_EXCEPTION("archive_register_type_name: slot/cookie already in use!", 0); \
     } \
     ATN[ATI< T >::cookie] = #T

#define ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(T) \
     ARCHIVE_REGISTER_TYPE_NAME(T); \
     ARCHIVE_REGISTER_TYPE_NAME(T*)

        ARCHIVE_REGISTER_TYPE_AND_PTR(unsigned char,0);
        ARCHIVE_REGISTER_TYPE_AND_PTR(unsigned short,1);
        ARCHIVE_REGISTER_TYPE_AND_PTR(unsigned int,2);
        ARCHIVE_REGISTER_TYPE_AND_PTR(unsigned long,3);
        ARCHIVE_REGISTER_TYPE_AND_PTR(unsigned long long,4);
        ARCHIVE_REGISTER_TYPE_AND_PTR(signed char,5);
        ARCHIVE_REGISTER_TYPE_AND_PTR(char,5);	// Needed, but why?
        ARCHIVE_REGISTER_TYPE_AND_PTR(signed short,6);
        ARCHIVE_REGISTER_TYPE_AND_PTR(signed int,7);
        ARCHIVE_REGISTER_TYPE_AND_PTR(signed long,8);
        ARCHIVE_REGISTER_TYPE_AND_PTR(signed long long,9);
        ARCHIVE_REGISTER_TYPE_AND_PTR(bool,10);
        ARCHIVE_REGISTER_TYPE_AND_PTR(float,11);
        ARCHIVE_REGISTER_TYPE_AND_PTR(double,12);
        ARCHIVE_REGISTER_TYPE_AND_PTR(long double,13);
        ARCHIVE_REGISTER_TYPE_AND_PTR(std::complex<float>,14);
        ARCHIVE_REGISTER_TYPE_AND_PTR(std::complex<double>,15);

        ARCHIVE_REGISTER_TYPE_AND_PTR(std::vector<char>,20);
        ARCHIVE_REGISTER_TYPE_AND_PTR(std::vector<unsigned char>,21);
        ARCHIVE_REGISTER_TYPE_AND_PTR(std::vector<short>,22);
        ARCHIVE_REGISTER_TYPE_AND_PTR(std::vector<unsigned short>,23);
        ARCHIVE_REGISTER_TYPE_AND_PTR(std::vector<int>,24);
        ARCHIVE_REGISTER_TYPE_AND_PTR(std::vector<unsigned int>,25);
        ARCHIVE_REGISTER_TYPE_AND_PTR(std::vector<long>,26);
        ARCHIVE_REGISTER_TYPE_AND_PTR(std::vector<unsigned long>,27);
        ARCHIVE_REGISTER_TYPE_AND_PTR(std::vector<bool>,28);
        ARCHIVE_REGISTER_TYPE_AND_PTR(std::vector<float>,29);
        ARCHIVE_REGISTER_TYPE_AND_PTR(std::vector<double>,30);

        ARCHIVE_REGISTER_TYPE_AND_PTR(std::string,31);

        ARCHIVE_REGISTER_TYPE_AND_PTR(Tensor<int>,32);
        ARCHIVE_REGISTER_TYPE_AND_PTR(Tensor<long>,33);
        ARCHIVE_REGISTER_TYPE_AND_PTR(Tensor<float>,34);
        ARCHIVE_REGISTER_TYPE_AND_PTR(Tensor<double>,35);
        ARCHIVE_REGISTER_TYPE_AND_PTR(Tensor< std::complex<float> >,36);
        ARCHIVE_REGISTER_TYPE_AND_PTR(Tensor< std::complex<double> >,37);

        /// Base class for all archives classes
        class BaseArchive {
        public:
            static const bool is_archive = true;
            static const bool is_input_archive = false;
            static const bool is_output_archive = false;
            static const bool is_parallel_archive = false;
            BaseArchive() {
                archive_initialize_type_names();
            }
        }; // class BaseArchive

        /// Base class for input archives classes
        class BaseInputArchive : public BaseArchive {
        public:
            static const bool is_input_archive = true;
        }; // class BaseInputArchive


        /// Base class for output archives classes
        class BaseOutputArchive : public BaseArchive {
        public:
            static const bool is_output_archive = true;
        }; // class BaseOutputArchive

        /// Checks that \c T is an archive type

        /// If \c T is an archive type then \c is_archive will be inherited from
        /// \c std::true_type , otherwise it is inherited from
        /// \c std::false_type .
        template <typename T>
        struct is_archive : public std::is_base_of<BaseArchive, T> {};

        /// Checks that \c T is an input archive type

        /// If \c T is an input archive type then \c is_archive will be
        /// inherited from \c std::true_type , otherwise it is inherited from
        /// \c std::false_type .
        template <typename T>
        struct is_input_archive : public std::is_base_of<BaseInputArchive, T> {};

        /// Checks that \c T is an output archive type

        /// If \c T is an output archive type then \c is_archive will be
        /// inherited from \c std::true_type , otherwise it is inherited from
        /// \c std::false_type .
        template <typename T>
        struct is_output_archive : public std::is_base_of<BaseOutputArchive, T> {};

        // Serialize an array of fundamental stuff
        template <class Archive, class T>
        typename enable_if_c< is_serializable<T>::value && is_output_archive<Archive>::value >::type
        serialize(const Archive& ar, const T* t, unsigned int n) {
            MAD_ARCHIVE_DEBUG(std::cout << "serialize fund array" << std::endl);
            ar.store(t,n);
        }


        // Deserialize an array of fundamental stuff
        template <class Archive, class T>
        typename enable_if_c< is_serializable<T>::value && is_input_archive<Archive>::value >::type
        serialize(const Archive& ar, const T* t, unsigned int n) {
            MAD_ARCHIVE_DEBUG(std::cout << "deserialize fund array" << std::endl);
            ar.load((T*) t,n);
        }


        // (de)Serialize an array of non-fundamental stuff
        template <class Archive, class T>
        typename enable_if_c< ! is_serializable<T>::value && is_archive<Archive>::value >::type
        serialize(const Archive& ar, const T* t, unsigned int n) {
            MAD_ARCHIVE_DEBUG(std::cout << "(de)serialize non-fund array" << std::endl);
            for (unsigned int i=0; i<n; ++i) ar & t[i];
        }


        /// Default implementation of pre/postamble
        template <class Archive, class T>
        struct ArchivePrePostImpl {
            static inline void preamble_load(const Archive& ar) {
                unsigned char ck = archive_typeinfo<T>::cookie;
                unsigned char cookie;
                ar.load(&cookie, 1); // cannot use >>
                if (cookie != ck) {
                    char msg[255];
                    std::sprintf(msg,"InputArchive type mismatch: expected cookie "
                                 "%u (%s) but got %u (%s) instead",
                                 ck, archive_type_names[ck],
                                 cookie,archive_type_names[cookie]);
                    std::cerr << msg << std::endl;
                    MADNESS_EXCEPTION(msg, static_cast<int>(cookie));
                }
                else {
                    MAD_ARCHIVE_DEBUG(std::cout << "read cookie " << archive_type_names[cookie] << std::endl);
                }
            }


            /// Serialize a cookie for type checking
            static inline void preamble_store(const Archive& ar) {
                unsigned char ck = archive_typeinfo<T>::cookie;
                ar.store(&ck, 1); // cannot use <<
                MAD_ARCHIVE_DEBUG(std::cout << "wrote cookie " << archive_type_names[ck] << std::endl);
            }


            /// By default there is no postamble
            static inline void postamble_load(const Archive& /*ar*/) {}

            /// By default there is no postamble
            static inline void postamble_store(const Archive& /*ar*/) {}
        };


        /// Default symmetric serialization of a non-fundamental thingy
        template <class Archive, class T>
        struct ArchiveSerializeImpl {
            static inline void serialize(const Archive& ar, T& t) {
                t.serialize(ar);
            }
        };


        // Redirect \c serialize(ar,t) to \c serialize(ar,&t,1) for fundamental types
        template <class Archive, class T>
        inline
        typename enable_if_c< is_serializable<T>::value && is_archive<Archive>::value >::type
        serialize(const Archive& ar, const T& t) {
            MAD_ARCHIVE_DEBUG(std::cout << "serialize(ar,t) -> serialize(ar,&t,1)" << std::endl);
            serialize(ar,&t,1);
        }


        // Redirect \c serialize(ar,t) to \c ArchiveSerializeImpl for non-fundamental types
        template <class Archive, class T>
        inline
        typename enable_if_c< !is_serializable<T>::value && is_archive<Archive>::value >::type
        serialize(const Archive& ar, const T& t) {
            MAD_ARCHIVE_DEBUG(std::cout << "serialize(ar,t) -> ArchiveSerializeImpl" << std::endl);
            ArchiveSerializeImpl<Archive,T>::serialize(ar,(T&) t);
        }


        /// Default store of a thingy via serialize(ar,t)
        template <class Archive, class T>
        struct ArchiveStoreImpl {
            static inline void store(const Archive& ar, const T& t) {
                MAD_ARCHIVE_DEBUG(std::cout << "store(ar,t) default" << std::endl);
                serialize(ar,t);
            }
        };


        /// Default load of a thingy via serialize(ar,t)
        template <class Archive, class T>
        struct ArchiveLoadImpl {
            static inline void load(const Archive& ar, const T& t) {
                MAD_ARCHIVE_DEBUG(std::cout << "load(ar,t) default" << std::endl);
                serialize(ar,t);
            }
        };


        /// Default implementation of wrap_store and wrap_load
        template <class Archive, class T>
        struct ArchiveImpl {
            static inline const Archive& wrap_store(const Archive& ar, const T& t) {
                MAD_ARCHIVE_DEBUG(std::cout << "wrap_store for default" << std::endl);
                ArchivePrePostImpl<Archive,T>::preamble_store(ar);
                ArchiveStoreImpl<Archive,T>::store(ar,t);
                ArchivePrePostImpl<Archive,T>::postamble_store(ar);
                return ar;
            }

            static inline const Archive& wrap_load(const Archive& ar, const T& t) {
                MAD_ARCHIVE_DEBUG(std::cout << "wrap_load for default" << std::endl);
                ArchivePrePostImpl<Archive,T>::preamble_load(ar);
                ArchiveLoadImpl<Archive,T>::load(ar,(T&) t);  // Loses constness here!
                ArchivePrePostImpl<Archive,T>::postamble_load(ar);
                return ar;
            }
        };


        // Redirect \c << to ArchiveImpl::wrap_store for output archives
        template <class Archive, class T>
        inline
        typename enable_if<is_output_archive<Archive>, const Archive&>::type
        operator<<(const Archive& ar, const T& t) {
            //PROFILE_FUNC;
            return ArchiveImpl<Archive,T>::wrap_store(ar,t);
        }

        // Redirect \c >> to ArchiveImpl::wrap_load for input archives
        template <class Archive, class T>
        inline
        typename enable_if<is_input_archive<Archive>, const Archive&>::type
        operator>>(const Archive& ar, const T& t) {
            //PROFILE_FUNC;
            return ArchiveImpl<Archive,T>::wrap_load(ar,t);
        }

        // Redirect \c & to ArchiveImpl::wrap_store for output archives
        template <class Archive, class T>
        inline
        typename enable_if<is_output_archive<Archive>, const Archive&>::type
        operator&(const Archive& ar, const T& t) {
            //PROFILE_FUNC;
            return ArchiveImpl<Archive,T>::wrap_store(ar,t);
        }

        // Redirect \c & to ArchiveImpl::wrap_load for input archives
        template <class Archive, class T>
        inline
        typename enable_if<is_input_archive<Archive>, const Archive&>::type
        operator&(const Archive& ar, const T& t) {
            //PROFILE_FUNC;
            return ArchiveImpl<Archive,T>::wrap_load(ar,t);
        }


        ///////////////////////////////////////////////////////////////

        /// Wrapper for opaque pointer ... bitwise copy of the pointer ... no remapping performed
        template <class T>
        class archive_ptr {
        public:
            T* ptr;

            archive_ptr(T* t = 0) : ptr(t) {}

            T& operator*(){return *ptr;}

            template <class Archive>
            void serialize(const Archive& ar) {ar & wrap_opaque(&ptr, 1);}
        };

	template <class T>
	inline
	archive_ptr<T> wrap_ptr(T* p) {return archive_ptr<T>(p);}

        /// Wrapper for dynamic arrays and pointers
        template <class T>
        class archive_array {
        public:
            const T* ptr;
            unsigned int n;

            archive_array(const T *ptr, unsigned int n) : ptr(ptr), n(n) {}

            archive_array() : ptr(0), n(0) {}
        };


        /// Factory function to wrap dynamically allocated pointer as typed archive_array
        template <class T>
        inline
        archive_array<T> wrap(const T* ptr, unsigned int n) {
            return archive_array<T>(ptr,n);
        }

        /// Factory function to wrap pointer to contiguous data as opaque (uchar) archive_array
        template <class T>
        inline
        archive_array<unsigned char> wrap_opaque(const T* ptr, unsigned int n) {
            return archive_array<unsigned char>((unsigned char*) ptr, n*sizeof(T));
        }

        /// Factory function to wrap contiguous scalar as opaque (uchar) archive_array
        template <class T>
        inline
        archive_array<unsigned char> wrap_opaque(const T& t) {
            return archive_array<unsigned char>((unsigned char*) &t,sizeof(t));
        }

        /// Partial specialization for archive_array

        /// This makes use of stuff that user specializations need not
        template <class Archive, class T>
        struct ArchiveImpl< Archive, archive_array<T> > {
            static inline const Archive& wrap_store(const Archive& ar, const archive_array<T>& t) {
                MAD_ARCHIVE_DEBUG(std::cout << "wrap_store for archive_array" << std::endl);
                ArchivePrePostImpl<Archive,T*>::preamble_store(ar);
                //ar << t.n;
                //ArchivePrePostImpl<Archive,T>::preamble_store(ar);
                serialize(ar,(T *) t.ptr,t.n);
                //ArchivePrePostImpl<Archive,T>::postamble_store(ar);
                ArchivePrePostImpl<Archive,T*>::postamble_store(ar);
                return ar;
            }

            static inline const Archive& wrap_load(const Archive& ar, const archive_array<T>& t) {
                MAD_ARCHIVE_DEBUG(std::cout << "wrap_load for archive_array" << std::endl);
                ArchivePrePostImpl<Archive,T*>::preamble_load(ar);
                //unsigned int n;
                //ar >> n;
                //if (n != t.n)
                //    MADNESS_EXCEPTION("deserializing archive_array: dimension mismatch", n);
                //ArchivePrePostImpl<Archive,T>::preamble_load(ar);
                serialize(ar,(T *) t.ptr,t.n);
                //ArchivePrePostImpl<Archive,T>::postamble_load(ar);
                ArchivePrePostImpl<Archive,T*>::postamble_load(ar);
                return ar;
            }
        };


        /// Partial specialization for fixed dimension array redirects to archive_array
        template <class Archive, class T, std::size_t n>
        struct ArchiveImpl<Archive, T[n]> {
            static inline const Archive& wrap_store(const Archive& ar, const T(&t)[n]) {
                MAD_ARCHIVE_DEBUG(std::cout << "wrap_store for array" << std::endl);
                ar << wrap(&t[0],n);
                return ar;
            }

            static inline const Archive& wrap_load(const Archive& ar, const T(&t)[n]) {
                MAD_ARCHIVE_DEBUG(std::cout << "wrap_load for array" << std::endl);
                ar >> wrap(&t[0],n);
                return ar;
            }
        };


        /// Serialize a complex number
        template <class Archive, typename T>
        struct ArchiveStoreImpl< Archive, std::complex<T> > {
            static inline void store(const Archive& ar, const std::complex<T>& c) {
                MAD_ARCHIVE_DEBUG(std::cout << "serialize complex number" << std::endl);
                ar & c.real() & c.imag();
            }
        };


        /// Deserialize a complex number
        template <class Archive, typename T>
        struct ArchiveLoadImpl< Archive, std::complex<T> > {
            static inline void load(const Archive& ar, std::complex<T>& c) {
                MAD_ARCHIVE_DEBUG(std::cout << "deserialize complex number" << std::endl);
                T r, i;
                ar & r & i;
                c = std::complex<T>(r,i);
            }
        };


        /// Serialize STL vector.
        template <class Archive, typename T>
        struct ArchiveStoreImpl< Archive, std::vector<T> > {
            static inline void store(const Archive& ar, const std::vector<T>& v) {
                MAD_ARCHIVE_DEBUG(std::cout << "serialize STL vector" << std::endl);
                ar & v.size();
                ar & wrap(&v[0],v.size());
            }
        };


        /// Deserialize STL vector. Clears & resizes as necessary.
        template <class Archive, typename T>
        struct ArchiveLoadImpl< Archive, std::vector<T> > {
            static void load(const Archive& ar, std::vector<T>& v) {
                MAD_ARCHIVE_DEBUG(std::cout << "deserialize STL vector" << std::endl);
                std::size_t n;
                ar & n;
                if (n != v.size()) {
                    v.clear();
                    v.resize(n);
                }
                ar & wrap((T *) &v[0],n);
            }
        };

        /// Serialize STL vector<bool> (as plain array of bool)
        template <class Archive>
        struct ArchiveStoreImpl< Archive, std::vector<bool> > {
            static inline void store(const Archive& ar, const std::vector<bool>& v) {
                MAD_ARCHIVE_DEBUG(std::cout << "serialize STL vector<bool>" << std::endl);
                std::size_t n = v.size();
                bool* b = new bool[n];
                for (std::size_t i=0; i<n; ++i) b[i] = v[i];
                ar & n & wrap(b,v.size());
                delete [] b;
            }
        };


        /// Deserialize STL vector<bool>. Clears & resizes as necessary.
        template <class Archive>
        struct ArchiveLoadImpl< Archive, std::vector<bool> > {
            static void load(const Archive& ar, std::vector<bool>& v) {
                MAD_ARCHIVE_DEBUG(std::cout << "deserialize STL vector" << std::endl);
                std::size_t n;
                ar & n;
                if (n != v.size()) {
                    v.clear();
                    v.resize(n);
                }
                bool* b = new bool[n];
                ar & wrap(b,v.size());
                for (std::size_t i=0; i<n; ++i) v[i] = b[i];
                delete [] b;
            }
        };

        /// Serialize STL string
        template <class Archive>
        struct ArchiveStoreImpl< Archive, std::string > {
            static void store(const Archive& ar, const std::string& v) {
                MAD_ARCHIVE_DEBUG(std::cout << "serialize STL string" << std::endl);
                ar & v.size();
                ar & wrap((const char*) &v[0],v.size());
            }
        };


        /// Deserialize STL string. Clears & resizes as necessary.
        template <class Archive>
        struct ArchiveLoadImpl< Archive, std::string > {
            static void load(const Archive& ar, std::string& v) {
                MAD_ARCHIVE_DEBUG(std::cout << "deserialize STL string" << std::endl);
                std::size_t n;
                ar & n;
                if (n != v.size()) {
                    v.clear();
                    v.resize(n);
                }
                ar & wrap((char*) &v[0],n);
            }
        };


        /// (de)Serialize an STL pair.
        template <class Archive, typename T, typename Q>
        struct ArchiveSerializeImpl< Archive, std::pair<T,Q> > {
            static inline void serialize(const Archive& ar, std::pair<T,Q>& t) {
                MAD_ARCHIVE_DEBUG(std::cout << "(de)serialize STL pair" << std::endl);
                ar & t.first & t.second;
            }
        };


        /// Serialize an STL map (crudely).
        template <class Archive, typename T, typename Q>
        struct ArchiveStoreImpl< Archive, std::map<T,Q> > {
            static void store(const Archive& ar, const std::map<T,Q>& t) {
                MAD_ARCHIVE_DEBUG(std::cout << "serialize STL map" << std::endl);
                ar << t.size();
                for (typename std::map<T,Q>::const_iterator p = t.begin();
                        p != t.end(); ++p) {
                    // Fun and games here since IBM's iterator (const or
                    // otherwise) gives us a const qualified key
                    // (p->first) which buggers up the type matching
                    // unless the user defines pair(T,Q) and pair(const
                    // T,Q) to have cookie (which is tedious).
                    std::pair<T,Q> pp = *p;
                    ar & pp;
                }
            }
        };


        /// Deserialize an STL map.  Map is NOT cleared; duplicate elements are replaced.
        template <class Archive, typename T, typename Q>
        struct ArchiveLoadImpl< Archive, std::map<T,Q> > {
            static void load(const Archive& ar, std::map<T,Q>& t) {
                MAD_ARCHIVE_DEBUG(std::cout << "deserialize STL map" << std::endl);
                std::size_t n;
                ar & n;
                while (n--) {
                    std::pair<T,Q> p;
                    ar & p;
                    t[p.first] = p.second;
                }
            }
        };



    }
}

#endif // MADNESS_WORLD_ARCHIVE_H__INCLUDED
