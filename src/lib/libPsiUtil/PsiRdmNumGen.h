/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
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
#ifndef SRC_LIB_LIBPSIUTIL_PSIRDMNUMGET_H_
#define SRC_LIB_LIBPSIUTIL_PSIRDMNUMGET_H_
#include <boost/shared_ptr.hpp>

#include <boost/random.hpp>//For our random stuff
#include <limits>//For maximum integer size
#include <ctime>//For our original seed
namespace psi{

///Wrapper class that ensures that a particular engine is only made once
template<typename T>
class PsiRdmNumEngine{
   public:
      ///The actual engine (G++ seems to think that this needs to be public);
      static T Engine_;
};

template<typename T>
T PsiRdmNumEngine<T>::Engine_(std::time(0));


/** \brief A class to provide better random numbers than srand() & rand()
 *
 *  Terms:
 *  -Generator: An object that can return a random number (is the result
 *              of an engine and a distribution)
 *  -Engine:    An algorithm that returns a number given a seed.
 *              Algorithms are chosen such that two seeds that are
 *              close together provide very different results
 *  -Distribution: An object that takes the value of an engine and
 *                 imposes certain restrictions on it (e.g. makes
 *                 it between 0 and 1, etc.), the result of the
 *                 distribution is the random number
 *
 *  If you want to get a bunch of random numbers one after another
 *  you can't do the traditional:
 *  \code
 *  srand(time(NULL));
 *  rand();
 *  \endcode
 *
 *  because time has a 1 second resolution, meaning the seed is set
 *  the same each time.  The traditional solution, you will find
 *  online, is to ensure srand() only gets called once.  This means
 *  you can't put it in an object's constructor and then make a
 *  bunch of those objects.  This is what bit me in the arse and
 *  caused me to make this object.  To avoid this I make global
 *  objects for each engine type (wrapped in the little class that
 *  proceeds this class).
 *
 *  Right now the object calls boost's random library, but in
 *  the future we may want to migrate this over to the c++11 random
 *  features.  Either library is capable of generating distributions, and
 *  also using multiple algorithms for creating the numbers.
 *
 *  Both the boost and c++11 random features work by pairing an engine
 *  with a distribution, the result is a generator.  Again, both
 *  standards say that a generator is a functor and calling the ()
 *  operator should return the number.  This class adheres to that
 *  standard as well, and that is how you will get your number out.
 *
 *  For the defaults, I don't care about being truly random so I have
 *  chosen the default engine as boost::mt19937, and my distribution as
 *  boost::uniform_int.  I don't particularly know what that engine is,
 *  but it's the one recommended by boost for being a balance between
 *  randomness, cost, and size. The distribution, on the other hand,
 *  was specifically chosen, as it generates integers uniformly between
 *  Max, and Min (inclusive, note all ranges are inclusive).
 *
 *
 */
template<typename T0=int,
         typename T1=boost::random::mt19937,
         typename T2=boost::random::uniform_int_distribution<> >
class PsiRdmNumGen{
   private:
      PsiRdmNumEngine<T1> Engine_;

      ///Convenient typedef of the current type
      typedef PsiRdmNumGen<T0,T1,T2> ThisType_;

      ///... of the generator
      typedef boost::variate_generator<T1&,T2> GenType_;

      ///.... of a shared pointer to the generator
      typedef boost::shared_ptr<GenType_> SharedGen_;

      ///The object that actually generates the random number
      SharedGen_ Generator_;

      ///Shallow copies the Generator
      void Copy(const ThisType_& other){
         this->Generator_=other.Generator_;
      }

      ///Abstraction of initialization, by type of returned number
      SharedGen_ Init(const T0& Max, const T0& Min){
         T2 dist(Min,Max);
        SharedGen_ temp(new GenType_(Engine_.Engine_,dist));
        return temp;
      }

   public:

      ///Creates a generator that will pick numbers between [0,Max]
      PsiRdmNumGen<T0,T1,T2>(
            const T0 Max=std::numeric_limits<T0>::max(),
            const T0 Min=0){
         Generator_=Init(Max,Min);
      }

      ///Creates by shallow copy
      PsiRdmNumGen<T0,T1,T2>(const ThisType_& other){
         this->Copy(other);
      }

      ///Assignment..
      const PsiRdmNumGen<T0,T1,T2> operator=(const ThisType_& other){
         if(this!=&other)this->Copy(other);
         return *this;
      }

      ///Returns a random number
      int operator()(){return (*Generator_)();}
};


}//End psi namespace
#endif /* SRC_LIB_LIBPSIUTIL_PSIRDMNUMGET_H_ */