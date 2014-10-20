/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */
#ifndef MBEFRAGSET_H_
#define MBEFRAGSET_H_

#include <boost/shared_ptr.hpp>
#include <vector>
#include "Fragmenter.h"

namespace psi{
namespace LibFrag{
class MBEFrag;
class Embedder;
class Capper;
class BSSEer;
class LibFragOptions;

/** \brief The MBEFragSet Class is intended to be the fundamental object
 *         of the LibFrag library.
 *
 *         This class is essentially a container class aimed at simplifying
 *         fragment manipulations.  One thing to note at this point is that
 *         it is friends with a couple classes, similar to how the Frag class
 *         is, and for essentially the same reason.  Basically, we have
 *         "factories", like the symmetrizer that are responsible for setting
 *         up the Frag Set.
 *
 */
class MBEFragSet{
   private:
      friend class Symmetrizer;
      typedef boost::shared_ptr<Fragmenter> SharedFragFac;
      typedef boost::shared_ptr<Capper> SharedCapFac;
      typedef boost::shared_ptr<BSSEer> SharedBSSEFac;
      typedef boost::shared_ptr<Embedder> SharedEmbedFac;
      typedef boost::shared_ptr<LibFragOptions> SharedOptions;

      SharedEmbedFac EmbedFactory_;
      SharedBSSEFac BSSEFactory_;
      SharedCapFac CapFactory_;
      SharedFragFac FragFactory_;

      void Copy(const MBEFragSet& other);
   protected:
      std::vector<boost::shared_ptr<MBEFrag> > Frags_;
      FragProps Properties_;


   public:


      //Hackzzz!!!
      boost::shared_ptr<Embedder> EmbedFactory();

      void Embed();

      bool Disjoint()const{return Properties_.Disjoint_;}

      bool Severed()const{return Properties_.Severed_;}

      int size()const{return Frags_.size();}

      ///Constructor for making fragments
      MBEFragSet(SharedOptions& Options,SharedMol& AMol);

      ///Constructor for making n-mers (note don't call with N=1)
      MBEFragSet(const MBEFragSet& Monomers,const int N);

      ///Copy constructor
      MBEFragSet(const MBEFragSet& other){this->Copy(other);}

      ///Assignment operator
      const MBEFragSet& operator=(const MBEFragSet& other){
         if(this!=&other)this->Copy(other);
         return *this;
      }

      boost::shared_ptr<const MBEFrag> operator[](const int i)const{
         return Frags_[i];
      }
};


}
}      //End namespaces



#endif /* MBEFRAGSET_H_ */
