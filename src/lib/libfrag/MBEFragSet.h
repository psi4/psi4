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
 */
class MBEFragSet{
   private:
      typedef boost::shared_ptr<Fragmenter> SharedFragFac;
      typedef boost::shared_ptr<Capper> SharedCapFac;
      typedef boost::shared_ptr<BSSEer> SharedBSSEFac;
      typedef boost::shared_ptr<Embedder> SharedEmbedFac;
      typedef boost::shared_ptr<LibFragOptions> SharedOptions;

      SharedEmbedFac EmbedFactory_;
      SharedBSSEFac BSSEFactory_;
      SharedCapFac CapFactory_;
      SharedFragFac FragFactory_;

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

      ///Constructor for making n-mers
      MBEFragSet(const MBEFragSet& Monomers,const int N);

      boost::shared_ptr<const MBEFrag> operator[](const int i){
         return Frags_[i];
      }
};


}
}      //End namespaces



#endif /* MBEFRAGSET_H_ */
