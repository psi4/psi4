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
#ifndef EMBEDDER_H_
#define EMBEDDER_H_
#include "LibFragTypes.h"
#include "libmints/molecule.h"

namespace psi{
namespace LibFrag{

typedef std::vector<Set<SharedCharge> > ChargeType;

class AtomSet;
class Embedder{
   protected:
      ///Are we iterating
      bool DoIterate_;
      ///A vector of all charges for the supersystem
      std::vector<double> Charges_;
      std::vector<double>Carts_;
      ///Did we just get the point charges
      bool NotFirstItr_;

      virtual void EmbedImpl(ChargeType& ChargesBySet, NMerSet& Set2Embed)=0;
      void CommonInit(ChargeType& ChargesBySet);
   public:
      Embedder(SharedMol& AMol,bool Iterating);
      bool HaveCharges(){return NotFirstItr_;}
      virtual ~Embedder(){}
      /** \brief Tells you whether you need to rerun the calculations
       *
       *  There are two reasons why one may need to rerun the calculations:
       *  first and the most common, we just ran the initial calculations
       *  that generated the density/charges.  Or two, we have already
       *  iterated for one, and the user wants us to iterate to convergence.
       *  Either way, this fxn returns the appropriate response.
       *
       */
      bool Iterate(const int itr);
      ///Only call set charge once you are ready to set all the charges
      void SetCharge(int i,double q);
      virtual void Embed(NMerSet& Set2Embed);
      void print_out();

};

class APCEmbedder:public Embedder{
   protected:
      void EmbedImpl(ChargeType& ChargesBySet, NMerSet& Set2Embed);

   public:
      APCEmbedder(SharedMol& AMol,bool Iterating):
         Embedder(AMol,Iterating){}


};

}}//End namespaces


#endif /* EMBEDDER_H_ */
