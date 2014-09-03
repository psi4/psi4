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
#ifndef BSSEER_H_
#define BSSEER_H_

#include "AtomSet.h"
#include "MBEFrag.h"

namespace psi{
namespace LibFrag{

typedef std::vector<Set<SharedGhost> > GhostType;

class BSSEer{
   protected:
      ///The total number of atoms in the entire system
      int NAtoms_;
      void CommonInit(GhostType& Ghosts,NMerSet& NMers);
      virtual void BSSEImpl(GhostType& Ghosts,NMerSet& NMers)=0;
   public:
      BSSEer(int natoms):NAtoms_(natoms){}
      virtual ~BSSEer(){}
      void AddBSSEJobs(NMerSet& NMers);
};

///Class that does nothing, removes check to see if factory exists
class NullBSSE:public BSSEer{
   protected:
      void BSSEImpl(GhostType& Ghosts,NMerSet& NMers){}
   public:
      NullBSSE():BSSEer(0){}
};

class FullBSSE:public BSSEer{
   protected:
      void BSSEImpl(GhostType& Ghosts,NMerSet& NMers);
   public:
      FullBSSE(int natoms):BSSEer(natoms){}
      ~FullBSSE(){}
};


}}//End namespaces



#endif /* BSSEER_H_ */
