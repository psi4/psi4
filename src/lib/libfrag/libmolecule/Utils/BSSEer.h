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
#ifndef SRC_LIB_LIBFRAG_LIBMOLECULE_UTILS_BSSEER_H_
#define SRC_LIB_LIBFRAG_LIBMOLECULE_UTILS_BSSEER_H_
#include "../LibFragFragment.h"

namespace psi{
namespace LibMolecule{
class Molecule;
class NMers;
class BSSEer{
   public:
      virtual void CalcBSSE(NMers& Sys,unsigned int Stop,unsigned int Start)const=0;
      ///Shut-up compiler....
      virtual ~BSSEer(){}
};


class FullBSSEer:public BSSEer{
   public:
      void CalcBSSE(NMers& Sys,unsigned int Stop,unsigned int Start)const;
};
/** \brief Little wrapper class to return the next SN for a VMFC(n) calculation
 *
 *   We assume we are given some serial number with m real
 *   fragments, and n-m ghosts (the top-level initialization should
 *   contain 0 ghosts; this description is for an arbitrary nested
 *   iterator within that top-level version).  This SN looks like:
 *   \verbatim
 *         (m+1)...n
 *   IJ...m
 *   \endverbatim
 *   where ghosts are on the superscript.  Internally ghosts
 *   are flagged as being negative fragments.  Upon creation
 *   we promote the m-th index to a ghost, reducing the number of
 *   real fragments by one:
 *   \verbatim
 *             m...n
 *   IJ...(m-1)
 *   \endverbatim
 *   This particular iterator then passes the resulting SN into
 *   a subiterator, where the entire process recurses until:
 *   \verbatim
 *    J...n
 *   I
 *   \endverbatim
 *   This gives us our initial state upon instance creation.
 *
 *   When the user calls operator()++ each iterator returns it's value of
 *   Current_, if it hasn't done so already, or increments it's iterator
 *   and returns that iterators value.  Hence we are always be given the
 *   shallowest iterator's value that hasn't been returned yet.  Once our
 *   sub iterator has finished iterating we implement our value, and
 *   re-initialize our sub-iterator; the process then repeats until all
 *   iterators return done.  As with most iterators attempting to
 *   increment this class once it has finished leads to undefined behavior.
 *
 *
 *   Ultimately this produces a sequence analogous to (assuming a trimer
 *   IJK):
 *   \verbatim
 *   I J -K
 *   I -J -K
 *   J -I -K
 *   I K -J
 *   I -J -K
 *   K -I -J
 *   J K -I
 *   J -I -K
 *   K -I -J
 *   \endverbatim
 *   note that some terms are repeated, if this is not the desired
 *   operation set the flag in the constructor to false.
 */
class VMFCnItr{
   private:
      ///The SN currently stored in this generator
      SerialNumber Mine_;
      ///The current SN, what will be returned
      SerialNumber Current_;
      ///The real fragments in the intialized SN
      SerialNumber Real_;
      ///The ghost fragments in the initialized SN
      SerialNumber Ghost_;
      ///An iterator to Current_ w/one additional ghost monomer
      boost::shared_ptr<VMFCnItr> MyItr_;
      ///A flag to symbolize when this iterator has completed
      bool IsDone_;
      ///A flag to symbolize whether the value in Current_ has been used
      bool UsedMyValue_;
      ///A flag symbolizing whether the user wants the entire list
      bool AllowDuplicates_;
      ///A map to keep track of the SNs we have found
      std::set<SerialNumber> FoundSNs_;
      ///The index we are updating
      unsigned int l_;
      ///Function that sets Current_ to it's next value
      void Next();
      ///Updates Current, so that index l of the original SN is now a ghost
      void UpdateMe(unsigned int l);
   public:
      ///Returns the current serial number
      const SerialNumber& operator*()const{return Current_;}
      ///Allows access to the current serial number's members
      const SerialNumber* operator->()const{return &Current_;}
      ///SN must contain no ghosts for the top-level iterator
      VMFCnItr(const SerialNumber& SN,const bool AllowDuplicates_=true);
      ///Returns true if this iterator has completed
      bool Done()const{return IsDone_;}
      ///Increments this iterator
      const VMFCnItr& operator++(){Next();return *this;}

};

class VMFCn:public BSSEer{
   public:
      void CalcBSSE(NMers& Sys,unsigned int Stop,unsigned int Start)const;
};

}}//End namespace


#endif /* SRC_LIB_LIBFRAG_LIBMOLECULE_UTILS_BSSEER_H_ */
