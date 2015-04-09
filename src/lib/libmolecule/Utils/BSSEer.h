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
#include "PowerSetItr.h"
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
 *   Assume we have some serial number, whose indices form the set:
 *   \f$S=\lbrace i\ j\ k\ l\rbrace\f$.  What we want is then to
 *   iterate over every element \f$p_i\f$ of the power set of \f$S\F$
 *   except for\f$S\$ itself and \f$\emptyset\f$.  The resulting serial
 *   number is then the union of \f$p_i\f$ and -1 times the complement of
 *   \f$p_i\f$ in \f$S\f$.  This generates all the unique calculations
 *   for a VMFC(n) job.
 *
 *   To put the BSSE corrected properties back together requires more
 *   finesse.  What we want to do here is start with some serial number
 *   that is all real monomers.  Then we generate it's power set using
 *   this iterator.  Each resulting set is an interaction that contributes
 *   to the all real SN's interaction.  Within each of thes subinteractions
 *   we need to repeat the process except that now our SN is not all
 *   real.
 *
 *   To generalize this iterator to these other interactions it is just a
 *   matter of keeping any ghost monomer that enters into this iterator
 *   constant.
 *
 *
 *
 */
class VMFCnItr{
   private:
      SerialNumber SN_;
      SerialNumber InitialGhosts_;
      SerialNumber Value_;
     boost::shared_ptr<PowerSetItr<SerialNumber> > It_;
     void Next();
   public:
      ///Returns the current serial number
      const SerialNumber& operator*()const{return Value_;}
      ///Allows access to the current serial number's members
      const SerialNumber* operator->()const{return &Value_;}
      ///SN must contain no ghosts for the top-level iterator
      VMFCnItr(const SerialNumber& SN);
      ///Returns true if this iterator has completed
      bool Done()const{return It_->Done();}
      ///Increments this iterator
      const VMFCnItr& operator++(){Next();return *this;}

};

class VMFCn:public BSSEer{
   private:
      Molecule& Mol_;
   public:
      void CalcBSSE(NMers& Sys,unsigned int Stop,unsigned int Start)const;
      VMFCn(Molecule& Mol);
};

}}//End namespace


#endif /* SRC_LIB_LIBFRAG_LIBMOLECULE_UTILS_BSSEER_H_ */
