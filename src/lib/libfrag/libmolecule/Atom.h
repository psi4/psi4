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
#ifndef SRC_LIB_LIBFRAG_LIBMOLECULE_ATOM_H_
#define SRC_LIB_LIBFRAG_LIBMOLECULE_ATOM_H_

#include <vector>
#include "LibMoleculeBase.h"
namespace psi{
namespace LibMolecule{

///Allows easy editing of the type
typedef std::vector<double> Carts_t;

class Atom: protected LibMoleculeBase{
   private:
      Carts_t Carts_;
      double Mass_;
      double Charge_;
      int Z_;
      int NElec_;
      std::string Label_;
   public:
      ///Accessors for the Cartesian coordinates
      const Carts_t& Carts()const{return Carts_;}
      double operator[](const int i)const{return Carts_[i];}
      double& operator[](const int i){return Carts_[i];}

      ///Accessors for other properties
      double Mass()const{return Mass_;}
      double Charge()const{return Charge_;}
      int Z()const{return Z_;}
      int NElec()const{return NElec_;}
      std::string Label()const{return Label_;}

      ///What type of atom is this?
      bool IsPointCharge()const{return (NElec()==0&&Mass()==0.0);}
      bool IsGhost()const{return (IsPointCharge()&&Charge()==0.0);}
      bool IsDummy()const{return (IsGhost()&&Z_==0);}

      /** \brief Prints out the atom, useful for debugging
       *
       *   \param[in] Bohr If true coordinates are printed in Bohr
       *              (default =false)
       *   \param[in] DebugLevel How much are we printing (Default =1)
       *              levels are cumulative and defined as:
       *              Level 0: <Atomic_Symbol> x y z
       *              Level 1: units
       *              Level 2: Mass (m) and Charge (q)
       *              Level 3: NElectrons (NE) and Atomic Number (Z)
       *   \return The atom as a string
       */
      std::string PrintOut(const bool Bohr=false,const int DebugLevel=1)const;

      /** \brief The main constructor
       *
       *  It is advisable that you set these values through the appropriate
       *  specialization, i.e. make a GhostAtom if you want a ghost atom, etc.
       *  otherwise the checks like IsGhost may fail.  If you are just
       *  trying to make a normal atom: Carts, Z, and possibly IsBohr, are all
       *  you should have to set unless you are doing something really fancy.
       *
       *  \param[in] Carts The Cartesian coordinates of this atom
       *  \param[in] IsBohr If true Carts is in atomic units
       *  \param[in] Z The atomic number (negative for ghost atoms,
       *                0 for point charges, and dummy atoms)
       *  \param[in] Label The thing that gets printed out, specifying a Z
       *                other than 0 and Label="" will cause it to be looked
       *                up by Z
       *  \param[in] Charge Defaults to 0, only non-zero for point charges
       *  \param[in] NElec Defaults to 0, if it equals 0 and Z is >=1
       *                   will be set to Z
       *  \param[in] Mass Should be in Daltons, and defaults to 0.
       *                   If Z is >=1 and Mass==0.0 will be set appropriately
       *
       *
       */
      Atom(const Carts_t& Carts, const int Z=0,const bool IsBohr=true,
            const std::string Label="",const double Charge=0.0,
           const int NElec=0,const double Mass=0.0);
};

}}//End namespaces



#endif /* SRC_LIB_LIBFRAG_LIBMOLECULE_ATOM_H_ */
