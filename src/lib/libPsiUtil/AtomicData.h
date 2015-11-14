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
#ifndef SRC_LIB_LIBPSIUTIL_ATOMICDATA_H_
#define SRC_LIB_LIBPSIUTIL_ATOMICDATA_H_

#include <vector>
#include <boost/shared_ptr.hpp>
#include "PsiMap.h"
namespace psi{

/** The objects in this header file are used like this:
 *
 *  1)Declare an AtomicData object:
 *  \code
 *  AtomicData CRC;
 *  \endcode
 *  2)Pick the element you want data for by Z not Z-1:
 *  \code
 *  const AtomData& HeData=CRC[2];//Not Li's data!!!
 *  \endcode
 *  3)Get the data pieces you want:
 *  \code
 *  double HeMass=HeData.Mass();
 *  \endcode
 *
 *  Note that steps 2 and 3 can be combined:
 *  \code
 *  double HeMass=CRC[2].Mass();
 *  \endcode
 *
 */


///Mini class to hold data for various isotopes
class Isotope{
      double Mass_;
      std::string Label_;
   public:
      double Mass()const {return Mass_;}
      std::string AtSym()const{return Label_;}
      Isotope(const double Mass,const std::string& Label):
         Mass_(Mass),Label_(Label){}
};

class AtomicData;
/** \brief A class to hold basic physical info about an atom
 *
 *  In order to allow for isotope data we give each atom an array
 *  of isotopes that holds the default symbol and mass as the 0-th element
 */
class AtomData{
   private:
      void AddIsotope(const std::string& Label,const double Mass);
      //Atomic Data does all the setup for this class at construction
      friend class AtomicData;
      std::vector<Isotope> Isotopes_;
      PsiMap<std::string,int> LookUp_;

      std::string FullName_;

      ///Covalent Radius (a.u.)
      double CovRad_;

      ///van Der Waal radius (a.u.)
      double VDWRad_;

   public:
      ///The number of unique isotopes recognized by Psi4
      int NIsotopes()const{return Isotopes_.size()-1;}

      ///Returns the i-th isotope's mass (0 is default isotope)
      double Mass(const std::string& Label)const{return Mass(LookUp_[Label]);}
      double Mass(const int i=0)const{return Isotopes_[i].Mass();}

      ///Returns the i-th isotope's AtSym (0 is default isotope, with no mass)
      std::string AtSym(const int i=0)const{return Isotopes_[i].AtSym();}
      ///Full element name
      std::string Name()const{return FullName_;}
      double CovRad()const{return CovRad_;}
      double VDWRad()const{return VDWRad_;}

      AtomData(const std::string& Label,const std::string& Full,
            const double Mass,const double CovRad,const double VDWRad):
               FullName_(Full),CovRad_(CovRad),VDWRad_(VDWRad){
         Isotopes_.push_back(Isotope(Mass,Label));
      }
};

class AtomicData{
      boost::shared_ptr<std::vector<AtomData> > Data_;
   public:
      const AtomData& operator[](const int i)const{
         return (*Data_)[i];
      }
      AtomicData();
};

}//End namespace


#endif /* SRC_LIB_LIBPSIUTIL_ATOMICDATA_H_ */
