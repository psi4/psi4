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
#ifndef SRC_LIB_LIBFRAG_LIBMOLECULE_UTILS_ASSIGNATOMFUNCTOR_H_
#define SRC_LIB_LIBFRAG_LIBMOLECULE_UTILS_ASSIGNATOMFUNCTOR_H_
#include <boost/shared_ptr.hpp>
#include "../Atom.h"
#include "OrganicGeom.h"
namespace psi{
namespace LibMolecule{

///A simple functor to assign atoms to fxnalgroups based on Z
template <int Z, typename T, typename...Args>
class AssignAtomFunctor:public AssignAtomFunctor<Z+1,Args...>{
   public:
      void operator()(int i, const boost::shared_ptr<const Atom> AtomIn,
                             ConnGroups& Groups)const{
         if(AtomIn->Z()==Z)Groups[i]=boost::shared_ptr<const T>(new T(&i));
         else AssignAtomFunctor<Z+1,Args...>::operator()(i,AtomIn,Groups);
      }
};

///Termination condition to previous functor, prevents infinite recursion
template<int Z,typename T>
class AssignAtomFunctor<Z,T>{
   public:
      void operator()(int i, const boost::shared_ptr<const Atom> AtomIn,
            ConnGroups& Groups)const{
         if(AtomIn->Z()==Z)
            Groups[i]=boost::shared_ptr<T>(new T(&i));
         else throw PSIEXCEPTION("Your atom is beyond 118...");
      }
};

///An AssignAtomFunctor that knows about the entire periodic table
typedef AssignAtomFunctor<1,Hydrogen,Helium,Lithium,Beryllium,Boron,Carbon,
      Nitrogen,Oxygen,Fluorine,Neon,Sodium,Magnesium,Aluminium,Silicon,
      Phosphorus,Sulfur,Chlorine,Argon,Potassium,Calcium,Scandium,Titanium,
      Vanadium,Chromium,Manganese,Iron,Cobalt,Nickel,Copper,Zinc,Gallium,
      Germanium,Arsenic,Selenium,Bromine,Krypton,Rubidium,Strontium,Yttrium,
      Zirconium,Niobium,Molybdenum,Technetium,Ruthenium,Rhodium,Palladium,
      Silver,Cadmium,Indium,Tin,Antimony,Tellurium,Iodine,Xenon,Cesium,
      Barium,Lanthanum,Cerium,Praseodymium,Neodymium,Promethium,Samarium,
      Europium,Gadolinium,Terbium,Dysprosium,Holmium,Erbium,Thulium,
      Ytterbium,Lutetium,Hafnium,Tantalum,Tungsten,Rhenium,Osmium,Iridium,
      Platinum,Gold,Mercury,Thallium,Lead,Bismuth,Polonium,Astatine,Radon,
      Francium,Radium,Actinium,Thorium,Protactinium,Uranium,Neptunium,
      Plutonium,Americium,Curium,Berkelium,Californium,Einsteinium,Fermium,
      Mendelevium,Nobelium,Lawrencium,Rutherfordium,Dubnium,Seaborgium,
      Bohrium,Hassium,Meitnerium,Darmstadtium,Roentgenium,Copernicium,
      Ununtrium,Flerovium,Ununpentium,Livermorium,Ununseptium,Ununoctium> PeriodicTable;



}}//End namespaces



#endif /* SRC_LIB_LIBFRAG_LIBMOLECULE_UTILS_ASSIGNATOMFUNCTOR_H_ */
