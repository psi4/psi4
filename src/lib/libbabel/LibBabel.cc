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

#include "LibBabel.h"
#include "openbabel/obconversion.h"
#include "openbabel/obiter.h"
#include "openbabel/mol.h"
#include "openbabel/shared_ptr.h"
#include "openbabel/base.h"
#include "Psi4OBUnitCell.h"
#include "masses.h"
#include <iostream>
#include <sstream>
namespace psi {
namespace LibBabel {

std::string OBMol2String(OpenBabel::OBMol& mol) {
   std::stringstream Molecule;
      FOR_ATOMS_OF_MOL(atom, mol){
      Molecule<<atomic_labels[atom->GetAtomicNum()]<<" "<<atom->GetX()<<" "
            <<atom->GetY()<<" "<<atom->GetZ()<<"\n";
      }
   return Molecule.str();
}

void MakeOBMolFromFile(const std::string& filename, OpenBabel::OBMol& mol) {
   // Open the file.
   std::ifstream ifs(filename.c_str());
   if (!ifs) {
      std::cout<<"Could not open "<<filename<<" for reading."<<std::endl;
   }
   // Create the OBConversion object.
   OpenBabel::OBConversion conv;
   OpenBabel::OBFormat *Format=conv.FormatFromExt(filename.c_str());
   if (!conv.SetInFormat(Format)) {
      std::cout<<"Could not find input format for file "<<filename<<std::endl;
   }
   // Read the molecule.
   if (!conv.Read(&mol, &ifs)) {
      std::cout<<"Could not read molecule from file "<<filename<<std::endl;
   }
}

void FormSuperCell(const std::vector<int>& dims,
      const std::vector<OpenBabel::vector3>& cellvecs,
      const OpenBabel::OBMol& UnitCell,int depth,
      OpenBabel::vector3& TVec,OpenBabel::OBMol& mol){
   if(depth<3){
      for (int i=0; i<=dims[depth]; i++){
            OpenBabel::vector3 vectori=cellvecs[depth];
            vectori*=((double)(i));
            TVec+=vectori;
            FormSuperCell(dims,cellvecs,UnitCell,depth+1,TVec,mol);
            vectori*=-1.0;
            TVec+=vectori;
      }
   }
   else{//End recursion
      OpenBabel::OBMol newmol(UnitCell);
      newmol.Translate(TVec);
      mol+=newmol;
   }
}

boost::python::str ParseFile(const std::string& filename, const int dim) {

   OpenBabel::OBMol mol;
   MakeOBMolFromFile(filename, mol);

   OpenBabel::OBMol unitcell(mol);
   if (mol.HasData(OpenBabel::OBGenericDataType::UnitCell)) {
      OpenBabel::OBUnitCell* uc=
            dynamic_cast<OpenBabel::OBUnitCell*>(mol.GetData(
                  OpenBabel::OBGenericDataType::UnitCell));
      OBUnitCellChild ouruc(*uc);
      ouruc.FillUnitCell2(&unitcell);
      std::vector<int> dims(3, dim);
      std::vector<OpenBabel::vector3>cellvecs=uc->GetCellVectors();
      OpenBabel::vector3 TVec;
      //assume dims are x,y,z
      FormSuperCell(dims,cellvecs,unitcell,0,TVec,mol);
   }
   std::string Molecule=OBMol2String(mol);
   boost::python::str Mol2Return(Molecule);
   return Mol2Return;
}

}
} //End namespaces

