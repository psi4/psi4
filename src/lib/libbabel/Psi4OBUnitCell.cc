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

#include "openbabel/math/spacegroup.h"
#include "openbabel/data.h"

#include "Psi4OBUnitCell.h"
namespace psi {

typedef OpenBabel::SpaceGroup OBSG;
typedef OpenBabel::transform3d OBTrans;
typedef OpenBabel::vector3 OBVec;

inline bool IsDone(const OBSG* sg,OBTrans const * value){
   OBTrans* aNULLptr=reinterpret_cast<OBTrans*>(NULL);
   return(value==aNULLptr);
}

std::list<OpenBabel::vector3> Transform2(const OBVec &v,const OBSG* sg) {
   static double prec=2e-5;
   std::list<OBVec> res;
   OpenBabel::transform3dIterator i;
   OBTrans const *value=sg->BeginTransform(i);
   bool done=IsDone(sg,value);
   while (!done) {
      OBVec t;
      t=(*value)*v;
      std::list<OBVec>::iterator j,jend=res.end();
      bool duplicate=false;
      for (j=res.begin(); j!=jend&&!duplicate; j++) {
         OBVec dt=t-(*j);
         if (fabs(dt.x())<prec && fabs(dt.y())<prec &&fabs(dt.z())<prec)
            duplicate=true;
      }
      if(!duplicate)res.push_back(t);
      value=sg->NextTransform(i);
      done=IsDone(sg,value);
   }
   return res;
}

inline bool areDuplicateAtoms(OBVec v1, OBVec v2) {
   double thresh=0.5;
   OBVec dr=v2-v1;
   for (int i=0; i<3; i++) {
      if (dr.x()<-1*thresh) dr.SetX(dr.x()+1.0);
      if (dr.x()>thresh) dr.SetX(dr.x()-1.0);
      if (dr.y()<-1*thresh) dr.SetY(dr.y()+1.0);
      if (dr.y()>thresh) dr.SetY(dr.y()-1.0);
      if (dr.z()<-1*thresh) dr.SetZ(dr.z()+1.0);
      if (dr.z()>thresh) dr.SetZ(dr.z()-1.0);
   }
   return (dr.length_2()<1e-6);
}

void OBUnitCellChild::FillUnitCell2(OpenBabel::OBMol *mol) {
   const OBSG *sg=GetSpaceGroup(); // the actual space group and transformations for this unit cell
   if(sg==NULL)return;

   // For each atom, we loop through: convert the coords back to
   //inverse space, apply the transformations and create new atoms
   OBVec uniqueV,newV,updatedCoordinate;

   // list of symmetry-defined copies of the atom
   std::list<OBVec>transformedVectors;

   std::list<OBVec>::iterator transformIterator,duplicateIterator;
   OpenBabel::OBAtom *newAtom;

   // keep the current list of unique atoms -- don't double-create
   std::list<OpenBabel::OBAtom*> atoms;
   // all coordinates to prevent duplicates
   std::list<OBVec>coordinates;

   bool foundDuplicate=false;
   FOR_ATOMS_OF_MOL(atom, *mol)
   atoms.push_back(&(*atom));

   std::list<OpenBabel::OBAtom*>::iterator i;
   for (i=atoms.begin(); i!=atoms.end(); ++i) {
      uniqueV=(*i)->GetVector();
      uniqueV=CartesianToFractional(uniqueV);
      coordinates.push_back(uniqueV);
      transformedVectors=Transform2(uniqueV, sg);
      int ii=0;
      for (transformIterator=transformedVectors.begin();
            transformIterator!=transformedVectors.end(); ++transformIterator) {
         ii++;
         updatedCoordinate=*transformIterator;
         foundDuplicate=false;

         // Check if the transformed coordinate is a duplicate of an atom
         for (duplicateIterator=coordinates.begin();
               duplicateIterator!=coordinates.end(); ++duplicateIterator) {
            if (areDuplicateAtoms(*duplicateIterator, updatedCoordinate)) {
               foundDuplicate=true;
               break;
            }
         }
         if (foundDuplicate) continue;

         coordinates.push_back(updatedCoordinate); // make sure to check the new atom for dupes
         newAtom=mol->NewAtom();
         newAtom->Duplicate(*i);
         newAtom->SetVector(FractionalToCartesian(updatedCoordinate));
      } // end loop of transformed atoms
      (*i)->SetVector(FractionalToCartesian(uniqueV)); // move the atom back into the unit cell
   } // end loop of atoms
   SetSpaceGroup(1); // We've now applied the symmetry, so we should act like a P1 unit cell
}
}//End namespace

