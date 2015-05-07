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

#include "../SuperCell.h"
#include "../Utils/GeomManipulator.h"
namespace psi{
namespace LibMolecule{

void SuperCell::FormSuperCell(const int dims[],const double cellvecs[],
      const UnitCell& UC,int depth){
   if(depth<3){
      //For a dim of d, we have d positive and d negative, plus no move
      //So 2d+1 moves, we'll make the even ones the negative shift
      for (int i=0; i<=2*dims[depth]; i++){
         double NewCellVecs[3];
         for(int j=0;j<3;j++)NewCellVecs[j]=cellvecs[j];
         NewCellVecs[depth]*=(double)(i%2==0?-i/2:(i+1)/2);
         FormSuperCell(dims,NewCellVecs,UC,depth+1);
      }
   }
   else{//End recursion
      double TVec[3];
      DGEMM_Int(&TVec[0],Frac2Cart_,&cellvecs[0],3,3,1);
      Molecule ShiftCell(UC);
      GeomManipulator Manip(&ShiftCell);
      Manip.Translate(TVec);
      Manip.Set();
      (*this)+=ShiftCell;
   }
}

SuperCell::SuperCell(const UnitCell& OrigCell,
      const std::vector<int>& NShells):
      OrigCell_(new UnitCell(OrigCell)){
   int NTimes[3];
   double cellvec[3];
   for(int i=0;i<3;i++){
      NTimes[i]=NShells[(NShells.size()==1?0:i)];
      cellvec[i]=1.0;
      sides_[i]=OrigCell.Sides()[i];
      angles_[i]=OrigCell.Angles()[i];
   }
   SetTrans();
   FormSuperCell(NTimes,cellvec,OrigCell);
   //Make the sides consistent with the new size
   for(int i=0;i<3;i++)sides_[i]*=(NTimes[i]+1);

}

}}//End namespaces

