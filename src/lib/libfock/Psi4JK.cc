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

#if HAVE_JK_FACTORY
#include "Psi4JK.h"
#include "Molecule.hpp"
#include "BasisSet.hpp"

//Psi dependencies
#include "psi4-dec.h"
#include "libmints/molecule.h"
#include "libmints/matrix.h"
#include "libmints/basisset.h"
#include "libmints/basisset_parser.h"
#include "libmints/wavefunction.h"
namespace psi {
   void Psi4JK::MakeMolecule(SharedPsiBasis &PsiBasis,pMol& JKMol) {
      boost::shared_ptr<psi::Molecule> PsiMol=PsiBasis->molecule();
      int natoms=PsiMol->natom();
      std::vector<double> Carts;
      std::vector<int> Zs;
      for(int i=0;i<natoms;i++) {
         Zs.push_back(PsiMol->Z(i));	//Capital Z is the atomic number
         Carts.push_back(PsiMol->x(i));
         Carts.push_back(PsiMol->y(i));
         Carts.push_back(PsiMol->z(i));
      }
      JKMol=pMol(new JKFactory::Molecule(Zs,Carts));
   }

   void Psi4JK::MakeBasis(SharedPsiBasis &PsiBasis,pMyBasis &JKBasis) {
      JKBasis=pMyBasis(new JKFactory::AOBasisSet(System));
      int NAtoms=System->NAtoms();
      for(int atom=0;atom<NAtoms;atom++) {
         boost::shared_ptr<JKFactory::AtomicBasisSet> ABS=(*JKBasis)[atom];
         int nshells=PsiBasis->nshell_on_center(atom);
         for(int shell=0;shell<nshells;shell++) {
            int l=(PsiBasis->shell(atom,shell)).am();
            bool isCart=(PsiBasis->shell(atom,shell)).is_cartesian();
            ABS->AddShell(l,isCart);
            int nprims=(PsiBasis->shell(atom,shell)).nprimitive();
            boost::shared_ptr<JKFactory::Shell> Dashell=(*ABS)[shell];
            for(int prim=0;prim<nprims;prim++) {
               double beta=(PsiBasis->shell(atom,shell)).exp(prim);
               double c0=(PsiBasis->shell(atom,shell)).erd_coef(prim);
               if(l>1) {
                  double prefactor=sqrt(pow(2.0,2*l)/df[2*l]);
                  c0/=prefactor;
               }
               Dashell->AddPrimitive(c0,beta);
            }
         }
      }
   }

   Psi4JK::Psi4JK(SharedPsiBasis &PsiBasis):
   JKFactory::Interface(1e-10,WorldComm->GetMPIComm()) {
      MakeMolecule(PsiBasis,System);
      MakeBasis(PsiBasis,Basis);
      Unitary=SharedMatrix(
            new Matrix(Basis->GetNBasis(),Basis->GetNBasis()));
      PsiReorderer Reorder;
      Reorder.BuildT(Basis);
      for(int i=0;i<Basis->GetNBasis();i++) {
         for(int j=0;j<Basis->GetNBasis();j++) {
            (*Unitary)(i,j)=((double *)Reorder)[i*Basis->GetNBasis()+j];
         }
      }
      IsIdentity=Reorder.IsIdentity;

   }

   void Psi4JK::UpdateDensity(SharedMatrix &Density,SharedMatrix &J,
         SharedMatrix &K) {
      JKFactory::Interface::UpdateJ(&((*J)(0,0)));
      JKFactory::Interface::UpdateK(&((*K)(0,0)));
      if(!IsIdentity) {
         SharedMatrix TempDensity(new Matrix(*Density));
         TempDensity->transform(Unitary);
         JKFactory::Interface::UpdateDensity(&(*TempDensity)(0,0));
      }
      else JKFactory::Interface::UpdateDensity(&(*Density)(0,0));
      if(!IsIdentity) {
         J->back_transform(Unitary);
         K->back_transform(Unitary);
      }
   }

   void Psi4JK::UpdateDensity(std::vector<SharedMatrix>& Density,
         std::vector<SharedMatrix>& J,std::vector<SharedMatrix>& K) {
      std::vector<double*> Rhos;
      std::vector<double*> Js;
      std::vector<double*> Ks;
      for(int i=0;i<Density.size();i++) {
         if(!IsIdentity) {
            SharedMatrix TempDensity=Density[i];
            TempDensity->transform(Unitary);
            Rhos.push_back(&(*TempDensity)(0,0));
         }
         else Rhos.push_back(&(*Density[i])(0,0));
         Js.push_back(&(*J[i])(0,0));
         if(K.size()>0)Ks.push_back(&(*K[i])(0,0));
      }
      JKFactory::Interface::UpdateJ(Js);
      JKFactory::Interface::UpdateK(Ks);
      JKFactory::Interface::UpdateDensity(Rhos);
      for(int i=0;i<Density.size()&&!IsIdentity;i++) {
         J[i]->back_transform(Unitary);
         if(K.size()>0)K[i]->back_transform(Unitary);
      }
   }

   void PsiReorderer::YourCartConverter(int l, int k, int& a, int& b,
         int &c) {
      int index=0;
      for(int i =0;i<=l && index>=0;i++) {
         for(int j=0;j<=i && index>=0;j++) {
            a = l - i;
            b = i - j;
            c = j;
            if(index==k)index=-1;
            else index++;
         }
      }
   }
   void PsiReorderer::YourPureConverter(int l, int k, int& m) {
      if(k==0)m=0;
      else m=(k%2==1 ? (k+1)/2: -k/2);
   }
}
#endif

