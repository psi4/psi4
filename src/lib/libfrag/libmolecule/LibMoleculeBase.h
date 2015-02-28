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
#ifndef SRC_LIB_LIBFRAG_LIBMOLECULE_LIBMOLECULEBASE_H_
#define SRC_LIB_LIBFRAG_LIBMOLECULE_LIBMOLECULEBASE_H_

#include <iostream>
#include "exception.h"
#include "Units.h"
#include "AtomicData.h"
#include "qt.h"
namespace psi{
namespace LibMolecule{

class LibMoleculeBase{
   public:
      void Print(const std::string& Message)const{
         std::cout<<Message<<std::endl;
      }
      void Error(const std::string& Message)const{
         std::string Error="Error: ";
         Error+=Message;
         throw PSIEXCEPTION(Message);
      }
      double AngToBohr()const{
         BaseUnitConverter Conv;
         return Conv(ANGSTROM,BOHR);
      }
      std::string Sym(const int i)const{
         AtomicData CRC;
         return (0<i?CRC[i].AtSym():"GH");
      }
      double Mass(const int i)const{
         AtomicData CRC;
         return CRC[i].Mass();
      }
      /** \brief C=a*A*B+b*C
       *
       *  \param[in/out] C The result
       *  \param[in] A The matrix A
       *  \param[in] B The matrix B
       *  \param[in] RowA The rows of the resulting matrix
       *  \param[in] ConDim The contracting dimension,after any transposes
       *  \param[in] ColB The columns of the resulting matrix
       *  \param[in] TransA True if we are transposing A
       *  \param[in] TransB True if we are transposing B
       *  \param[in] a The value A.B is being scaled by
       *  \param[in] b The value C is being scaled by before being added
       *               to A.B
       */
      void DGEMM_Int(double* C,const double* A, const double* B,
                     int RowA,int ConDim,int ColB,
                     bool TransA=false,bool TransB=false,
                     double a=1.0,double b=0.0){
         char transA=(TransA?'T':'N');
         char transB=(TransB?'T':'N');
         C_DGEMM(transA,transB,
                 RowA,ColB,ConDim,
                 a,const_cast<double*>(A),(TransA?RowA:ConDim),
                 const_cast<double*>(B),(TransB?ConDim:ColB),
                 b,C,ColB);
      }
      ///Diagonalizes the matrix C, that is N by N
      void Diagonalize(double* C, const int N,double* EVals,bool EVec=true){
         int NewN=N;
         double NWork;
         C_DSYEV((EVec?'V':'N'),'U',NewN,C,N,EVals,&NWork, -1);
         double *Work=new double[(int)NWork];
         C_DSYEV((EVec?'V':'N'),'U',NewN,C,N,EVals,Work,NWork);
         delete [] Work;
      }
};

}}//End namespaces




#endif /* SRC_LIB_LIBFRAG_LIBMOLECULE_LIBMOLECULEBASE_H_ */
