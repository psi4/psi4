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
#ifndef MOFILE_H_
#define MOFILE_H_
#include <boost/shared_ptr.hpp>
#include "BinaryFile.h"
#include "psi4-dec.h"
namespace psi{
class Matrix;
typedef boost::shared_ptr<Matrix> SharedMatrix;
class MOFile:public BinaryFile{
   private:
      ///The alpha coefficients
      SharedMatrix Ca_;
      ///The beta coefficients
      SharedMatrix Cb_;
      void Copy(const MOFile& other);
   public:
      MOFile();
      MOFile(const MOFile& other):BinaryFile(other){this->Copy(other);}

      void FillFile(const int BasisNameLength,const char *BasisName,
                            const int PureAm,const double Energy,
                            const int NIrrep,const int NBf,
                            const int* NSOPI,const int* NAlphaI,
                            const int* NBetaI,SharedMatrix Ca,SharedMatrix Cb);
      void Read();

      ///Note Broadcast and Receive are written lazily and do not check if
      ///the destination is the same as source, so don't call them in
      ///non-MPI runs

      ///Wrapper fxn that broadcasts this file
      void Broadcast(std::string& Comm,const int proc);

      ///Wrapper fxn that receives the broadcast of this file
      void Receive(std::string& Comm,const int proc);

      const MOFile& operator=(const MOFile& other){
         BinaryFile::operator=(other);
         if(this!=&other)this->Copy(other);return *this;
      }

      ///Return the alpha or beta coefficient matrices
      SharedMatrix GetCa()const{return Ca_;}
      SharedMatrix GetCb()const{return Cb_;}

      ///Return the energy
      double GetEnergy()const;

      ///Returns the number of Irreps
      int GetNIrreps()const;

      ///Returns the number of Spin Orbitals per Irrep
      boost::shared_ptr<int[]> GetNSOPI()const;

      ///Returns the number of occupied alpha spin orbitals per Irrep
      boost::shared_ptr<int[]> GetNAlpha()const;

      ///Returns the number of occupied beta spin orbitals per Irrep
      boost::shared_ptr<int[]> GetNBeta()const;
      ~MOFile();
      void print_out();
      void Write();

};
}//End namespace psi


#endif /* MOFILE_H_ */
