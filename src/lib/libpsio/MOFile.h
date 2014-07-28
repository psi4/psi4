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
namespace psi{
class Matrix;
typedef boost::shared_ptr<Matrix> SharedMatrix;
class MOFile{
   private:
      int BasisNameLength;
      char* BasisName;
      int puream;
      double energy;
      int NIrrep;
      int* NSOPI;
      int* NAlpha;
      int* NBeta;
      SharedMatrix Ca;
      SharedMatrix Cb;
      void Copy(const MOFile& other);
   public:
      MOFile();
      MOFile(const MOFile& other){this->Copy(other);}
      void Init();
      ///Wrapper fxn that broadcasts this file
      void Broadcast(std::string& Comm,const int proc);
      ///Wrapper fxn that receives the broadcast of this file
      void Receive(std::string& Comm,const int proc);

      const MOFile& operator=(const MOFile& other){
         if(this!=&other)this->Copy(other);return *this;
      }
      ~MOFile();
      ///Returns an MOFile that is the direct sum of the current file and other
      MOFile DirectSum(const MOFile& other)const;
      void WriteFile();

};
}//End namespace psi


#endif /* MOFILE_H_ */
