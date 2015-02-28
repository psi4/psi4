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
#ifndef BINARYFILE_H_
#define BINARYFILE_H_


#include <vector>
#include <string>
#include <boost/shared_ptr.hpp>
#include <boost/pointer_cast.hpp>
#include "libpsio/psio.h"
#include "libpsio/psio.hpp"
#include "exception.h"
#include "StringArray.h"
namespace psi {

///Types that are recognized
enum CTypes{INT,DOUBLE,CHAR};

/** \brief The purpose of this class is to form an abstract base class for
 *         the various binary files used in Psi.
 *
 *  When each binary file is readed and wroted certain data is expected to
 *  be in that binary file.  As of right now every time someone wants to use
 *  that file they need to take care to ensure that all of the data is inside
 *  the file for other parts of the program.  Consequentially, the act of
 *  opening a file means large chunks of code are being dumped all over the
 *  place.  In an ideal C++ world those large chunks of code would be
 *  mitigated to their own classes, so that if the structure of the file
 *  changes, it only changes in one place: the instance of this class.
 *
 *  This binary file class also has a nifty feature in that it is designed
 *  to save code.  Most binary files in Psi4 have a series of "Table of
 *  Contents" (ToC) entries associated with basic variables.  I have made it
 *  so that these entries will be read'n and written for you if you tell
 *  me about them at file creation time.  See the MOFile class for an
 *  example, but basically you need to use the AddVariable and Fill
 *  functions.  Call AddVariable in your constructor, and then Fill, in
 *  whatever fxn you use to add the actual values to your file.
 *
 *
 *
 *
 */
class BinaryFile {
   private:
      ///Makes this copy of other
      void Copy(const BinaryFile& other);

      /** \brief This a map of all the various members we are responsible
       *         for.
       *
       *  To understand how this works it helps to think of boost::any
       *  as a typesafe void* or a type that can be any C++ type.
       *  Basically it says we are mapping CTypes, to StringArrays,
       *  that will have different types, but we want them in the same
       *  array.
       */
      std::map<CTypes,boost::shared_ptr<PsiIOStringArray> > Members_;

   protected:
      ///The instance of psio we are using
      boost::shared_ptr<PSIO> psio_;

      ///The file we are working on
      int FileNumber_;

      ///The names of the variables the base class is responsible for
      std::vector<std::string> VariableNames_;

      ///The type associated with each variable in the base class
      std::map<std::string,CTypes> VariableTypes_;

      boost::shared_ptr<double[]> GetDouble(const std::string& Name)const;
      boost::shared_ptr<int[]> GetInt(const std::string& Name)const;
      boost::shared_ptr<char[]> GetChar(const std::string& Name)const;

      /** \brief Function designed to facilitate filling a binary file
       *
       * \param[in] Name The ASCII string to denote the variable in the file
       * \param[in] Length The size of the variable
       * \param[in] Values An array of the values to be copied into the file
       */
      template<typename T>
      void Fill(const std::string& Name, const T* Values);

      ///Declares a variable
      void AddVariable(const std::string& Name,const CTypes& Type);

      ///Declares a variable, whose length depends on "Coupled"'s value
      void AddCoupledVariable(const std::string& Name,
            const CTypes& Type,const std::string& Coupled);

   public:
      ///Sets pointer to psio and sets FileNumber appropriately
      BinaryFile(const int FileN=0);

      ///No memory to free
      virtual ~BinaryFile();

      ///Initializes this BinaryFile by copying other
      BinaryFile(const BinaryFile& other) {this->Copy(other);}

      ///Sets this BinaryFile equal to other
      const BinaryFile& operator=(const BinaryFile& other) {
         if (this!=&other) this->Copy(other);
         return *this;
      }

      ///Reads all member variables into class, doesn't close file
      virtual void Read();

      /** Writes all member variables from class, doesn't close file
       *  Should probably be const, but that won't mesh with matrix
       *  class
       */
      virtual void Write();

      ///Broadcasts member variables
      void Broadcast(const std::string& Comm, const int proc)const;

      ///Receives member variables
      void Receive(const std::string& Comm, const int proc);

};


template<typename T>
void BinaryFile::Fill(const std::string& Name,const T* values){
   boost::dynamic_pointer_cast<PsiIOSAImpl<T> >
      (Members_[VariableTypes_[Name]])->SetValue(Name,values);
}

}

#endif /* BINARYFILE_H_ */
