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

#include <boost/any.hpp>
#include <boost/shared_ptr.hpp>
#include <map>
#include <vector>
#include <string>
#include "libpsio/psio.h"
#include "libpsio/psio.hpp"
namespace psi {

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
 */
class BinaryFile {
   private:
      ///Makes this copy of other
      void Copy(const BinaryFile& other);

   protected:
      ///The instance of psio we are using
      boost::shared_ptr<PSIO> psio;

      ///The file we are working on
      int FileNumber;
   public:
      ///Sets pointer to psio and sets FileNumber appropriately
      BinaryFile(const int FileN=0);

      ///No memory to free
      virtual ~BinaryFile() {}

      ///Initializes this BinaryFile by copying other
      BinaryFile(const BinaryFile& other) {
         this->Copy(other);
      }

      ///Sets this BinaryFile equal to other
      const BinaryFile& operator=(const BinaryFile& other) {
         if (this!=&other) this->Copy(other);
         return *this;
      }

      virtual void Read()=0;
      virtual void Write()=0;
};
}

#endif /* BINARYFILE_H_ */
