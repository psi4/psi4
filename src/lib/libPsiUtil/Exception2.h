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
#ifndef SRC_LIB_LIBPSIUTIL_EXCEPTION2_H_
#define SRC_LIB_LIBPSIUTIL_EXCEPTION2_H_

#include<sstream>
#include<vector>
#include<execinfo.h>
#include<exception>

namespace psi{


/** \brief A basic exception class that tells you helpful info about the
 *         error.
 *
 *  Much like the other classes in this library, this class
 *  reinvents another class.  This time around I've taken it upon myself
 *  to change how errors are thrown.  Now if you use this error class,
 *  via the macro defined below it, you will get an output like:
 *
 *  \verbatim
 *  Fatal Error: <Error Message>
 *  Error occurred in file: <File Error Occurred In> on line: <da line number>
 *  The most recent 3 function calls were:
 *  ../objdir/bin/psi4(<Current Fxn>) [0x1f16ae5]
 *  ../objdir/bin/psi4(<Fxn That Called Current Fxn>) [0x1f15a1c]
 *  ../objdir/bin/psi4(<Fxn That Called The Fxn That Called Current Fxn>) [0x1f15fc2]
 *  \endverbatim
 *
 *  That is your exception will automatically add the current file, line
 *  number, and three most recent function calls to the error message.
 *
 *  The code you now use looks like:
 *  \code
 *  PSIERROR(<Error Message>)
 *  //Not this:
 *  //throw PSIEXCEPTION(<Error Message>)
 *  \endcode
 *
 *  An error is still thrown (of type std::exception) so feel free to
 *  catch it.
 *
 */
class PsiException2: public std::exception{
   private:
      std::string Error_;
      virtual const char* what() const throw(){
         return Error_.c_str();
      }
   public:
      PsiException2(const std::string& arg,
                    const char* file,
                    const int line,
                    const size_t NCalls=5){
         std::stringstream Error;
         Error<<std::endl<<"Fatal Error: "<<arg<<std::endl
             <<"Error occurred in file: "<<file<<" on line: "<<line<<std::endl;
         std::vector<void *> Stack(NCalls);
         char **strings;
         int size=backtrace(&Stack[0],NCalls);
         Error<<"The most recent "<<(size<NCalls?size:NCalls)<<" function calls were:"<<std::endl;
         strings=backtrace_symbols(&Stack[0],size);
         for(int i=0;i<size;i++)
            Error<<strings[i]<<std::endl;
         Error_=Error.str();
      }
};


#define PSIERROR(arg) throw PsiException2(arg,__FILE__,__LINE__);

}



#endif /* SRC_LIB_LIBPSIUTIL_EXCEPTION2_H_ */
