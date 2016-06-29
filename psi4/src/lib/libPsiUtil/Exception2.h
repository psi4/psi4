/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
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
 * @END LICENSE
 */
#ifndef SRC_LIB_LIBPSIUTIL_EXCEPTION2_H_
#define SRC_LIB_LIBPSIUTIL_EXCEPTION2_H_

#include<sstream>
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
 *  <Current Fxn>
 *  <Fxn That Called Current Fxn>
 *  <Fxn That Called The Fxn That Called Current Fxn>
 *  \endverbatim
 *
 *  That is your exception will automatically add the current file, line
 *  number, and three (now 5) most recent function calls to the error message.
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
      ///Our actual error message
      std::string Error_;
      ///The function called when this exception is tipped
      virtual const char* what() const throw(){return Error_.c_str();}
   public:
      /** \brief Makes our exception
       *
       *   The usage of this class should be pretty straightforward
       *   from the class documentation so here I just specify what the
       *   arguments are:
       *
       *  \param[in] arg The message we want printed
       *  \param[in] file A string containing the full path to the source
       *                  file.  Added automatically by the macro.
       *  \param[in] line The line number our error occurred on.  Added
       *                  automatically by the macro.
       *  \param[in] NCalls The number of functions to print in the call
       *                    back tree.  Defaults to 5.
       *
       */
      PsiException2(const std::string& arg,
                    const char* file,
                    const int line,
                    const size_t NCalls=5);
      ///Compiler can now STFU
      virtual ~PsiException2()throw(){}
};

///The macro mentioned in the PsiException2 documentation
#define PSIERROR(arg) throw PsiException2(arg,__FILE__,__LINE__);

}



#endif /* SRC_LIB_LIBPSIUTIL_EXCEPTION2_H_ */