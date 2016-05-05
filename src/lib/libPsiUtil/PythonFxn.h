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
#ifndef SRC_LIB_LIBPSIUTIL_PYTHONFXN_H_
#define SRC_LIB_LIBPSIUTIL_PYTHONFXN_H_
#include <vector>
#include <iostream>
#include <boost/python.hpp>
#include <boost/python/object.hpp>
namespace psi {
/**** The objects in this file are designed to be used as wrappers
 * to python functions.  You probably just want to read this block
 * and not look at the code below (it's nasty template stuff)
 *
 * The syntax is:
 * PythonFxn<T> MyPythonFxn("ModuleItIsIn","NameOfFuntion",
 *     "WhatTheArgumentListLooksLike");
 *
 * Where T is the return type of the function.  For functions with no
 * return omit the template parameter, i.e.
 * PythonFxn<> MyPythonFxn....
 *
 * The last argument is a printf-like syntax string, for full details
 * consult boost's documentation.  The general gist of it is, the c++
 * type is given in ():
 * - s is a string (char *)
 * - u a null-terminated buffer of unicode data (Py_UNICODE *)
 * - i an integer (int)
 * - l long integer (long int)
 * - c a character (char)
 * - d a double (double)
 * - f a float (float)
 * - ( x y z ) makes a tuple with values of types x, y, and z
 * - [ x y z ] makes a list  with values of types x, y, and z
 * - { x y z a } makes a dictionary with values of types x, y, z, and a
 *           x is the key, y is the value of entry one, and z is the key
 *           and a is the value of entry two
 *
 * Note: If you really want to, you can make the syntax more python-like
 * for tuples, lists, and dictionaries, e.g.
 * (x,y,z), [x,y,z], and {x:y,z:a}.  The result will be the same as above.
 *
 * Some examples to cement the idea (note you need to pass the quotes in
 * these examples).  If your function takes an integer
 * and a double your declaration would be:
 *
 * PythonFxn<T> MyPythonFxn("ModuleItIsIn","NameOfFuntion","i d")
 *
 * a dictionary of two doubles, with integer keys:
 *
 * PythonFxn<T> MyPythonFxn("ModuleItIsIn","NameOfFuntion","{i d i d}")
 *
 * Finally, if it takes no argument:
 *
 * PythonFxn<T> MyPythonFxn("ModuleItIsIn","NameOfFuntion","")
 *
 * To call the function the syntax is:
 * T ReturnValue=MyPythonFxn(Parm1,Parm2,Parm3,etc...)
 *
 * or (for the void case):
 *
 * MyPythonFxn(Param1,Param2,Param3,etc....)
 *
 *
 */

///The basic implementation shared between the void, and non-void fxn
class PythonFxnGuts {
   protected:
      ///The Module
      PyObject* Module_;
      ///The function
      PyObject* Fxn_;
      ///The arguments we are passing to the function
      PyObject* Args_;
      ///The syntax the user told us it has
      std::string Syntax_;
      ///A vector of the parameter types
      std::vector<char> ParamTypes_;
   public:
      ///Ensures our memory is taken care of correctly
      ~PythonFxnGuts();
      ///Makes a new function with syntax given in the c string
      PythonFxnGuts(const char* Module, const char* FxnName,
                    const char* Syntax=NULL);
};
///Stupid class to flag that our function "returns" void
class VoidReturn{};

///Fxn that crashes our program if we couldn't bind a python object
void Try(PyObject* &Obj,PyObject* Command);

static std::string HandleError(){
   using namespace boost::python;
      using namespace boost;

      PyObject *exc,*val,*tb;
      object formatted_list, formatted;
      PyErr_Fetch(&exc,&val,&tb);
      handle<> hexc(exc),hval(allow_null(val)),htb(allow_null(tb));
      object traceback(import("traceback"));
      if (!tb) {
          object format_exception_only(traceback.attr("format_exception_only"));
          formatted_list = format_exception_only(hexc,hval);
      } else {
          object format_exception(traceback.attr("format_exception"));
          formatted_list = format_exception(hexc,hval,htb);
      }
      formatted = boost::python::str("\n").join(formatted_list);
      return extract<std::string>(formatted);
}

///Actual defintion of our fxn (if it returns a value)
template <typename T=VoidReturn>
class PythonFxn:public PythonFxnGuts{
   public:
      ///Wrapper to the base class's constructor
      PythonFxn(const char* Module,const char* FxnName,
                const char* Syntax):PythonFxnGuts(Module,FxnName,Syntax){}
      ///Calls our function and returns a value.
      template<typename...Args>
      T operator()(Args...args){
         Try(Args_, Py_BuildValue(Syntax_.c_str(), args...));

           PyObject* ret;
           try{ret=PyEval_CallObject(Fxn_, Args_);}
           catch(boost::python::error_already_set &){
              if(PyErr_Occurred())std::cout<<HandleError();
           }
           return T(boost::python::handle<>(ret));
      }
};

///Our actual fxn if we don't return a value
template <>
class PythonFxn<VoidReturn>:public PythonFxnGuts{
   public:
      ///Wrapper to base class's constructor
      PythonFxn(const char* Module,const char* FxnName,
                const char* Syntax):PythonFxnGuts(Module,FxnName,Syntax){}
      ///Calls our function and returns a value.  Can have up to 10 params
      template<typename...Args>
      void operator()(Args...args){
         Try(Args_, Py_BuildValue(Syntax_.c_str(), args...));
         try{ PyEval_CallObject(Fxn_, Args_);}
         catch(boost::python::error_already_set &){
            if(PyErr_Occurred())std::cerr<<HandleError();
         }
      }
};

} //End namespace
#endif /* SRC_LIB_LIBPSIUTIL_PYTHONFXN_H_ */