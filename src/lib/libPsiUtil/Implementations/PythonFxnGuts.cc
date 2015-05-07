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
#include "PythonFxn.h"
#include <iostream>
namespace psi{

void Try(PyObject*& Obj,PyObject* Command) {
   Obj=Command;
   if(Obj==NULL){
      PyErr_Print();
      exit(1);
   }
}
PythonFxnGuts::~PythonFxnGuts() {
   if (Args_!=NULL) Py_DECREF(Args_);
   if (Fxn_!=NULL) Py_DECREF(Fxn_);
   if (Module_!=NULL) Py_DECREF(Module_);
}

PythonFxnGuts::PythonFxnGuts(const char* Module, const char* Fxn,
      const char* Syntax) :
      Syntax_(""),Module_(NULL), Fxn_(NULL), Args_(NULL) {
   Syntax_+="(";
   Syntax_+=Syntax;
   Syntax_+=")";
   Try(Module_, PyImport_ImportModule(Module));
   Try(Fxn_, PyObject_GetAttrString(Module_, Fxn));
   int i=0;
   if (Syntax!=NULL) {
      while (Syntax[i]!='\0') {
         if (Syntax[i]=='s'||Syntax[i]=='z'||Syntax[i]=='u'||Syntax[i]=='i'
               ||Syntax[i]=='b'||Syntax[i]=='h'||Syntax[i]=='l'||Syntax[i]=='c'
               ||Syntax[i]=='d'||Syntax[i]=='f'||Syntax[i]=='D'||Syntax[i]=='O'
               ||Syntax[i]=='S'||Syntax[i]=='U'||Syntax[i]=='N')
            ParamTypes_.push_back(Syntax[i]);
         i++;
      }
   }
}


}//End namespace


