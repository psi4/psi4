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

#include "LibParallelHelper.h"
#include "exception.h"
#include "MPIJob.h"

namespace psi{
namespace LibParallel{

typedef boost::python::str PyStr;
typedef boost::python::list PyList;

LibParallelHelper::LibParallelHelper():Job_(NULL){}
LibParallelHelper::~LibParallelHelper(){if(Job_!=NULL)delete Job_;}

void LibParallelHelper::AddTask(PyStr& name,const int Prior){
   Tasks_.push_back(MPITask<PyStr>(name,Prior));
}

void LibParallelHelper::MakeJob(){
   if(Job_==NULL)Job_=new MPIJob<PyStr>(Tasks_);
   else throw PSIEXCEPTION("You already made a job...");
}

PyStr LibParallelHelper::Begin(){return Job_->Begin();}
PyStr LibParallelHelper::Next(){return Job_->Next();}
bool LibParallelHelper::Done(){return Job_->Done();}
PyList LibParallelHelper::Synch(PyList& Local,const int N){
   int size=len(boost::python::extract<PyList>(Local));
   std::vector<double> temp(size);
   for(int i=0;i<size;i++)
      temp[i]=boost::python::extract<double>(Local[i]);
   std::vector<double> results=Job_->Synch(temp,N);
   PyList ReturnList;
   for(int i=0;i<results.size();i++)ReturnList.append(results[i]);
   return ReturnList;

}
}}//End namesapces

