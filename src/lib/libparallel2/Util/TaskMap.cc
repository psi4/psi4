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

#include "TaskMap.h"
namespace psi{
namespace LibParallel{

int TaskMap::Rel(const int Abs)const{return Abs2Rel_[Abs];}
int TaskMap::Abs(const int Rel)const{return Rel2Abs_[Rel];}

TaskMap::TaskMap(const int NSlaves):TasksPerSlave_(NSlaves,0),
      TaskQueue_(NSlaves){}

int TaskMap::Slave(const int Abs)const{
   int Rel=Abs2Rel_[Abs];
   int slave=Task2Slave_[Rel];
   return slave;
}

int TaskMap::Index(const int Abs,const int slave)const{
   int slavei=slave;
   if(slave==-1)slavei=Slave(Abs);
   int size=TaskQueue_[slavei].size();
   int index=-1;
   for(int i=0;i<size&&index==-1;i++){
      if(TaskQueue_[slavei][i]==Abs)index=i;
   }
   if(index==-1)Error("This slave wasn't in charge of this task");
   return index;
}

void TaskMap::AddTask(const int Rel,const int Abs,const int Slave){
   Task2Slave_[Rel]=Slave;
   Rel2Abs_[Rel]=Abs;
   Abs2Rel_[Abs]=Rel;
   TasksPerSlave_[Slave]++;
   TaskQueue_[Slave].push_back(Abs);
}

int TaskMap::Offset(const int slave)const{
   int offset=0;
   for(int i=0;i<slave;i++)offset+=TasksPerSlave_[i];
   return offset;
}

}}
