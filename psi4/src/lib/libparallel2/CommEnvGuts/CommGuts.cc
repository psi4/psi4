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

#include "CommGuts.h"
#include "ParallelEnvironmentGuts.h"
#include "../Communicator.h"
#include "../../libPsiUtil/PsiRdmNumGen.h"
namespace psi{
namespace LibParallel{
typedef boost::shared_ptr<CommGuts> SharedThis;
typedef boost::shared_ptr<Comm::CommType_> SharedBase;

void CommGuts::FreeComm(){
   Active_=false;
   Env_->UpdateComms();
}

void CommGuts::Copy(const CommGuts& other){
   this->Env_=other.Env_;
   this->DaComm_=other.DaComm_;
   this->Name_=other.Name_;
   this->Active_=other.Active_;
}

CommGuts::CommGuts(const CommGuts& other){
   this->Copy(other);
}

const CommGuts& CommGuts::operator=(const CommGuts& other){
   if(this!=&other)this->Copy(other);
   return *this;
}

CommGuts::CommGuts(const std::string& Name,ParallelEnvironmentGuts* Env):
      Name_(Name),Env_(Env),DaComm_(new Comm::CommType_()),Active_(true){

}

void CommGuts::RegisterComm(boost::shared_ptr<Communicator> Comm)const{
   Env_->AddComm(Comm);
}

SharedThis CommGuts::MakeComm(const std::string& Name,const int Color)const{
   SharedThis temp(new CommGuts(*this));
   PsiRdmNumGen<> Rdm;
   std::stringstream tempName;
   tempName<<"COMM"<<Rdm();
   temp->Name_=(Name!=""?Name:tempName.str());
   temp->DaComm_=this->DaComm_->MakeComm(Color);
   return temp;
}
CommGuts::~CommGuts(){
   if(Active())FreeComm();
}
void CommGuts::Barrier()const{DaComm_->Barrier();}
int CommGuts::Probe(const int Sender,const int MessageTag,const bool Block)const{
   return DaComm_->Probe(Sender,MessageTag,Block);
}
int CommGuts::Me()const{
   return DaComm_->Me();
}
int CommGuts::NProc()const{
   return DaComm_->NProc();
}
}}//End namespaces