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

#include "BinaryFile.h"
#include "psi4-dec.h"

namespace psi {

void BinaryFile::Copy(const BinaryFile& other) {
   this->FileNumber_=other.FileNumber_;
   this->psio_=other.psio_;
   this->VariableNames_=other.VariableNames_;
   this->VariableTypes_=other.VariableTypes_;
   this->Members_=other.Members_;
}

BinaryFile::BinaryFile(const int FileN) :
      FileNumber_(FileN), psio_(psi::_default_psio_lib_) {
      Members_[INT]=boost::shared_ptr<PsiIOSAImpl<int> >
                        (new PsiIOSAImpl<int>());
      Members_[DOUBLE]=boost::shared_ptr<PsiIOSAImpl<double> >
                        (new PsiIOSAImpl<double>());
      Members_[CHAR]=boost::shared_ptr<PsiIOSAImpl<char> >
                         (new PsiIOSAImpl<char>());
}
typedef std::map<CTypes,boost::shared_ptr<PsiIOStringArray> > MemBase;

boost::shared_ptr<double[]> BinaryFile::GetDouble(const std::string& Name)const{
   return boost::dynamic_pointer_cast<PsiIOSAImpl<double> >
               (MapDeRef(Members_,DOUBLE))->GetValue(Name);
}

boost::shared_ptr<int[]> BinaryFile::GetInt(const std::string& Name)const{
   return boost::dynamic_pointer_cast<PsiIOSAImpl<int> >
           (MapDeRef(Members_,INT))->GetValue(Name);
}

boost::shared_ptr<char[]> BinaryFile::GetChar(const std::string& Name)const{
      return boost::dynamic_pointer_cast<PsiIOSAImpl<char> >
               (MapDeRef(Members_,CHAR))->GetValue(Name);
}

void BinaryFile::Read(){
   psio_->open(FileNumber_,PSIO_OPEN_OLD);
   for(int Var=0;Var<VariableNames_.size();Var++){
      std::string Name=VariableNames_[Var];
      Members_[VariableTypes_[Name]]->Read(psio_,Name,FileNumber_);
   }
}

void BinaryFile::Write(){
   psio_->open(FileNumber_,PSIO_OPEN_NEW);
   for(int Var=0;Var<VariableNames_.size();Var++){
      std::string Name=VariableNames_[Var];
      MapDeRef(Members_,MapDeRef(VariableTypes_,Name))->Write(psio_,Name,FileNumber_);
   }
}

BinaryFile::~BinaryFile(){}

void BinaryFile::Broadcast(const std::string& CommIn, const int proc)const{
   for(int i=0;i<VariableNames_.size();i++){
      std::string Name=VariableNames_[i];
      MapDeRef(Members_,MapDeRef(VariableTypes_,Name))->Broadcast(CommIn,proc,Name);
   }
}

void BinaryFile::Receive(const std::string& Comm, const int proc){
   for(int i=0;i<VariableNames_.size();i++){
      std::string Name=VariableNames_[i];
      Members_[VariableTypes_[Name]]->Receive(Comm,proc,Name);
   }
}


void BinaryFile::AddVariable(const std::string& Name,const CTypes& Type){
   VariableNames_.push_back(Name);
   VariableTypes_[Name]=Type;
   boost::shared_ptr<int[]> temp(new int[1]);
   temp[0]=1;
   Members_[Type]->SetLength(Name,temp);
   ///This is just a b.s. set, so the memory is allocated incase we need
   ///to couple it
   if(Type==INT)boost::dynamic_pointer_cast<PsiIOSAImpl<int> >(
         Members_[Type])->SetValue(Name,temp.get());
}

void BinaryFile::AddCoupledVariable(const std::string& Name,
      const CTypes& Type,const std::string& Coupled){
   VariableNames_.push_back(Name);
   VariableTypes_[Name]=Type;
   boost::shared_ptr<int[]> temp=GetInt(Coupled);
   Members_[Type]->SetLength(Name,temp);
}

} //End namespace psi
