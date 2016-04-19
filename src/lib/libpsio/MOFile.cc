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
#include "MOFile.h"
#include "libpsio/psio.h"
#include "libpsio/psio.hpp"
#include "psifiles.h"
#include "libmints/matrix.h"
#include "psi4-dec.h"
#
#include<cstring>

namespace psi {

double MOFile::GetEnergy()const{
   return GetDouble("SCF ENERGY")[0];
}

int MOFile::GetNIrreps()const{
   return GetInt("NIRREP")[0];
}

boost::shared_ptr<int[]> MOFile::GetNSOPI()const{
   return GetInt("NSOPI");
}

boost::shared_ptr<int[]> MOFile::GetNAlpha()const{
   return GetInt("NALPHAPI");
}

boost::shared_ptr<int[]> MOFile::GetNBeta()const{
   return GetInt("NBETAPI");
}


MOFile::MOFile() :BinaryFile(PSIF_SCF_MOS){
   AddVariable("BASIS NAME LENGTH",INT);
   AddVariable("SCF ENERGY",DOUBLE);
   AddVariable("PUREAM",INT);
   AddVariable("NIRREP",INT);
   AddVariable("NUMBER OF BASIS FUNCTIONS",INT);
   AddCoupledVariable("BASIS NAME",CHAR,"BASIS NAME LENGTH");
   AddCoupledVariable("NSOPI",INT,"NIRREP");
   AddCoupledVariable("NALPHAPI",INT,"NIRREP");
   AddCoupledVariable("NBETAPI",INT,"NIRREP");
}

void MOFile::FillFile(const int BasisNameLength,const char *BasisName,
                      const int PureAm,const double Energy,
                      const int NIrrep,const int NBf,
                      const int* NSOPI,const int* NAlphaI,
                      const int* NBetaI,SharedMatrix Ca,SharedMatrix Cb){
   Fill("BASIS NAME LENGTH",&BasisNameLength);
   Fill("PUREAM",&PureAm);
   Fill("SCF ENERGY",&Energy);
   Fill("NIRREP",&NIrrep);
   Fill("NUMBER OF BASIS FUNCTIONS",&NBf);
   Fill("BASIS NAME",BasisName);
   Fill("NSOPI",NSOPI);
   Fill("NALPHAPI",NAlphaI);
   Fill("NBETAPI",NBetaI);
   Ca_=Ca;
   Cb_=Cb;
}

void MOFile::Broadcast(std::string& CommIn, const int proc) {
   BinaryFile::Broadcast(CommIn,proc);
   double* buffer=&(*Ca_)(0,0);
   int nirrep=GetNIrreps();
   if(nirrep!=1)throw PSIEXCEPTION("Your matrix has too high of symmetry");
   int nsopi=GetNSOPI()[0];
   int nalpha=GetNAlpha()[0],nbeta=GetNBeta()[0];
   int length=nsopi*nalpha;
   boost::shared_ptr<const LibParallel::Communicator> Comm=
         WorldComm->GetComm();
   Comm->Bcast(buffer,length,proc);
   buffer=&(*Cb_)(0,0);
   length=nsopi*nbeta;
   Comm->Bcast(buffer,length,proc);
}

void MOFile::Receive(std::string& CommIn,const int proc){
   BinaryFile::Receive(CommIn,proc);
   int nirrep=GetNIrreps();
   if(nirrep!=1)throw PSIEXCEPTION("Your matrix has too high of symmetry");
   int nalpha=GetNAlpha()[0],nbeta=GetNBeta()[0];
   int nsopi=GetNSOPI()[0];
   Ca_=SharedMatrix(new Matrix("ALPHA MOS", nirrep,
         GetNSOPI().get(),GetNAlpha().get()));
   Cb_=SharedMatrix(new Matrix("BETA MOS", nirrep,
         GetNSOPI().get(),GetNBeta().get()));
   int length1=nsopi*nalpha,length2=nsopi*nbeta;
   double* abuffer=new double[length1];
   double* bbuffer=new double[length2];
   boost::shared_ptr<const LibParallel::Communicator> Comm=
         WorldComm->GetComm();
   Comm->Bcast(abuffer,length1,proc);
   Comm->Bcast(bbuffer,length2,proc);
   for(int i=0;i<nsopi;i++){
      for(int j=0;j<nalpha;j++)
         Ca_->set(i,j,abuffer[i*nalpha+j]);
      for(int j=0;j<nbeta;j++)
         Cb_->set(i,j,bbuffer[i*nbeta+j]);
   }
   delete [] abuffer;delete [] bbuffer;
}

MOFile::~MOFile() {}

void MOFile::Copy(const MOFile& other) {
   if(other.Ca_){
      this->Ca_=SharedMatrix(new Matrix((*other.Ca_)));
      this->Ca_->set_name(other.Ca_->name());
   }
   if(other.Cb_){
      this->Cb_=SharedMatrix(new Matrix((*other.Cb_)));
      this->Cb_->set_name(other.Cb_->name());
   }
}


void MOFile::Read() {
   BinaryFile::Read();
   int nirreps=GetNIrreps();
   Ca_=SharedMatrix(new Matrix("ALPHA MOS",nirreps,
       GetNSOPI().get(), GetNAlpha().get()));
   Ca_->load(psio_, FileNumber_, psi::Matrix::SubBlocks);
   Cb_=SharedMatrix(new psi::Matrix("BETA MOS",nirreps,
         GetNSOPI().get(), GetNBeta().get()));
   Cb_->load(psio_, FileNumber_, psi::Matrix::SubBlocks);
   psio_->close(FileNumber_, 1);
}



/*MOFile MOFile::DirectSum(const MOFile& other) const {
   if (this->NIrrep!=1&&other.NIrrep!=1)
      throw PSIEXCEPTION("I don't want to code direct sums up for symmetry");

   int NrowsCa=this->Ca->nrow()+other.Ca->nrow();
   int NrowsCb=this->Cb->nrow()+other.Cb->nrow();
   int NColsCa=this->Ca->ncol()+other.Ca->ncol();
   int NColsCb=this->Cb->ncol()+other.Cb->ncol();
   MOFile NewFile(*this);
   NewFile.energy+=other.energy;
   for (int i=0; i<NIrrep; i++) {
      NewFile.NAlpha[i]=this->NAlpha[i]+other.NAlpha[i];
      NewFile.NBeta[i]=this->NBeta[i]+other.NBeta[i];
      NewFile.NSOPI[i]=this->NSOPI[i]+other.NSOPI[i];
      NewFile.NBf=this->NBf+other.NBf;
   }
   NewFile.Ca->init(1, &NrowsCa, &NColsCa, "ALPHA MOS");
   NewFile.Cb->init(1, &NrowsCa, &NColsCb, "BETA MOS");
   NewFile.Ca->set(0.0);
   NewFile.Cb->set(0.0);
   for (int i=0; i<this->NSOPI[0]; i++) {
      for (int j=0; j<this->NAlpha[0]; j++)
         NewFile.Ca->set(i, j, (*this->Ca)(i, j));
      for (int j=0; j<this->NBeta[0]; j++)
         NewFile.Cb->set(i, j, (*this->Cb)(i, j));
   }
   int rowoff=this->NSOPI[0];
   int acoloff=this->NAlpha[0];
   int bcoloff=this->NBeta[0];
   for (int i=0; i<other.NSOPI[0]; i++) {
      for (int j=0; j<other.NAlpha[0]; j++)
         NewFile.Ca->set(i+rowoff, j+acoloff, (*(other.Ca))(i, j));

      for (int j=0; j<other.NBeta[0]; j++)
         NewFile.Cb->set(i+rowoff, j+bcoloff, (*(other.Cb))(i, j));

   }
   return NewFile;
}*/


void MOFile::print_out(){Ca_->print_out();}

void MOFile::Write(){
   BinaryFile::Write();
   //There's no reason why Write can't be const, but matrix's save isn't
   //Also matrix's copy, doesn't copy the name, which is annoying....
   Matrix CaCopy((*Ca_));
   CaCopy.set_name(Ca_->name());
   Matrix CbCopy((*Cb_));
   CbCopy.set_name(Cb_->name());
   CaCopy.save(psio_, FileNumber_, Matrix::SubBlocks);
   CbCopy.save(psio_, FileNumber_, Matrix::SubBlocks);
   psio_->close(FileNumber_, 1);
}
}      //End namespace psi