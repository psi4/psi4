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
#include "MOFile.h"
#include "libpsio/psio.h"
#include "libpsio/psio.hpp"
#include "psifiles.h"
#include "libmints/matrix.h"
#include "psi4-dec.h"
#include<cstring>

namespace psi {

void MOFile::Broadcast(std::string& Comm, const int proc) {
   WorldComm->bcast(&BasisNameLength, 1, proc, Comm);
   WorldComm->bcast(BasisName,BasisNameLength,proc,Comm);
   WorldComm->bcast(&puream,1,proc,Comm);
   WorldComm->bcast(&energy,1,proc,Comm);
   WorldComm->bcast(&NIrrep,1,proc,Comm);
   WorldComm->bcast(NSOPI,NIrrep,proc,Comm);
   WorldComm->bcast(NAlpha,NIrrep,proc,Comm);
   WorldComm->bcast(NBeta,NIrrep,proc,Comm);
   double* buffer=&(*Ca)(0,0);
   WorldComm->bcast(buffer,NSOPI[0]*NAlpha[0],proc,Comm);
   buffer=&(*Cb)(0,0);
   WorldComm->bcast(buffer,NSOPI[0]*NBeta[0],proc,Comm);
}

void MOFile::Receive(std::string& Comm,const int proc){
   WorldComm->bcast(&BasisNameLength, 1, proc, Comm);
   BasisName=new char[BasisNameLength];
   WorldComm->bcast(BasisName,BasisNameLength,proc,Comm);
   WorldComm->bcast(&puream,1,proc,Comm);
   WorldComm->bcast(&energy,1,proc,Comm);
   WorldComm->bcast(&NIrrep,1,proc,Comm);
   NSOPI=new int[NIrrep];
   NAlpha=new int[NIrrep];
   NBeta=new int[NIrrep];
   WorldComm->bcast(NSOPI,NIrrep,proc,Comm);
   WorldComm->bcast(NAlpha,NIrrep,proc,Comm);
   WorldComm->bcast(NBeta,NIrrep,proc,Comm);
   double* abuffer=new double[NSOPI[0]*NAlpha[0]];
   double* bbuffer=new double[NSOPI[0]*NBeta[0]];
   WorldComm->bcast(abuffer,NSOPI[0]*NAlpha[0],proc,Comm);
   WorldComm->bcast(bbuffer,NSOPI[0]*NBeta[0],proc,Comm);
   Ca=SharedMatrix(new Matrix("ALPHA MOS", NIrrep, NSOPI, NAlpha));
   Cb=SharedMatrix(new Matrix("BETA MOS", NIrrep, NSOPI, NAlpha));
   for(int i=0;i<NSOPI[0];i++){
      for(int j=0;j<NAlpha[0];j++){
         Ca->set(i,j,abuffer[i*NAlpha[0]+j]);
      }
      for(int j=0;j<NBeta[0];j++){
         Cb->set(i,j,bbuffer[i*NBeta[0]+j]);
      }
   }
   delete [] abuffer;delete [] bbuffer;
}
MOFile::~MOFile() {
   if (BasisName!=NULL) delete[] BasisName;
   if (NSOPI!=NULL) delete[] NSOPI;
   if (NAlpha!=NULL) delete[] NAlpha;
   if (NBeta!=NULL) delete[] NBeta;
   BasisName=NULL;
   NSOPI=NAlpha=NBeta=NULL;
}

void MOFile::Copy(const MOFile& other) {
   this->BasisNameLength=other.BasisNameLength;
   this->BasisName=new char[BasisNameLength];
   std::memcpy(this->BasisName, other.BasisName,
         this->BasisNameLength*sizeof(char));
   this->puream=other.puream;
   this->energy=other.energy;
   this->NIrrep=other.NIrrep;
   this->NSOPI=new int[this->NIrrep];
   this->NAlpha=new int[this->NIrrep];
   this->NBeta=new int[this->NIrrep];
   std::memcpy(this->NSOPI, other.NSOPI, this->NIrrep*sizeof(int));
   std::memcpy(this->NAlpha, other.NAlpha, this->NIrrep*sizeof(int));
   std::memcpy(this->NBeta, other.NBeta, this->NIrrep*sizeof(int));
   if(other.Ca){
      this->Ca=SharedMatrix(new Matrix((*other.Ca)));
      this->Ca->set_name(other.Ca->name());
   }
   if(other.Cb){
      this->Cb=SharedMatrix(new Matrix((*other.Cb)));
      this->Cb->set_name(other.Cb->name());
   }
}

void MOFile::Read() {
   psio->open(FileNumber, PSIO_OPEN_OLD);
   psio->read_entry(FileNumber,"BASIS NAME LENGTH",(char *)&BasisNameLength,
         sizeof(int));
   BasisName=new char[BasisNameLength];
   psio->read_entry(FileNumber, "BASIS NAME", BasisName,
         BasisNameLength*sizeof(char));
   psio->read_entry(FileNumber, "PUREAM", (char *)(&puream), sizeof(int));
   psio->read_entry(FileNumber, "SCF ENERGY", (char *)&(energy),
         sizeof(double));
   psio->read_entry(FileNumber, "NIRREP", (char *)&(NIrrep), sizeof(int));

   std::size_t size=NIrrep*sizeof(int);
   NSOPI=new int[size];
   NAlpha=new int[size];
   NBeta=new int[size];
   psio->read_entry(FileNumber, "NSOPI", (char *)NSOPI, size);
   psio->read_entry(FileNumber, "NALPHAPI", (char *)NAlpha, size);
   psio->read_entry(FileNumber, "NBETAPI", (char *)NBeta, size);
   Ca=SharedMatrix(new psi::Matrix("ALPHA MOS", NIrrep, NSOPI, NAlpha));
   Ca->load(psio, FileNumber, psi::Matrix::SubBlocks);
   Cb=SharedMatrix(new psi::Matrix("BETA MOS", NIrrep, NSOPI, NBeta));
   Cb->load(psio, FileNumber, psi::Matrix::SubBlocks);
   psio->close(FileNumber, 1);
}

MOFile::MOFile() :
      BasisName(NULL), NSOPI(NULL), NAlpha(NULL), NBeta(NULL),
      BasisNameLength(0),puream(-1),energy(0.0),NIrrep(0),
      BinaryFile(PSIF_SCF_MOS){
}

MOFile MOFile::DirectSum(const MOFile& other) const {
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
}
void MOFile::print_out(){Ca->print_out();}
void MOFile::Write() {
   psio->open(FileNumber, PSIO_OPEN_NEW);
   psio->write_entry(FileNumber, "SCF ENERGY", (char *)&(energy),
         sizeof(double));
   psio->write_entry(FileNumber, "NIRREP", (char *)&(NIrrep), sizeof(int));
   std::size_t size=NIrrep*sizeof(int);
   psio->write_entry(FileNumber, "NSOPI", (char *)NSOPI, size);
   psio->write_entry(FileNumber, "NALPHAPI", (char *)NAlpha, size);
   psio->write_entry(FileNumber, "NBETAPI", (char *)NBeta, size);
   psio->write_entry(FileNumber, "BASIS NAME LENGTH",
         (char *)(&BasisNameLength), sizeof(int));
   psio->write_entry(FileNumber, "BASIS NAME", BasisName,
         BasisNameLength*sizeof(char));
   psio->write_entry(FileNumber, "PUREAM", (char *)(&puream), sizeof(int));
   Ca->save(psio, FileNumber, Matrix::SubBlocks);
   Cb->save(psio, FileNumber, Matrix::SubBlocks);
   psio->close(FileNumber, 1);

}
}      //End namespace psi

