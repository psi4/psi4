/*
 * JKFactory: Interface and code for highly parallel J and K
 *             builds.
 *
 *  Copyright (c) 2014 Ryan M. Richard
 *
 *  This file is part of JKFactory.
 *
 *  JKFactory is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "Matrix.hpp"
#include "MPIManager.hpp"
#include<cmath>
extern "C" {
int numroc_(int *, int *, int*, int*, int*);
}
namespace JKFactory {

void MatrixMPIData::SplitProc() {
   std::string Comm=MPI->Comm();
   int NProc=MPI->NProc(Comm);
   NPRow=(int)floor(std::sqrt((double)NProc));
   bool done=false;
   while (!done) {
      if (NProc%NPRow==0) {
         NPCol=NProc/NPRow;
         done=true;
      }
      else NPRow--;
      if (NPRow==0)
         Error("NProc%(nrows=1)==0, so how did we trigger else nrows--?");
   }
}

void MatrixMPIData::DivvyMatrix(int NBasis) {
   std::string Comm=MPI->Comm();
   int myrank=MPI->Me(Comm);
   int rowid=myrank/NPCol;
   int colid=myrank%NPCol;
   MPI->MakeComm(Comm, rowid, colid, "ROW_COMM");
   MPI->MakeComm(Comm, colid, rowid, "COL_COMM");
   int mycol=MPI->Me("ROW_COMM"),myrow=MPI->Me("COL_COMM");
   int izero=0,nb=std::min(NBasis/NPRow, NBasis/NPCol);
   nBlkRow=numroc_(&NBasis, &nb, &myrow, &izero, &NPRow);
   nBlkCol=numroc_(&NBasis, &nb, &mycol, &izero, &NPCol);
   RowPerBlock=boost::shared_ptr<int[]>(new int[NPRow]);
   ColPerBlock=boost::shared_ptr<int[]>(new int[NPCol]);
   if (RowPerBlock==NULL||ColPerBlock==NULL)
      Error("Could not allocate memory for RowStarts or ColStarts");
   MPI->AllGather(&nBlkRow, 1, &RowPerBlock[0], "COL_COMM");
   MPI->AllGather(&nBlkCol, 1, &ColPerBlock[0], "ROW_COMM");
   for (int i=0; i<myrow; i++)
      startrow+=RowPerBlock[i];
   endrow=startrow+nBlkRow-1;
   for (int i=0; i<mycol; i++)
      startcol+=ColPerBlock[i];
   endcol=startcol+nBlkCol-1;
}

MatrixMPIData::MatrixMPIData(int NBasis) :
      NPRow(0), NPCol(0),startrow(0),startcol(0),endrow(0),endcol(0),
      nBlkRow(0), nBlkCol(0) {
   SplitProc();
   DivvyMatrix(NBasis);
   int stride=endcol-startcol+1;

   MyBuffer=boost::shared_ptr<double[]>(new double[stride*(endrow-startrow+1)]);
}

double* Matrix::GetMyBlock() {
   for (int i=MyData->startrow; i<=MyData->endrow; i++) {
      for (int j=MyData->startcol; j<=MyData->endcol; j++) {
         (*MyData)(i,j)=(*this)(i,j);
      }
   }
   return &(MyData->MyBuffer[0]);
}

Matrix::Matrix(double* Matrix2Store, int Dim) :
      GatheredData(Matrix2Store), nbasis(Dim), MyData(new MatrixMPIData(Dim))
{
}

void Matrix::PrintOut() const {
   for (int row=0; row<nbasis; row++) {
      for (int col=0; col<nbasis; col++) {
         (*this)<<(*this)(row, col);
         this->Print(" ");
      }
      this->Print("\n");
   }
}

void Matrix::Gather(){
   std::string Comm=MPI->Comm();
   int NElem=0;
   double* buffer=NULL;
   int startr=0,endr=0;
   int startc=0,endc=0;
   for(int i=0;i<MPI->NProc(Comm);i++){
      if(MPI->Me(Comm)==i){
         NElem=MyData->nelements();
         buffer=&(MyData->MyBuffer[0]);
         startr=MyData->startrow;
         endr=MyData->endrow;
         startc=MyData->startcol;
         endc=MyData->endcol;
      }
      MPI->Bcast(&NElem,1,i,Comm);
      MPI->Bcast(&startr,1,i,Comm);
      MPI->Bcast(&endr,1,i,Comm);
      MPI->Bcast(&startc,1,i,Comm);
      MPI->Bcast(&endc,1,i,Comm);
      MPI->Sync(Comm);
      if(MPI->Me(Comm)!=i)
         buffer=new double[NElem];
      MPI->Bcast(buffer,NElem,i,Comm);
      int stride=endc-startc+1;
      for (int j=startr; j<=endr; j++) {
         for (int k=startc; k<=endc; k++) {
            int ii=j-startr,jj=k-startc;
            (*this)(j,k)=buffer[ii*stride+jj];
         }
      }
      if(MPI->Me(Comm)!=i)
         delete [] buffer;
      MPI->Sync(Comm);
   }
}

} //End namespace JKFactory

