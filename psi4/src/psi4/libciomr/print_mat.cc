/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
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

/*!
** \file
** \brief Print a matrix of doubles
** \ingroup CIOMR
*/

#include <cstdio>
#include "psi4/psi4-dec.h"
#include "psi4/libparallel/ParallelPrinter.h"
namespace psi {

/*!
** print_mat: Print a matrix a of dimensions mxn to file pointer out.
**
** \param a   = matrix to print
** \param m   = number of rows in matrix
** \param n   = number of columns in matrix
** \param out = file pointer for output
**
** Returns: none
**
** \ingroup CIOMR
*/
void print_mat(double **a, int m, int n, std::string out)
{
   std::shared_ptr<psi::PsiOutStream> printer=(out=="outfile"?outfile:
         std::shared_ptr<OutFile>(new OutFile(out)));
int num_frames = int(n/10);
int num_frames_rem = n%10; //adding one for changing 0->1 start
int num_frame_counter = 0;
  //for each frame
  for(num_frame_counter=0;num_frame_counter<num_frames;num_frame_counter++){
    printer->Printf("\n");
    for(int j=10*num_frame_counter+1;j<10*num_frame_counter+11;j++){
       if(j==10*num_frame_counter+1){ printer->Printf("%18d",j); }
       else{ printer->Printf("        %5d",j); }
    }
    printer->Printf("\n\n");

    for(int k=1; k<=m; ++k){
      for(int j=10*num_frame_counter+1;j<10*num_frame_counter+12;j++){
         if(j==10*num_frame_counter+1){ printer->Printf("%5d",k);}
         else{ printer->Printf(" %12.7f",a[k-1][j-2]); }
      }
      printer->Printf("\n");
    }
  }

// ALREADY DID THE FULL FRAMES BY THIS POINT
// NEED TO TAKE CARE OF THE REMAINDER
if(num_frames_rem != 0){
  printer->Printf("\n");
  for(int j=10*num_frame_counter+1;j<=n;j++){
       if(j==10*num_frame_counter+1){ printer->Printf("%18d",j); }
       else{ printer->Printf("        %5d",j); }
  }
  printer->Printf("\n\n");

  for(int k=1; k<=m; ++k){
    for(int j=10*num_frame_counter+1;j<n+2;j++){
         if(j==10*num_frame_counter+1){ printer->Printf("%5d",k); }
         else{ printer->Printf(" %12.7f",a[k-1][j-2]); }
      }
      printer->Printf("\n");
  }
}
printer->Printf("\n\n");
//R.I.P. goto statements - Aug 4th 2010 - MSM

}

}
