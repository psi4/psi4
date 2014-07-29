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

/*!
** \file
** \brief Print a matrix of doubles
** \ingroup CIOMR
*/

#include <cstdio>
#include "psi4-dec.h"
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
void print_mat(double **a, int m, int n, FILE *out)
{

int num_frames = int(n/10);
int num_frames_rem = n%10; //adding one for changing 0->1 start
int num_frame_counter = 0;
  //for each frame
  for(num_frame_counter=0;num_frame_counter<num_frames;num_frame_counter++){
    psi::fprintf(out,"\n");
    for(int j=10*num_frame_counter+1;j<10*num_frame_counter+11;j++){
       if(j==10*num_frame_counter+1){ psi::fprintf(out,"%18d",j); }
       else{ psi::fprintf(out,"        %5d",j); }
    }
    psi::fprintf(out,"\n\n");

    for(int k=1; k<=m; ++k){
      for(int j=10*num_frame_counter+1;j<10*num_frame_counter+12;j++){
         if(j==10*num_frame_counter+1){ psi::fprintf(out,"%5d",k);}
         else{ psi::fprintf(out," %12.7f",a[k-1][j-2]); }
      }
      psi::fprintf(out,"\n");
    }
  }

// ALREADY DID THE FULL FRAMES BY THIS POINT
// NEED TO TAKE CARE OF THE REMAINDER
if(num_frames_rem != 0){
  psi::fprintf(out,"\n");
  for(int j=10*num_frame_counter+1;j<=n;j++){
       if(j==10*num_frame_counter+1){ psi::fprintf(out,"%18d",j); }
       else{ psi::fprintf(out,"        %5d",j); }
  }
  psi::fprintf(out,"\n\n");

  for(int k=1; k<=m; ++k){
    for(int j=10*num_frame_counter+1;j<n+2;j++){
         if(j==10*num_frame_counter+1){ psi::fprintf(out,"%5d",k); }
         else{ psi::fprintf(out," %12.7f",a[k-1][j-2]); }
      }
      psi::fprintf(out,"\n");
  }
}
psi::fprintf(out,"\n\n");
//R.I.P. goto statements - Aug 4th 2010 - MSM

}

}

