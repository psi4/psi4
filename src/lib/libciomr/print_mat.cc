/*!
** \file
** \brief Print a matrix of doubles
** \ingroup CIOMR
*/

#include <cstdio>

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
    fprintf(out,"\n");
    for(int j=10*num_frame_counter+1;j<10*num_frame_counter+11;j++){
       if(j==10*num_frame_counter+1){ fprintf(out,"%18d",j); }
       else{ fprintf(out,"        %5d",j); }
    }
    fprintf(out,"\n\n");

    for(int k=1; k<=m; ++k){
      for(int j=10*num_frame_counter+1;j<10*num_frame_counter+12;j++){
         if(j==10*num_frame_counter+1){ fprintf(out,"%5d",k);}
         else{ fprintf(out," %12.7f",a[k-1][j-2]); }
      }
      fprintf(out,"\n");
    }
  }

// ALREADY DID THE FULL FRAMES BY THIS POINT
// NEED TO TAKE CARE OF THE REMAINDER
if(num_frames_rem != 0){
  fprintf(out,"\n");
  for(int j=10*num_frame_counter+1;j<=n;j++){
       if(j==10*num_frame_counter+1){ fprintf(out,"%18d",j); }
       else{ fprintf(out,"        %5d",j); }
  }
  fprintf(out,"\n\n");

  for(int k=1; k<=m; ++k){
    for(int j=10*num_frame_counter+1;j<n+2;j++){
         if(j==10*num_frame_counter+1){ fprintf(out,"%5d",k); }
         else{ fprintf(out," %12.7f",a[k-1][j-2]); }
      }
      fprintf(out,"\n");
  }

fprintf(out,"\n\n");
//R.I.P. goto statements - Aug 4th 2010 - MSM

}

}

