#include "includes.h"

extern "C" {

void sort_buffer_o(buffer,size)
   struct pkints {
      int ij;
      int kl;
      double p;
      double k;
      } *buffer;
   int size;

{
   int *indices;
   int ijkl,i,ij,kl;

   indices = (int *) malloc(sizeof(int)*size);

   for(i=0; i < size ; i++) {
      ij = buffer[i].ij;
      kl = buffer[i].kl;
      indices[i] = ij*(ij+1)/2 + kl;
      }

   qsort_o(indices,buffer,0,size-1);

   free(indices);
   }

void qsort_o(indices,buffer,left,right)
   struct pkints {
      int ij;
      int kl;
      double p;
      double k;
      } *buffer;
   int *indices,left,right;

{
   register int i,j,x,y;
   struct pkints bx,by;

   i = left;
   j = right;
   x = indices[(left+right)/2];
   bx = buffer[(left+right)/2];

   do {
      while(indices[i] < x && i < right) i++;
      while(x < indices[j] && j > left) j--;

      if(i <= j) {
         y = indices[i]; by = buffer[i];
         indices[i] = indices[j];
         buffer[i] = buffer[j];
         indices[j] = y;
         buffer[j] = by;
         i++; j--;
         }
      } while(i <= j);

   if(left < j) qsort_o(indices,buffer,left,j);
   if(i < right) qsort_o(indices,buffer,i,right);
}

} /* extern "C" */
