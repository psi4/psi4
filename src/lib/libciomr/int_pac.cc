#include "includes.h"

extern "C" {

void int_pac(int* i,int* ib, int* j, int* jb, unsigned int* ipak, unsigned int* jpak)
   {
      *ipak = *ib;
      *ipak <<= 4;
      *ipak += *i;


      *jpak = *jb;
      *jpak <<= 4;
      *jpak += *j;
      }

void int_unpac(int* i, int* ib, int* j, int* jb, unsigned int* ipak, unsigned int* jpak)
   {
      *ib = *ipak >> 4;
      *i = *ipak & 15;

      *jb = *jpak >> 4;
      *j = *jpak & 15;
      }

} /* extern "C" */
