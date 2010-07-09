/*! \file
    \ingroup OPTKING
    \brief This file contains functions which provide information
     from the group character tables as given in Cotton 
*/

#define EXTERN
#include "globals.h"
#undef EXTERN
#include "cartesians.h"
#include "simples.h"
#include "salc.h"
#include "opt.h"

namespace psi { //namespace optking {

int **get_char_table(char *ptgrp) {

   int nirreps,i,j;
   int **table;
   int count;

   static int C1[1] = {1};

/* This works for CS, CI and C2 */
   static int CS[4] = { 1, 1,
                        1,-1};

/* This works for C2V and D2 */
   static int D2[16] = { 1, 1, 1, 1,
                         1, 1,-1,-1,
                         1,-1, 1,-1,
                         1,-1,-1, 1};

   static int C2H[16] = { 1, 1, 1, 1,
                          1,-1, 1,-1,
                          1, 1,-1,-1,
                          1,-1,-1, 1};

/* This works for C3V and D3 */
   static int D3[9] = { 1, 1, 1,
                        1, 1,-1,
                        2,-1, 0};

/* This works for C4V, D2D and D4 */
   static int D4[25] = { 1, 1, 1, 1, 1,
                         1, 1, 1,-1,-1,
                         1,-1, 1, 1,-1,
                         1,-1, 1,-1, 1,
                         2, 0,-2, 0, 0};

/* This works for O and TD point groups */
   static int TD[25] = { 1, 1, 1, 1, 1,
                           1, 1, 1,-1,-1,
                           2,-1, 2, 0, 0,
                           3, 0,-1, 1,-1,
                           3, 0,-1,-1, 1};

/* This works for C6V and D6 */
   static int D6[36] = { 1, 1, 1, 1, 1, 1,
                         1, 1, 1, 1,-1,-1,
                         1,-1, 1,-1, 1,-1,
                         1,-1, 1,-1,-1, 1,
                         2, 1,-1,-2, 0, 0,
                         2,-1,-1, 2, 0, 0};

/* This works for D3D and D3H */
   static int D3D[36] = { 1, 1, 1, 1, 1, 1,
                          1, 1,-1, 1, 1,-1,
                          2,-1, 0, 2,-1, 0,
                          1, 1, 1,-1,-1,-1,
                          1, 1,-1,-1,-1, 1,
                          2,-1, 0,-2, 1, 0};

   static int D2H[64] = { 1, 1, 1, 1, 1, 1, 1, 1,
                          1, 1,-1,-1, 1, 1,-1,-1,
                          1,-1, 1,-1, 1,-1, 1,-1,
                          1,-1,-1, 1, 1,-1,-1, 1,
                          1, 1, 1, 1,-1,-1,-1,-1,
                          1, 1,-1,-1,-1,-1, 1, 1,
                          1,-1, 1,-1,-1, 1,-1, 1,
                          1,-1,-1, 1,-1, 1, 1,-1};


   static int D4H[100] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                           1, 1, 1,-1,-1, 1, 1, 1,-1,-1,
                           1,-1, 1, 1,-1, 1,-1, 1, 1,-1,
                           1,-1, 1,-1, 1, 1,-1, 1,-1, 1,
                           2, 0,-2, 0, 0, 2, 0,-2, 0, 0,
                           1, 1, 1, 1, 1,-1,-1,-1,-1,-1,
                           1, 1, 1,-1,-1,-1,-1,-1, 1, 1,
                           1,-1, 1, 1,-1,-1, 1,-1,-1, 1,
                           1,-1, 1,-1, 1,-1, 1,-1, 1,-1,
                           2, 0,-2, 0, 0,-2, 0, 2, 0, 0};

   static int OH[100] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                          1, 1,-1,-1, 1, 1,-1, 1, 1,-1,
                          2,-1, 0, 0, 2, 2, 0,-1, 2, 0,
                          3, 0,-1, 1,-1, 3, 1, 0,-1,-1,
                          3, 0, 1,-1,-1, 3,-1, 0,-1, 1,
                          1, 1, 1, 1, 1,-1,-1,-1,-1,-1,
                          1, 1,-1,-1, 1,-1, 1,-1,-1, 1,
                          2,-1, 0, 0, 2,-2, 0, 1,-2, 0,
                          3, 0,-1, 1,-1,-3,-1, 0, 1, 1,
                          3, 0, 1,-1,-1,-3, 1, 0, 1,-1};

   static int D6H[144] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                           1, 1, 1, 1,-1,-1, 1, 1, 1, 1,-1,-1,
                           1,-1, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1,
                           1,-1, 1,-1,-1, 1, 1,-1, 1,-1,-1, 1,
                           2, 1,-1,-2, 0, 0, 2, 1,-1,-2, 0, 0,
                           2,-1,-1, 2, 0, 0, 2,-1,-1, 2, 0, 0,
                           1, 1, 1, 1, 1, 1,-1,-1,-1,-1,-1,-1,
                           1, 1, 1, 1,-1,-1,-1,-1,-1,-1, 1, 1,
                           1,-1, 1,-1, 1,-1,-1, 1,-1, 1,-1, 1,
                           1,-1, 1,-1,-1, 1,-1, 1,-1, 1, 1,-1,
                           2, 1,-1,-2, 0, 0,-2,-1, 1, 2, 0, 0,
                           2,-1,-1, 2, 0, 0,-2, 1, 1,-2, 0, 0};

   nirreps = get_nirreps(ptgrp);
   table = init_int_matrix(nirreps,nirreps);

   count = 0;

if (strcmp(ptgrp,"C1 ") == 0)
   table[0][0] = 1;

else if ((strcmp(ptgrp,"CS ") == 0) || (strcmp(ptgrp,"CI ") == 0) || (strcmp(ptgrp,"C2 ") == 0)) {
   for (i=0;i<nirreps;++i)
      for (j=0;j<nirreps;++j)
         table[i][j] = CS[count++];
   }
else if ((strcmp(ptgrp,"D2 ") == 0) || (strcmp(ptgrp,"C2V") == 0)) {
   for (i=0;i<nirreps;++i)
      for (j=0;j<nirreps;++j)
         table[i][j] = D2[count++];
   }
else if (strcmp(ptgrp,"C2H") == 0) {
   for (i=0;i<nirreps;++i)
      for (j=0;j<nirreps;++j)
         table[i][j] = C2H[count++];
   }
else if ((strcmp(ptgrp,"C3V") == 0) || (strcmp(ptgrp,"D3 ") == 0)) {
   for (i=0;i<nirreps;++i)
      for (j=0;j<nirreps;++j)
         table[i][j] = D3[count++];
   }
else if ((strcmp(ptgrp,"C4V") == 0) || (strcmp(ptgrp,"D2D") == 0) || (strcmp(ptgrp,"D4 ") == 0)) {
   for (i=0;i<nirreps;++i)
      for (j=0;j<nirreps;++j)
         table[i][j] = D4[count++];
   }
else if ((strcmp(ptgrp,"O  ") == 0) || (strcmp(ptgrp,"TD ") == 0)) {
   for (i=0;i<nirreps;++i)
      for (j=0;j<nirreps;++j)
         table[i][j] = TD[count++];
   }
else if ((strcmp(ptgrp,"C6V") == 0) || (strcmp(ptgrp,"D6 ") == 0)) {
   for (i=0;i<nirreps;++i)
      for (j=0;j<nirreps;++j)
         table[i][j] = D6[count++];
   }
else if ((strcmp(ptgrp,"D3D") == 0) || (strcmp(ptgrp,"D3H") == 0)) {
   for (i=0;i<nirreps;++i)
      for (j=0;j<nirreps;++j)
         table[i][j] = D3D[count++];
   }
else if (strcmp(ptgrp,"D2H") == 0) {
   for (i=0;i<nirreps;++i)
      for (j=0;j<nirreps;++j)
         table[i][j] = D2H[count++];
   }
else if (strcmp(ptgrp,"D4H") == 0) {
   for (i=0;i<nirreps;++i)
      for (j=0;j<nirreps;++j)
         table[i][j] = D4H[count++];
   }
else if (strcmp(ptgrp,"OH ") == 0) {
   for (i=0;i<nirreps;++i)
      for (j=0;j<nirreps;++j)
         table[i][j] = OH[count++];
   }
else if (strcmp(ptgrp,"D6H") == 0) {
   for (i=0;i<nirreps;++i)
      for (j=0;j<nirreps;++j)
         table[i][j] = D6H[count++];
   }
else
   table[0][0] = 1;



  return table;
}




int get_nirreps(char *ptgrp) {

if (strcmp(ptgrp,"C1 ") == 0)
   return 1;
else if ((strcmp(ptgrp,"CS ") == 0) || (strcmp(ptgrp,"CI ") == 0) || (strcmp(ptgrp,"C2 ") == 0))
   return 2;
else if ((strcmp(ptgrp,"C2V") == 0) || (strcmp(ptgrp,"D2 ") == 0))
   return 4;
else if (strcmp(ptgrp,"C2H") == 0)
   return 4;
else if ((strcmp(ptgrp,"C3V") == 0) || (strcmp(ptgrp,"D3 ") == 0))
   return 3;
else if ((strcmp(ptgrp,"C5V") == 0) || (strcmp(ptgrp,"D5") == 0))
   return 4;
else if ((strcmp(ptgrp,"C4V") == 0) || (strcmp(ptgrp,"D2D") == 0) || (strcmp(ptgrp,"D4 ") == 0))
   return 5;
else if ((strcmp(ptgrp,"O  ") == 0) || (strcmp(ptgrp,"TD ") == 0))
   return 5;
else if ((strcmp(ptgrp,"C6V") == 0) || (strcmp(ptgrp,"D6 ") == 0))
   return 6;
else if ((strcmp(ptgrp,"D3D") == 0) || (strcmp(ptgrp,"D3H") == 0))
   return 6;
else if (strcmp(ptgrp,"D4D") == 0)
   return 7;
else if (strcmp(ptgrp,"D2H") == 0)
   return 8;
else if ((strcmp(ptgrp,"D5D") == 0) || (strcmp(ptgrp,"D5H") == 0))
   return 8;
else if (strcmp(ptgrp,"D6D") == 0)
   return 9;
else if (strcmp(ptgrp,"D4H") == 0)
   return 10;
else if (strcmp(ptgrp,"OH ") == 0)
   return 10;
else if (strcmp(ptgrp,"D6H") == 0)
   return 12;
else
   return 1;
}

const char **get_symm_ops(char *ptgrp) {

   int i,j,nirreps;
   char **irreps;

   static const char *C1[] =  {"E"};
   static const char *CS[] =  {"E","SGH"};
   static const char *CI[] =  {"E","I"};
   static const char *C2[] =  {"E","C2"};
   static const char *D2[] =  {"E","C2Z","C2Y","C2X"};
   static const char *C2V[] = {"E","C2","SGV","SGV"};
   static const char *C2H[] = {"E","C2","I","SGH"};
   static const char *D3[] =  {"E","C3","C2"};
   static const char *C3V[] = {"E","C3","SGV"};
   static const char *D4[] =  {"E","C4","C2","C2'","C2\""};
   static const char *D2D[] = {"E","S4","C2","C2'","SGD"};
   static const char *C4V[] = {"E","C4","C2","SGV","SGD"};
   static const char *TD[] =  {"E","C3","C2","S4","SGD"};
   static const char *O[] =   {"E","C3","C2","C4","C2"};
   static const char *D6[] =  {"E","C6","C3","C2","C2'","C2\""};
   static const char *C6V[] = {"E","C6","C3","C2","SGV","SGD"};
   static const char *D3D[] = {"E","C3","C2","I","S6","SGD"};
   static const char *D3H[] = {"E","C3","C2","SGH","S3","SGV"};
   static const char *D2H[] = {"E","C2Z","C2Y","C2X","I","SGXY","SGXZ","SGYZ"}; 
   static const char *D4H[] = {"E","C4","C2","C2'","C2\"","I","S4","SGH","SGV","SGD"}; 
   static const char *OH[] =  {"E","C3","C2","C4","C2","I","S4","S6","SGH","SGD"};
   static const char *D6H[] = {"E","C6","C3","C2","C2'","C2\"","I","S3","S6","SGH","SGD","SGV"};

if (strcmp(ptgrp,"C1 ") == 0)
   return C1;
else if (strcmp(ptgrp,"CS ") == 0)
   return CS;
else if (strcmp(ptgrp,"CI ") == 0)
   return CI;
else if (strcmp(ptgrp,"C2 ") == 0)
   return C2;
else if (strcmp(ptgrp,"C2V") == 0)
   return C2V;
else if (strcmp(ptgrp,"D2 ") == 0)
   return D2;
else if (strcmp(ptgrp,"C2H") == 0)
   return C2H;
else if (strcmp(ptgrp,"C3V") == 0) 
   return C3V;
else if (strcmp(ptgrp,"D3 ") == 0)
   return D3;
else if (strcmp(ptgrp,"C4V") == 0) 
   return C4V;
else if (strcmp(ptgrp,"D2D") == 0)
   return D2D;
else if (strcmp(ptgrp,"D4 ") == 0)
   return D4;
else if (strcmp(ptgrp,"O  ") == 0) 
   return O;
else if (strcmp(ptgrp,"TD ") == 0)
   return TD;
else if (strcmp(ptgrp,"C6V") == 0) 
   return C6V;
else if (strcmp(ptgrp,"D6 ") == 0)
   return D6;
else if (strcmp(ptgrp,"D3D") == 0) 
   return D3D;
else if (strcmp(ptgrp,"D3H") == 0)
   return D3H;
else if (strcmp(ptgrp,"D2H") == 0)
   return D2H;
else if (strcmp(ptgrp,"D4H") == 0)
   return D4H;
else if (strcmp(ptgrp,"OH ") == 0)
   return OH;
else if (strcmp(ptgrp,"D6H") == 0)
   return D6H;
else
   return C1;
}




const char **get_irrep_labels(char *ptgrp) {

   int i,j,nirreps;
   char **irreps;

   static const char *C1[] = {"A"};
   static const char *CS[] = {"Ap","App"};
   static const char *CI[] = {"Ag","Au"};
   static const char *C2[] = {"A","B"};
   static const char *D2[] = {"A","B1","B2","B3"};
   static const char *C2V[] = {"A1","A2","B1","B2"};
   static const char *C2H[] = {"Ag","Bg","Au","Bu"};
   static const char *D3[] = {"A1","A2","E"};
   static const char *C3V[] = {"A1","A2","E"};
   static const char *D4[] = {"A1","A2","B1","B2","E"};
   static const char *D2D[] = {"A1","A2","B1","B2","E"};
   static const char *C4V[] = {"A1","A2","B1","B2","E"};
   static const char *TD[] = {"A1","A2","E","T1","T2"};
   static const char *O[] = {"A1","A2","E","T1","T2"};
   static const char *D6[] = {"A1","A2","B1","B2","E1","E2"};
   static const char *C6V[] = {"A1","A2","B1","B2","E1","E2"};
   static const char *D3D[] = {"A1g","A2g","Eg","A1u","A2u","Eu"};
   static const char *D3H[] = {"A1p","A2p","Ep","A1pp","A2pp","Epp"};
   static const char *D2H[] = {"Ag","B1g","B2g","B3g","Au","B1u","B2u","B3u"};
   static const char *D4H[] = {"A1g","A2g","B1g","B2g","Eg","A1u","A2u","B1u","B2u","Eu"};
   static const char *OH[] = {"A1g","A2g","Eg","T1g","T2g","A1u","A2u","Eu","T1u","T2u"};
   static const char *D6H[] = {"A1g","A2g","B1g","B2g","E1g","E2g","A1u","A2u","B1u","B2u","E1u","E2u"};

if (strcmp(ptgrp,"C1 ") == 0)
   return C1;
else if (strcmp(ptgrp,"CS ") == 0)
   return CS;
else if (strcmp(ptgrp,"CI ") == 0)
   return CI;
else if (strcmp(ptgrp,"C2 ") == 0)
   return C2;
else if (strcmp(ptgrp,"C2V") == 0)
   return C2V;
else if (strcmp(ptgrp,"D2 ") == 0)
   return D2;
else if (strcmp(ptgrp,"C2H") == 0)
   return C2H;
else if (strcmp(ptgrp,"C3V") == 0)
   return C3V;
else if (strcmp(ptgrp,"D3 ") == 0)
   return D3;
else if (strcmp(ptgrp,"C4V") == 0)
   return C4V;
else if (strcmp(ptgrp,"D2D") == 0)
   return D2D;
else if (strcmp(ptgrp,"D4 ") == 0)
   return D4;
else if (strcmp(ptgrp,"O  ") == 0)
   return O;
else if (strcmp(ptgrp,"TD ") == 0)
   return TD;
else if (strcmp(ptgrp,"C6V") == 0)
   return C6V;
else if (strcmp(ptgrp,"D6 ") == 0)
   return D6;
else if (strcmp(ptgrp,"D3D") == 0)
   return D3D;
else if (strcmp(ptgrp,"D3H") == 0)
   return D3H;
else if (strcmp(ptgrp,"D2H") == 0)
   return D2H;
else if (strcmp(ptgrp,"D4H") == 0)
   return D4H;
else if (strcmp(ptgrp,"OH ") == 0)
   return OH;
else if (strcmp(ptgrp,"D6H") == 0)
   return D6H;
else
   return C1;
}




int *get_ops_coeffs(char *ptgrp) {

   static int C1[]  = {1}; 
   static int CS[]  = {1,1};
   static int CI[]  = {1,1};
   static int C2[]  = {1,1};
   static int D2[]  = {1,1,1,1};
   static int C2V[] = {1,1,1,1};
   static int C2H[] = {1,1,1,1};
   static int D3[]  = {1,2,3};
   static int C3V[] = {1,2,3};
   static int D4[]  = {1,2,1,2,2};
   static int D2D[] = {1,2,1,2,2};
   static int C4V[] = {1,2,1,2,2};
   static int TD[]  = {1,8,3,6,6};
   static int O[]   = {1,8,3,6,6};
   static int D6[]  = {1,2,2,1,3,3};
   static int C6V[] = {1,2,2,1,3,3};
   static int D3D[] = {1,2,3,1,2,3};
   static int D3H[] = {1,2,3,1,2,3};
   static int D2H[] = {1,1,1,1,1,1,1,1};
   static int D4H[] = {1,2,1,2,2,1,2,1,2,2};
   static int OH[]  = {1,8,6,6,3,1,6,8,3,6};
   static int D6H[] = {1,2,2,1,3,3,1,2,2,1,3,3};

if (strcmp(ptgrp,"C1 ") == 0)
   return C1;
else if (strcmp(ptgrp,"CS ") == 0)
   return CS;
else if (strcmp(ptgrp,"CI ") == 0)
   return CI;
else if (strcmp(ptgrp,"C2 ") == 0)
   return C2;
else if (strcmp(ptgrp,"C2V") == 0)
   return C2V;
else if (strcmp(ptgrp,"D2 ") == 0)
   return D2;
else if (strcmp(ptgrp,"C2H") == 0)
   return C2H;
else if (strcmp(ptgrp,"C3V") == 0)
   return C3V;
else if (strcmp(ptgrp,"D3 ") == 0)
   return D3;
else if (strcmp(ptgrp,"C4V") == 0)
   return C4V;
else if (strcmp(ptgrp,"D2D") == 0)
   return D2D;
else if (strcmp(ptgrp,"D4 ") == 0)
   return D4;
else if (strcmp(ptgrp,"O  ") == 0)
   return O;
else if (strcmp(ptgrp,"TD ") == 0)
   return TD;
else if (strcmp(ptgrp,"C6V") == 0)
   return C6V;
else if (strcmp(ptgrp,"D6 ") == 0)
   return D6;
else if (strcmp(ptgrp,"D3D") == 0)
   return D3D;
else if (strcmp(ptgrp,"D3H") == 0)
   return D3H;
else if (strcmp(ptgrp,"D2H") == 0)
   return D2H;
else if (strcmp(ptgrp,"D4H") == 0)
   return D4H;
else if (strcmp(ptgrp,"OH ") == 0)
   return OH;
else if (strcmp(ptgrp,"D6H") == 0)
   return D6H;
else
   return C1;
}




int get_num_ops(char *ptgrp) {

if (strcmp(ptgrp,"C1 ") == 0)
   return 1;
else if (strcmp(ptgrp,"CS ") == 0)
   return 2;
else if (strcmp(ptgrp,"CI ") == 0)
   return 2;
else if (strcmp(ptgrp,"C2 ") == 0)
   return 2;
else if (strcmp(ptgrp,"C2V") == 0)
   return 4;
else if (strcmp(ptgrp,"D2 ") == 0)
   return 4;
else if (strcmp(ptgrp,"C2H") == 0)
   return 4;
else if (strcmp(ptgrp,"C3V") == 0)
   return 6;
else if (strcmp(ptgrp,"D3 ") == 0)
   return 6;
else if (strcmp(ptgrp,"C4V") == 0)
   return 8;
else if (strcmp(ptgrp,"D2D") == 0)
   return 8;
else if (strcmp(ptgrp,"D4 ") == 0)
   return 8;
else if (strcmp(ptgrp,"O  ") == 0)
   return 24;
else if (strcmp(ptgrp,"TD ") == 0)
   return 24;
else if (strcmp(ptgrp,"C6V") == 0)
   return 12;
else if (strcmp(ptgrp,"D6 ") == 0)
   return 12;
else if (strcmp(ptgrp,"D3D") == 0)
   return 12;
else if (strcmp(ptgrp,"D3H") == 0)
   return 12;
else if (strcmp(ptgrp,"D2H") == 0)
   return 8;
else if (strcmp(ptgrp,"D4H") == 0)
   return 16;
else if (strcmp(ptgrp,"OH ") == 0)
   return 48;
else if (strcmp(ptgrp,"D6H") == 0)
   return 24;
else
   return 1;
}




int get_num_classes(char *ptgrp) {

if (strcmp(ptgrp,"C1 ") == 0)
   return 1;
else if (strcmp(ptgrp,"CS ") == 0)
   return 2;
else if (strcmp(ptgrp,"CI ") == 0)
   return 2;
else if (strcmp(ptgrp,"C2 ") == 0)
   return 2;
else if (strcmp(ptgrp,"C2V") == 0)
   return 4;
else if (strcmp(ptgrp,"D2 ") == 0)
   return 4;
else if (strcmp(ptgrp,"C2H") == 0)
   return 4;
else if (strcmp(ptgrp,"C3V") == 0)
   return 3;
else if (strcmp(ptgrp,"D3 ") == 0)
   return 3;
else if (strcmp(ptgrp,"C4V") == 0)
   return 5;
else if (strcmp(ptgrp,"D2D") == 0)
   return 5;
else if (strcmp(ptgrp,"D4 ") == 0)
   return 5;
else if (strcmp(ptgrp,"O  ") == 0)
   return 5;
else if (strcmp(ptgrp,"TD ") == 0)
   return 5;
else if (strcmp(ptgrp,"C6V") == 0)
   return 6;
else if (strcmp(ptgrp,"D6 ") == 0)
   return 6;
else if (strcmp(ptgrp,"D3D") == 0)
   return 6;
else if (strcmp(ptgrp,"D3H") == 0)
   return 6;
else if (strcmp(ptgrp,"D2H") == 0)
   return 8;
else if (strcmp(ptgrp,"D4H") == 0)
   return 10;
else if (strcmp(ptgrp,"OH ") == 0)
   return 10;
else if (strcmp(ptgrp,"D6H") == 0)
   return 12;
else
   return 1;
}


/*** GET_OPS_IN_CLASS : finds numbers of operations in the classes
 * of the molecular point group. nirreps is for checking */

int *get_ops_in_class(char *ptgrp, int nirreps) {

  int num_class, *num_ops, error_var = 0; 

  static int class_C1[] = {1},
             num_class_C1 = 1,

             /* for CS,CI,C2 */
             class_CS[] = {1,1},
             num_class_CS = 2,

             /* for C2V,D2,C2h */
             class_D2[] = {1,1,1,1},
             num_class_D2 = 4,

             /* for C3V,D3 */
             class_D3[] = {1,2,3},
             num_class_D3 = 3,

             /* for C4V,D2D,D4 */
             class_D4[] = {1,2,1,2,2},
             num_class_D4 = 5,

             /* for O,TD */
             class_TD[] = {1,8,3,6,6},
             num_class_TD = 5,

             /* for C6V,D6 */
             class_D6[] = {1,2,2,1,3,3},
             num_class_D6 = 6,

             /*for D3D,D3H */
             class_D3D[] = {1,2,3,1,2,3},
             num_class_D3D = 6,

             class_D2H[] = {1,1,1,1,1,1,1,1},
             num_class_D2H = 8,

             class_D4H[] = {1,2,1,2,2,1,2,1,2,2},
             num_class_D4H = 10,

             class_OH[] = {1,8,6,6,3,1,6,8,3,6},
             num_class_OH = 10,

             class_D6H[] = {1,2,2,1,3,3,1,2,2,1,3,3},
             num_class_D6H = 12;


  if (strcmp(ptgrp,"C1 ") == 0) {
    if(nirreps == num_class_C1)
      num_ops = class_C1;
    else error_var = 1;
  }
  else if (strcmp(ptgrp,"CS ")==0 || strcmp(ptgrp,"CI ")==0 ||
          strcmp(ptgrp,"C2 ")==0) {
    if(nirreps == num_class_CS)
      num_ops = class_CS;
    else error_var = 1;
  }
  else if (strcmp (ptgrp, "C2V")==0 || strcmp (ptgrp, "D2 ")==0 ||
          strcmp(ptgrp,"C2H")==0) {
    if(nirreps == num_class_D2)
      num_ops = class_D2;
    else error_var = 1;
  }
  else if(strcmp(ptgrp,"C3V")==0 || strcmp(ptgrp,"D3 ")==0) {
    if(nirreps == num_class_D3)
      num_ops = class_D3;
    else error_var = 1;
  }  

  else if(strcmp(ptgrp,"C4V")==0 || strcmp(ptgrp,"D2D")==0 || 
          strcmp(ptgrp,"D4 ")==0) {
    if(nirreps == num_class_D4)
      num_ops = class_D4;
    else error_var = 1;
  }

  else if(strcmp(ptgrp,"O  ")==0 || strcmp(ptgrp,"TD ")==0) { 
    if(nirreps == num_class_TD)
      num_ops = class_TD;
    else error_var = 1;
  }

  else if(strcmp(ptgrp,"C6V")==0 || strcmp(ptgrp,"D6 ")==0) {
    if(nirreps == num_class_D6)
      num_ops = class_D6;
    else error_var = 1;
  }  

  else if(strcmp(ptgrp,"D3D")==0 || strcmp(ptgrp,"D3H")==0) {
    if(nirreps == num_class_D3D)
      num_ops = class_D3D;
    else error_var = 1;
  }

  else if(strcmp(ptgrp,"D2H")==0) {
    if(nirreps == num_class_D2H)
      num_ops = class_D2H;
    else error_var = 1;
  }

  else if(strcmp(ptgrp,"D4H")==0) {
    if(nirreps == num_class_D4H)
      num_ops = class_D4H;
    else error_var = 1;
  }

  else if(strcmp(ptgrp,"OH ")==0) {
    if(nirreps == num_class_OH)
      num_ops = class_OH;
    else error_var = 1;
  }

  else if(strcmp(ptgrp,"D6H")==0) {
    if(nirreps == num_class_D6H)
      num_ops = class_D6H;
    else error_var = 1;
  }

  else if (num_ops[0] == 0) error_var = 1;

  if(error_var != 0) {
    punt("problem assigning number of operations per class");
  } 

return num_ops;

}      

}//} /* namespace psi::optking */

