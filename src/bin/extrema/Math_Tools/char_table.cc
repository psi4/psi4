/*###########################################################################*/
/*! 
** \file
** \ingroup EXTREMA
** \brief Character table class. 
**
** This file contains functions which provide information
**  from the group character tables as given in Cotton 
**             Rollin King       1996
*/
/*##########################################################################*/

#include<cstdlib>
#include<cstring>
#include<cstdio>
#include <cctype>
#include"math_tools.h"

using namespace psi::extrema;

/*--------------------------------------------------------------------------*/
/*! \fn char_table::char_table(char *ptgrp_name) 
  \brief Initializes character table object.
  \param ptgrp_name name of point group
/*--------------------------------------------------------------------------*/

char_table::char_table(char *ptgrp_name) {

    int i;
    for(i=0;i<3;++i)
	ptgrp_name[i] = toupper(ptgrp_name[i]);
    ptgrp = ptgrp_name;
    
    num_irreps = get_num_irreps();
    ctable = (int**) malloc(num_irreps*sizeof(int*));
    for(i=0;i<num_irreps;++i)
	ctable[i] = (int*) malloc(num_irreps*sizeof(int));
    get_char_table();
    sym_ops = get_sym_ops();
    irrep_labels = get_irrep_labels();
    ops_coeffs = get_ops_coeffs();
    num_ops = get_num_ops();
    num_classes = get_num_classes();
}

char_table :: ~char_table() {
    int i;
    for(i=0;i<num_irreps;++i) 
	free(ctable[i]); 
    free(ctable);
    return;
}

/*--------------------------------------------------------------------------*/
/*! \fn char_table::get_num_irreps()
  \brief Returns number of irreps in point group. */
/*--------------------------------------------------------------------------*/

int char_table::get_num_irreps() {

    if (strcmp(ptgrp,"C1 ") == 0)
	return 1;
    else if ((strcmp(ptgrp,"CS ") == 0) || (strcmp(ptgrp,"CI ") == 0) 
	     || (strcmp(ptgrp,"C2 ") == 0))
	return 2;
    else if ((strcmp(ptgrp,"C2V") == 0) || (strcmp(ptgrp,"D2 ") == 0))
	return 4;
    else if (strcmp(ptgrp,"C2H") == 0)
	return 4;
    else if ((strcmp(ptgrp,"C3V") == 0) || (strcmp(ptgrp,"D3 ") == 0))
	return 3;
    else if ((strcmp(ptgrp,"C5V") == 0) || (strcmp(ptgrp,"D5") == 0))
	return 4;
    else if ((strcmp(ptgrp,"C4V") == 0) || (strcmp(ptgrp,"D2D") == 0) 
	     || (strcmp(ptgrp,"D4 ") == 0))
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



/*--------------------------------------------------------------------------*/
/*! \fn char_table::get_char_table()
  /brief Returns character table.
/*--------------------------------------------------------------------------*/

void char_table::get_char_table() {

   int i,j,count;
 
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

   
   count =0;

if (strcmp(ptgrp,"C1 ") == 0)
   ctable[0][0] = 1;

 else if ((strcmp(ptgrp,"CS ") == 0) || (strcmp(ptgrp,"CI ") == 0) 
	  || (strcmp(ptgrp,"C2 ") == 0)) {
     for (i=0;i<num_irreps;++i)
	 for (j=0;j<num_irreps;++j)
	     ctable[i][j] = CS[count++];
   }
 else if ((strcmp(ptgrp,"D2 ") == 0) || (strcmp(ptgrp,"C2V") == 0)) {
     for (i=0;i<num_irreps;++i)
	 for (j=0;j<num_irreps;++j)
	     ctable[i][j] = D2[count++];
   }
 else if (strcmp(ptgrp,"C2H") == 0) {
     for (i=0;i<num_irreps;++i)
	 for (j=0;j<num_irreps;++j)
	     ctable[i][j] = C2H[count++];
   }
 else if ((strcmp(ptgrp,"C3V") == 0) || (strcmp(ptgrp,"D3 ") == 0)) {
     for (i=0;i<num_irreps;++i)
	 for (j=0;j<num_irreps;++j)
	     ctable[i][j] = D3[count++];
   }
 else if ((strcmp(ptgrp,"C4V") == 0) || (strcmp(ptgrp,"D2D") == 0) 
	  || (strcmp(ptgrp,"D4 ") == 0)) {
     for (i=0;i<num_irreps;++i)
	 for (j=0;j<num_irreps;++j)
	     ctable[i][j] = D4[count++];
 }
 else if ((strcmp(ptgrp,"O  ") == 0) || (strcmp(ptgrp,"TD ") == 0)) {
     for (i=0;i<num_irreps;++i)
	 for (j=0;j<num_irreps;++j)
	     ctable[i][j] = TD[count++];
 }
 else if ((strcmp(ptgrp,"C6V") == 0) || (strcmp(ptgrp,"D6 ") == 0)) {
     for (i=0;i<num_irreps;++i)
	 for (j=0;j<num_irreps;++j)
	     ctable[i][j] = D6[count++];
 }
 else if ((strcmp(ptgrp,"D3D") == 0) || (strcmp(ptgrp,"D3H") == 0)) {
     for (i=0;i<num_irreps;++i)
	 for (j=0;j<num_irreps;++j)
	     ctable[i][j] = D3D[count++];
 }
 else if (strcmp(ptgrp,"D2H") == 0) {
     for (i=0;i<num_irreps;++i)
	 for(j=0;j<num_irreps;++j)
	     ctable[i][j] = D2H[count++];
 }
 else if (strcmp(ptgrp,"D4H") == 0) {
     for (i=0;i<num_irreps;++i)
	 for (j=0;j<num_irreps;++j)
	     ctable[i][j] = D4H[count++];
 }
 else if (strcmp(ptgrp,"OH ") == 0) {
     for (i=0;i<num_irreps;++i)
	 for (j=0;j<num_irreps;++j)
	     ctable[i][j] = OH[count++];
 }
 else if (strcmp(ptgrp,"D6H") == 0) {
     for (i=0;i<num_irreps;++i)
	 for (j=0;j<num_irreps;++j)
	     ctable[i][j] = D6H[count++];
 }
 else
     ctable[0][0] = 1;

  return;
}



/*--------------------------------------------------------------------------*/
/*! char_table::get_sym_ops()
  \brief Returns labels for symmetry operations. */
/*---------------------------------------------------------------------------*/

const char **char_table::get_sym_ops() {

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
   static const char *D4H[] = {"E","C4","C2","C2'","C2\"","I",
			 "S4","SGH","SGV","SGD"}; 
   static const char *OH[] =  {"E","C3","C2","C4","C2","I","S4","S6","SGH","SGD"};
   static const char *D6H[] = {"E","C6","C3","C2","C2'","C2\"","I",
			 "S3","S6","SGH","SGD","SGV"};

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



/*--------------------------------------------------------------------------*/
/*! \fn get_irrep_labels()
  \brief Returns irrep labels.
/*---------------------------------------------------------------------------*/
const char **char_table::get_irrep_labels() {

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
   static const char *D4H[] = {"A1g","A2g","B1g","B2g","Eg",
			 "A1u","A2u","B1u","B2u","Eu"};
   static const char *OH[] = {"A1g","A2g","Eg","T1g","T2g",
			"A1u","A2u","Eu","T1u","T2u"};
   static const char *D6H[] = {"A1g","A2g","B1g","B2g","E1g","E2g",
			 "A1u","A2u","B1u","B2u","E1u","E2u"};

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



/*--------------------------------------------------------------------------*/
/*! \fn char_table::get_ops_coeffs()
  \brief Returns array containing numbers of operations in each class.
/*--------------------------------------------------------------------------*/

int *char_table::get_ops_coeffs() {

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



/*--------------------------------------------------------------------------*/
/*! \fn char_table::get_num_ops()
  \brief Returns total number of symmetry operations. */
/*---------------------------------------------------------------------------*/

int char_table::get_num_ops() {
    
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



/*--------------------------------------------------------------------------*/
/*! \fn char_table::get_num_classes()
  \brief Returns the number of classes in a point group. */
/*--------------------------------------------------------------------------*/

int char_table::get_num_classes() {
    
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




    





















