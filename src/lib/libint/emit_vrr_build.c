/*! \file 
    \ingroup (INT)
    \brief Enter brief description of file here 
*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "build_libint.h"
#include <libint/constants.h>

extern FILE *outfile, *vrr_header;
extern void punt(char *);
extern LibintParams_t Params;

static void declare_localv();
static void define_localv();

static char **k1, **k2, **k3;

int emit_vrr_build()
{
  int old_am = Params.old_am;
  int new_am = Params.opt_am;
  int max_class_size = Params.max_class_size;
  int am_to_inline = Params.max_am_to_inline_vrr_worker;

  FILE *code;
  int i, j, k, l, f;
  int a, b;
  int flag;
  int am[2][3];
  int am_in[2];
  int current_highest_am, to_inline;
  int nflip = 0;
  int t1, t2, t3, t4;
  int class_size;
  int type;
  int max1 = 0;
  int max2 = 0;
  int foo;
  int la, lc, lc_min, lc_max;
  int k1max, k2max, k3max;
  int split,num_subfunctions,subbatch_length;
  int curr_count,curr_subfunction;
  static char *k4[] = {"lpoz","lpon"};
  static const char *k1_suff = "o2z";
  static const char *k2_suff = "o2zn";
  static const char *k3_suff = "o2n";
  char *code_name;
  char *function_name;
  char **subfunction_name;

  k1 = (char **) malloc(new_am*sizeof(char *));
  k2 = (char **) malloc(new_am*sizeof(char *));
  k3 = (char **) malloc(new_am*sizeof(char *));
  for(i=1;i<=new_am;i++) {
    j = strlen(number[i]);
    k1[i-1] = (char*) malloc((4+j)*sizeof(char));
    k2[i-1] = (char*) malloc((5+j)*sizeof(char));
    k3[i-1] = (char*) malloc((4+j)*sizeof(char));
    strcpy(k1[i-1],number[i]);
    strcpy(k2[i-1],number[i]);
    strcpy(k3[i-1],number[i]);
    strcat(k1[i-1],k1_suff);
    strcat(k2[i-1],k2_suff);
    strcat(k3[i-1],k3_suff);
  }
  code_name = (char *) malloc(sizeof(char)*21);
  function_name = (char *) malloc(sizeof(char)*18);

  for(la=0;la<=new_am;la++) {
    lc_min = (la >= old_am + 1) ? 0 : old_am + 1;
    lc_max = new_am;
    for(lc=lc_min;lc<=lc_max;lc++) {

      /* Is this function to be made inline */
      current_highest_am = (la > lc) ? la : lc;
      to_inline = (current_highest_am <= am_to_inline) ? 1 : 0;
      
      fprintf(outfile,"  AM_a = %c  AM_c = %c\n",am_letter[la],am_letter[lc]);
      am_in[0] = la;
      am_in[1] = lc;
      if (la == 0) {
	a = 1;
	k2max = la;
	k3max = lc - 1;
      }
      else {
	a = 0;
	k2max = lc;
	k1max = la - 1;
      }
      foo = 5;
      if(a==0) foo = 4;

      class_size = ((am_in[a]+1)*(am_in[a]+2)*(am_in[a^1]+1)*(am_in[a^1]+2))/4;

      fprintf(vrr_header,"#define _BUILD_%c0%c0(Data,vp,i0,i1,i2,i3,i4) {",am_letter[la],am_letter[lc]);
      /* Decide if the routine has to be split into several routines producing "subbatches" */
      if (class_size > max_class_size) {
	split = 1;
	num_subfunctions = ceil((double)class_size/max_class_size);
	subbatch_length = 1 + class_size/num_subfunctions;
	fprintf(vrr_header," tmp = _build_%c0%c0_0(Data,vp,i0,i1,i2,i3,i4); \\\n",am_letter[la],am_letter[lc]);
	for(f=1;f<num_subfunctions;f++)
	  fprintf(vrr_header," tmp = _build_%c0%c0_%d(Data,tmp,i0,i1,i2,i3,i4); \\\n",am_letter[la],am_letter[lc],f);
	fprintf(vrr_header,"}\n");
	if (to_inline)
	  fprintf(vrr_header, "#ifndef INLINE_VRR_WORKER\n");
	for(f=0;f<num_subfunctions;f++) {
	  fprintf(vrr_header, " REALTYPE *_build_%c0%c0_%d(prim_data *Data, REALTYPE *, const REALTYPE *, const REALTYPE *, const REALTYPE *, const REALTYPE *, const REALTYPE *);\n",
		  am_letter[la],am_letter[lc],f);
	}
	if (to_inline)
	  fprintf(vrr_header, "#endif\n");
      }
      else {
	split = 0;
	fprintf(vrr_header," _build_%c0%c0(Data,vp,i0,i1,i2,i3,i4);}\n",am_letter[la],am_letter[lc]);
	if (to_inline)
	  fprintf(vrr_header,"#ifndef INLINE_VRR_WORKER\n");
	fprintf(vrr_header," void _build_%c0%c0(prim_data *, REALTYPE *, const REALTYPE *, const REALTYPE *, const REALTYPE *, const REALTYPE *, const REALTYPE *);\n",am_letter[la],am_letter[lc]);
	if (to_inline)
	  fprintf(vrr_header, "#endif\n");
      }

      if(a==0) foo = 4;

      sprintf(function_name,"build_%c0%c0",am_letter[la],am_letter[lc]);
      sprintf(code_name,"build_%c0%c0.cc",am_letter[la],am_letter[lc]);
      code = fopen(code_name,"w");

      /*target
        |I0[],I1[]
        |    |I2[],I3[]
        |    |    |   I4[]
        |    |    |    |     */
        t1 = t2 = t3 = t4 = 0;

      /* print local variable declarations */

      fprintf(code,"  /* These machine-generated functions compute a quartet of (%cs|%cs) integrals */\n\n",am_letter[la],am_letter[lc]);
      if (split) {
	subfunction_name = (char **) malloc (num_subfunctions*sizeof(char *));
	for(i=0;i<num_subfunctions;i++) {
	  subfunction_name[i] = (char *) malloc(20*sizeof(char));
	  sprintf(subfunction_name[i],"_build_%c0%c0_%d",am_letter[la],am_letter[lc],i);
	}
      }

      fprintf(code,"#include \"libint.h\"\n\n");

      if (split == 1) {
	curr_subfunction = 0;
	curr_count = 0;
	fprintf(code,"REALTYPE *%s(prim_data *Data, REALTYPE *vp, const REALTYPE *I0, const REALTYPE *I1, const REALTYPE *I2, const REALTYPE *I3, const REALTYPE *I4)\n{\n",
		subfunction_name[0]);
      }
      else
	fprintf(code,"void _%s(prim_data *Data, REALTYPE *vp, const REALTYPE *I0, const REALTYPE *I1, const REALTYPE *I2, const REALTYPE *I3, const REALTYPE *I4)\n{\n",function_name);
      declare_localv(a,k1max,k2max,k3max,code);
      define_localv(a,foo,k1max,k2max,k3max,code);
      fprintf(code,"\n");

      for(i = 0; i <= am_in[0]; i++){
	am[0][0] = am_in[0] - i;
	for(j = 0; j <= i; j++){
	  am[0][1] = i - j;
	  am[0][2] = j;

	  for(k = 0; k <= am_in[1]; k++){
	    am[1][0] = am_in[1] - k;
	    for(l = 0; l <= k; l++){
	      am[1][1] = k - l;
	      am[1][2] = l;

	      if(am[a][2]) b = 2;
	      if(am[a][1]) b = 1;
	      if(am[a][0]) b = 0;

          
	      am[a][b] = am[a][b] - 1;
	      am_in[a] = am_in[a] - 1;
	      t2 = hash(am,am_in);
	      fprintf(code, "*(vp++) = U%d%d*I0[%d] + U%d%d*I1[%d]",
		      a*2, b, t2, foo, b , t2); 
	      if(am[a][b]){
		am[a][b] = am[a][b] - 1;
		am_in[a] = am_in[a] - 1;
		t3 = hash(am,am_in);
		fprintf(code, "\n           + (%s)*(I2[%d] - (%s)*I3[%d])", 
			(a==0 ? k1[am[a][b]] : k3[am[a][b]]), 
			t3, (k4[a]), t3);
		max1 = (max1>am[a][b]+1) ? max1 : am[a][b]+1;
		am[a][b] = am[a][b] + 1;
		am_in[a] = am_in[a] + 1;
	      }
	      if(am[a^1][b]){
		am[a^1][b] = am[a^1][b] - 1;
		am_in[a^1] = am_in[a^1] - 1;
		t4 = hash(am,am_in);
		fprintf(code, "\n           + (%s)*I4[%d]", k2[am[a^1][b]], t4);
		max2 = (max2>am[a^1][b]+1) ? max2 : am[a^1][b]+1;
		am[a^1][b] = am[a^1][b] + 1;
		am_in[a^1] = am_in[a^1] + 1;
	      }
	      fprintf(code, ";\n");
	      am[a][b] = am[a][b] + 1;
	      am_in[a] = am_in[a] + 1;
		
	      t1++;
	      curr_count++;
	      if (curr_count == subbatch_length && split == 1) {
		curr_count = 0;
		curr_subfunction++;
		fprintf(code,"return vp;\n}\n\n");
		fprintf(code,"REALTYPE *%s(prim_data *Data, REALTYPE *vp, const REALTYPE *I0, const REALTYPE *I1, const REALTYPE *I2, const REALTYPE *I3, const REALTYPE *I4)\n{\n",
			subfunction_name[curr_subfunction]);
		declare_localv(a,k1max,k2max,k3max,code);
		define_localv(a,foo,k1max,k2max,k3max,code);
		fprintf(code,"\n");
	      }
	    }
	  }
	}
      }
      if (split == 1)
	fprintf(code,"return vp;\n}\n");
      else
	fprintf(code,"\n}\n");
      fclose(code);
      if (split == 1) {
	for(i=0;i<num_subfunctions;i++)
	  free(subfunction_name[i]);
	free(subfunction_name);
      }
      printf("Done with %s\n",code_name);
    }
  }
  free(function_name);
  free(code_name);
}


void declare_localv(int a, int k1max, int k2max, int k3max, FILE *code)
{
  int i;

  fprintf(code,"  REALTYPE U00, U01, U02, U10, U11, U12, U20, U21, U22;\n");
  fprintf(code,"  REALTYPE U30, U31, U32, U40, U41, U42, U50, U51, U52;\n");
  fprintf(code,"  REALTYPE lpoz = Data->poz;\n  REALTYPE lpon = Data->pon;\n");
  for(i=0;i<k2max;i++)
    fprintf(code,"  REALTYPE %s;\n",k2[i]);
  if (a==0)
    for(i=0;i<k1max;i++)
      fprintf(code,"  REALTYPE %s;\n",k1[i]);
  else
    for(i=0;i<k3max;i++)
      fprintf(code,"  REALTYPE %s;\n",k3[i]);
}

void define_localv(int a, int foo, int k1max, int k2max, int k3max, FILE *code)
{
  int i;
  
  for(i=0;i<k2max;i++)
    fprintf(code,"  %s = %.1lf*Data->oo2zn;\n",k2[i],(double)(i+1));
  if(a==0)
    for(i=0;i<k1max;i++)
      fprintf(code,"  %s = %.1lf*Data->oo2z;\n",k1[i],(double)(i+1));
  else
    for(i=0;i<k3max;i++)
      fprintf(code,"  %s = %.1lf*Data->oo2n;\n",k3[i],(double)(i+1));
  fprintf(code,"  U%d0 = Data->U[%d][0];\n", a*2, a*2);
  fprintf(code,"  U%d1 = Data->U[%d][1];\n", a*2, a*2);
  fprintf(code,"  U%d2 = Data->U[%d][2];\n", a*2, a*2);
  fprintf(code,"  U%d0 = Data->U[%d][0];\n", foo, foo);
  fprintf(code,"  U%d1 = Data->U[%d][1];\n", foo, foo);
  fprintf(code,"  U%d2 = Data->U[%d][2];\n\n", foo, foo);

}
