/*! \file
    \ingroup R12
    \brief Enter brief description of file here 
*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "build_libr12.h"
#include <libint/constants.h>

extern FILE *outfile, *vrr_header;
extern Libr12Params_t Params;

extern void punt(char *);
static void declare_localv(int la, FILE *code);
static void define_localv(int la, FILE *code);

static char **k1;

int emit_vrr_t1_build()
{
  int old_am = Params.old_am;
  int new_am = Params.opt_am;
  int max_class_size = Params.max_class_size;

  FILE *code;
  int i, j, k, l, f;
  int dec_C;             /* Decrease AM on C */
  int xyz;               /* Cartesian direction along which to decrease AM */
  int flag;
  int am[2][3];
  int am_in[2];
  int nflip = 0;
  int t0, t1, t2, t3, t4;
  int class_size;
  int type;
  int max1 = 0;
  int max2 = 0;
  int la, lc, lc_min, lc_max;
  int k1max;
  int split,num_subfunctions,subbatch_length;
  int curr_count,curr_subfunction;
  static const char *k1_suff = "zboz";
  char *code_name;
  char *function_name;
  char **subfunction_name;


  k1 = (char **) malloc(new_am*sizeof(char *));
  for(i=1;i<=new_am;i++) {
    j = strlen(number[i]);
    k1[i-1] = (char*) malloc((5+j)*sizeof(char));
    strcpy(k1[i-1],number[i]);
    strcat(k1[i-1],k1_suff);
  }
  code_name = (char *) malloc(sizeof(char)*21);
  function_name = (char *) malloc(sizeof(char)*18);

  for(la=0;la<=new_am;la++) {
    lc_min = (la >= old_am + 1) ? 0 : old_am + 1;
    if (la == 0 && old_am == 0)
      lc_min = 0;
    lc_max = new_am;
    for(lc=lc_min;lc<=lc_max;lc++) {
      am_in[0] = la;
      am_in[1] = lc;

      class_size = ((am_in[0]+1)*(am_in[0]+2)*(am_in[1]+1)*(am_in[1]+2))/4;

      fprintf(vrr_header,"#define _T1_BUILD_%c0%c0(Data,ShellData,vp,i0,i1,i2,i3,i4) {",am_letter[la],am_letter[lc]);
      fprintf(outfile,"  # of integrals in the (%cs|%cs) class - %d\n",am_letter[la],am_letter[lc],class_size);
      /* Decide if the routine has to be split into several routines producing "subbatches" */
      if (class_size > max_class_size) {
	split = 1;
	num_subfunctions = ceil((double)class_size/max_class_size);
	subbatch_length = 1 + class_size/num_subfunctions;
	fprintf(outfile,"  Each function for this quartet split into %d sub_functions\n\n",num_subfunctions);
	fprintf(vrr_header," tmp = _t1_build_%c0%c0_0(Data,ShellData,vp,i0,i1,i2,i3,i4); \\\n",am_letter[la],am_letter[lc]);
	for(f=1;f<num_subfunctions;f++)
	  fprintf(vrr_header," tmp = _t1_build_%c0%c0_%d(Data,ShellData,tmp,i0,i1,i2,i3,i4); \\\n",am_letter[la],am_letter[lc],f);
	fprintf(vrr_header,"}\n");
	for(f=0;f<num_subfunctions;f++)
	  fprintf(vrr_header, "REALTYPE *_t1_build_%c0%c0_%d(prim_data *, contr_data *, REALTYPE *, const REALTYPE *, const REALTYPE *, const REALTYPE *, const REALTYPE *, const REALTYPE *);\n",
		  am_letter[la],am_letter[lc],f);
      }
      else {
	split = 0;
	fprintf(vrr_header," _t1_build_%c0%c0(Data,ShellData,vp,i0,i1,i2,i3,i4);}\n",am_letter[la],am_letter[lc]);
	fprintf(vrr_header,"void _t1_build_%c0%c0(prim_data *, contr_data *, REALTYPE *, const REALTYPE *, const REALTYPE *, const REALTYPE *, const REALTYPE *, const REALTYPE *);\n",
		am_letter[la],am_letter[lc]);
      }

      sprintf(function_name,"t1_build_%c0%c0",am_letter[la],am_letter[lc]);
      sprintf(code_name,"t1_build_%c0%c0.cc",am_letter[la],am_letter[lc]);
      code = fopen(code_name,"w");

  /* target,I0[]
        |  I1[]
        |    |   I2[]
        |    |    |   I3[]
        |    |    |    |   I4[]
        |    |    |    |    |    */
        t0 = t1 = t2 = t3 = t4 = 0;

      /* print local variable declarations */

      fprintf(code,"  /* These machine-generated functions compute a quartet of (%cs|[r12,T1]|%cs) integrals */\n\n",
	      am_letter[la],am_letter[lc]);
      if (split) {
	subfunction_name = (char **) malloc (num_subfunctions*sizeof(char *));
	for(i=0;i<num_subfunctions;i++) {
	  subfunction_name[i] = (char *) malloc(20*sizeof(char));
	  sprintf(subfunction_name[i],"_t1_build_%c0%c0_%d",am_letter[la],am_letter[lc],i);
	}
      }

      fprintf(code,"#include <libint/libint.h>\n");
      fprintf(code,"#include \"libr12.h\"\n\n");
      if (split == 1) {
	curr_subfunction = 0;
	curr_count = 0;
	fprintf(code,"REALTYPE *%s(prim_data *Data, contr_data *ShellData, REALTYPE *vp, const REALTYPE *I0, const REALTYPE *I1, const REALTYPE *I2, const REALTYPE *I3, const REALTYPE *I4)\n{\n",
		subfunction_name[0]);
      }
      else
	fprintf(code,"void _%s(prim_data *Data, contr_data *ShellData, REALTYPE *vp, const REALTYPE *I0, const REALTYPE *I1, const REALTYPE *I2, const REALTYPE *I3, const REALTYPE *I4)\n{\n",function_name);
      declare_localv(la,code);
      define_localv(la,code);
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

	      /*------------
		Add (a0|c0)
	       ------------*/
	      fprintf(code, "*(vp++) = U0*I0[%d]", t0);

	      /*-----------------------------
		Add (a+10|c0) and (a0|c+10)
	       -----------------------------*/
	      for(xyz=0;xyz<3;xyz++) {
		am[0][xyz] = am[0][xyz] + 1;
		am_in[0] = am_in[0] + 1;
		t1 = hash(am,am_in);
		am[0][xyz] = am[0][xyz] - 1;
		am_in[0] = am_in[0] - 1;
		
		am[1][xyz] = am[1][xyz] + 1;
		am_in[1] = am_in[1] + 1;
		t2 = hash(am,am_in);
		am[1][xyz] = am[1][xyz] - 1;
		am_in[1] = am_in[1] - 1;
		fprintf(code, "\n          - U1%d*(I1[%d] - I2[%d])",
			xyz, t1, t2);
	      }
	      
	      /*-----------------------------------------
		Add (a-10|c+10) and (a-10|c0) if possible
	       -----------------------------------------*/
	      for(xyz=0;xyz<3;xyz++) {
		if(am[0][xyz]){
		  am[0][xyz] = am[0][xyz] - 1;
		  am_in[0] = am_in[0] - 1;
		  am[1][xyz] = am[1][xyz] + 1;
		  am_in[1] = am_in[1] + 1;
		  t3 = hash(am,am_in);
		  am[1][xyz] = am[1][xyz] - 1;
		  am_in[1] = am_in[1] - 1;
		  t4 = hash(am,am_in);
		  am[0][xyz] = am[0][xyz] + 1;
		  am_in[0] = am_in[0] + 1;

		  fprintf(code, "\n          + (%s)*(I3[%d] - AC%d*I4[%d])", 
			  k1[am[0][xyz]-1], t3, xyz, t4);
		}
	      }
	      fprintf(code,";\n");

	      t0++;
	      curr_count++;
	      if (curr_count == subbatch_length && split == 1) {
		curr_count = 0;
		curr_subfunction++;
		fprintf(code,"return vp;\n}\n\n");
		fprintf(code,"REALTYPE *%s(prim_data *Data, contr_data *ShellData, REALTYPE *vp, const REALTYPE *I0, const REALTYPE *I1, const REALTYPE *I2, const REALTYPE *I3, const REALTYPE *I4)\n{\n",
			subfunction_name[curr_subfunction]);
		declare_localv(la,code);
		define_localv(la,code);
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


void declare_localv(int la, FILE *code)
{
  int i;

  fprintf(code,"  REALTYPE U0, U10, U11, U12;\n");
  fprintf(code,"  REALTYPE AC0, AC1, AC2;\n");
  for(i=0;i<la;i++)
    fprintf(code,"  REALTYPE %s;\n",k1[i]);

  return;

}

void define_localv(int la, FILE *code)
{
  int i;

  fprintf(code,"  AC0 = ShellData->AC[0];\n");
  fprintf(code,"  AC1 = ShellData->AC[1];\n");
  fprintf(code,"  AC2 = ShellData->AC[2];\n");
  for(i=0;i<la;i++)
    fprintf(code,"  %s = %.1lf*(Data->twozeta_b*Data->oo2z);\n",k1[i],(double)(i+1));
  fprintf(code,"  U10 = ShellData->AB[0]*(Data->twozeta_a*Data->twozeta_b*Data->oo2z);\n");
  fprintf(code,"  U11 = ShellData->AB[1]*(Data->twozeta_a*Data->twozeta_b*Data->oo2z);\n");
  fprintf(code,"  U12 = ShellData->AB[2]*(Data->twozeta_a*Data->twozeta_b*Data->oo2z);\n");
  fprintf(code,"  U0  = (Data->twozeta_a - Data->twozeta_b*(Data->twozeta_a*ShellData->ABdotAC + %lf))*Data->oo2z;\n",(double)(la+1));

  return;
}
