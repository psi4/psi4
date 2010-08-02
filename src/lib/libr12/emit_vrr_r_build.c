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
static void declare_localv(int dec_C, int k1max, int k2max, int k3max, FILE *code);
static void define_localv(int dec_C, int k1max, int k2max, int k3max, FILE *code);

static char **k1, **k2, **k3;

void emit_vrr_r_build()
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
  int t1, t2, t3, t4;
  int class_size;
  int type;
  int max1 = 0;
  int max2 = 0;
  int la, lc, lc_min, lc_max;
  int k1max, k2max, k3max;
  int split,num_subfunctions,subbatch_length;
  int curr_count,curr_subfunction;
  static const char *k1_suff = "o2z";
  static const char *k2_suff = "o4zn";
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
      fprintf(outfile,"  AM_a = %c  AM_c = %c\n",am_letter[la],am_letter[lc]);
      am_in[0] = la;
      am_in[1] = lc;
      if (la == 0) {
	dec_C = 1;
	k2max = lc - 1;
	k3max = lc - 1;
      }
      else {
	dec_C = 0;
	k2max = (la-1 > lc) ? la-1 : lc;
	k1max = la - 1;
      }

      class_size = ((am_in[dec_C]+1)*(am_in[dec_C]+2)*(am_in[dec_C^1]+1)*(am_in[dec_C^1]+2))/4;

      fprintf(vrr_header,"#define _R_BUILD_%c0%c0(Data,vp,i0,i1,i2,i3,i4,i5) {",am_letter[la],am_letter[lc]);
      fprintf(outfile,"  # of integrals in the (%cs|%cs) class - %d\n",am_letter[la],am_letter[lc],class_size);
      /* Decide if the routine has to be split into several routines producing "subbatches" */
      if (class_size > max_class_size) {
	split = 1;
	num_subfunctions = ceil((double)class_size/max_class_size);
	subbatch_length = 1 + class_size/num_subfunctions;
	fprintf(outfile,"  Each function for this quartet split into %d sub_functions\n\n",num_subfunctions);
	fprintf(vrr_header," tmp = _r_build_%c0%c0_0(Data,vp,i0,i1,i2,i3,i4,i5); \\\n",am_letter[la],am_letter[lc]);
	for(f=1;f<num_subfunctions;f++)
	  fprintf(vrr_header," tmp = _r_build_%c0%c0_%d(Data,tmp,i0,i1,i2,i3,i4,i5); \\\n",am_letter[la],am_letter[lc],f);
	fprintf(vrr_header,"}\n");
	for(f=0;f<num_subfunctions;f++)
	  fprintf(vrr_header, 
	  "REALTYPE *_r_build_%c0%c0_%d(prim_data *, REALTYPE *, const REALTYPE *, const REALTYPE *, const REALTYPE *, const REALTYPE *, const REALTYPE *, const REALTYPE *);\n",
		  am_letter[la],am_letter[lc],f);
      }
      else {
	split = 0;
	fprintf(vrr_header," _r_build_%c0%c0(Data,vp,i0,i1,i2,i3,i4,i5);}\n",am_letter[la],am_letter[lc]);
	fprintf(vrr_header,
	"void _r_build_%c0%c0(prim_data *, REALTYPE *, const REALTYPE *, const REALTYPE *, const REALTYPE *, const REALTYPE *, const REALTYPE *, const REALTYPE *);\n",
		am_letter[la],am_letter[lc]);
      }

      sprintf(function_name,"r_build_%c0%c0",am_letter[la],am_letter[lc]);
      sprintf(code_name,"r_build_%c0%c0.cc",am_letter[la],am_letter[lc]);
      code = fopen(code_name,"w");

      /*target,I2[]
        |I0[],I3[]
        |    |I1[],I4[]
        |    |    |   I5[]
        |    |    |    |     */
        t1 = t2 = t3 = t4 = 0;

      /* print local variable declarations */

      fprintf(code,"  /* These machine-generated functions compute a quartet of (%cs||%cs) integrals */\n\n",
	      am_letter[la],am_letter[lc]);
      if (split) {
	subfunction_name = (char **) malloc (num_subfunctions*sizeof(char *));
	for(i=0;i<num_subfunctions;i++) {
	  subfunction_name[i] = (char *) malloc(20*sizeof(char));
	  sprintf(subfunction_name[i],"_r_build_%c0%c0_%d",am_letter[la],am_letter[lc],i);
	}
      }

      fprintf(code,"#include <libint/libint.h>\n");
      fprintf(code,"#include \"libr12.h\"\n\n");
      if (split == 1) {
	curr_subfunction = 0;
	curr_count = 0;
	fprintf(code,
	"REALTYPE *%s(prim_data *Data, REALTYPE *vp, const REALTYPE *I0, const REALTYPE *I1, const REALTYPE *I2, const REALTYPE *I3, const REALTYPE *I4, const REALTYPE *I5)\n{\n",
		subfunction_name[0]);
      }
      else
	fprintf(code,
	"void _%s(prim_data *Data, REALTYPE *vp, const REALTYPE *I0, const REALTYPE *I1, const REALTYPE *I2, const REALTYPE *I3, const REALTYPE *I4, const REALTYPE *I5)\n{\n",function_name);
      declare_localv(dec_C,k1max,k2max,k3max,code);
      define_localv(dec_C,k1max,k2max,k3max,code);
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

	      if(am[dec_C][2]) xyz = 2;
	      if(am[dec_C][1]) xyz = 1;
	      if(am[dec_C][0]) xyz = 0;

	      /*-----------------------------
		Add (a-10||c0) and (a-10|c0)
	       -----------------------------*/
	      am[dec_C][xyz] = am[dec_C][xyz] - 1;
	      am_in[dec_C] = am_in[dec_C] - 1;
	      t2 = hash(am,am_in);
	      fprintf(code, "*(vp++) = U%d%d*I0[%d] - U1%d*I3[%d]",
		      dec_C*2, xyz, t2, xyz, t2);
	      
	      /*------------
		Add (a0|c0)
	       ------------*/
	      fprintf(code, " + loo2p*I2[%d]", t1);

	      /*-----------------------------------------
		Add (a-20||c0) and (a-20|c0) if possible
	       -----------------------------------------*/
	      if(am[dec_C][xyz]){
		am[dec_C][xyz] = am[dec_C][xyz] - 1;
		am_in[dec_C] = am_in[dec_C] - 1;
		t3 = hash(am,am_in);
		fprintf(code, "\n        + (%s)*I1[%d] - (%s)*I4[%d]", 
			(dec_C==0 ? k1[am[dec_C][xyz]] : k3[am[dec_C][xyz]]), 
			t3, (k2[am[dec_C][xyz]]), t3);
		am[dec_C][xyz] = am[dec_C][xyz] + 1;
		am_in[dec_C] = am_in[dec_C] + 1;
	      }

	      /*----------------------------
		Add (a-10|c-10) if possible
	       ----------------------------*/
	      if(am[dec_C^1][xyz]){
		am[dec_C^1][xyz] = am[dec_C^1][xyz] - 1;
		am_in[dec_C^1] = am_in[dec_C^1] - 1;
		t4 = hash(am,am_in);
		fprintf(code, " - (%s)*I5[%d]", k2[am[dec_C^1][xyz]], t4);
		am[dec_C^1][xyz] = am[dec_C^1][xyz] + 1;
		am_in[dec_C^1] = am_in[dec_C^1] + 1;
	      }
	      fprintf(code, ";\n");
	      am[dec_C][xyz] = am[dec_C][xyz] + 1;
	      am_in[dec_C] = am_in[dec_C] + 1;
		
	      t1++;
	      curr_count++;
	      if (curr_count == subbatch_length && split == 1) {
		curr_count = 0;
		curr_subfunction++;
		fprintf(code,"return vp;\n}\n\n");
		fprintf(code,
		"REALTYPE *%s(prim_data *Data, REALTYPE *vp, const REALTYPE *I0, const REALTYPE *I1, const REALTYPE *I2, const REALTYPE *I3, const REALTYPE *I4, const REALTYPE *I5)\n{\n",
			subfunction_name[curr_subfunction]);
		declare_localv(dec_C,k1max,k2max,k3max,code);
		define_localv(dec_C,k1max,k2max,k3max,code);
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


void declare_localv(int dec_C, int k1max, int k2max, int k3max, FILE *code)
{
  int i;

  fprintf(code,"  REALTYPE U00, U01, U02, U10, U11, U12, U20, U21, U22;\n");
  fprintf(code,"  REALTYPE U30, U31, U32, U40, U41, U42, U50, U51, U52;\n");
  fprintf(code,"  REALTYPE loo2p = Data->oo2p;\n");
  for(i=0;i<k2max;i++)
    fprintf(code,"  REALTYPE %s;\n",k2[i]);
  if (dec_C == 0)
    for(i=0;i<k1max;i++)
      fprintf(code,"  REALTYPE %s;\n",k1[i]);
  else
    for(i=0;i<k3max;i++)
      fprintf(code,"  REALTYPE %s;\n",k3[i]);
}

void define_localv(int dec_C, int k1max, int k2max, int k3max, FILE *code)
{
  int i;
  
  for(i=0;i<k2max;i++)
    fprintf(code,"  %s = %.1lf*Data->oo2z*Data->oo2n;\n",k2[i],(double)(i+1));
  if(dec_C == 0)
    for(i=0;i<k1max;i++)
      fprintf(code,"  %s = %.1lf*Data->oo2z;\n",k1[i],(double)(i+1));
  else
    for(i=0;i<k3max;i++)
      fprintf(code,"  %s = %.1lf*Data->oo2n;\n",k3[i],(double)(i+1));
  fprintf(code,"  U%d0 = Data->U[%d][0];\n", dec_C*2, dec_C*2);
  fprintf(code,"  U%d1 = Data->U[%d][1];\n", dec_C*2, dec_C*2);
  fprintf(code,"  U%d2 = Data->U[%d][2];\n", dec_C*2, dec_C*2);
  if (dec_C == 0) {
      /*                     PA_i      /    2\eta  +    QA_i      /    2\zeta    */
    fprintf(code,"  U10 = Data->U[0][0]*Data->oo2n + Data->U[1][0]*Data->oo2z;\n");
    fprintf(code,"  U11 = Data->U[0][1]*Data->oo2n + Data->U[1][1]*Data->oo2z;\n");
    fprintf(code,"  U12 = Data->U[0][2]*Data->oo2n + Data->U[1][2]*Data->oo2z;\n");
  }
  else {
      /*                     QC_i      /    2\zeta  +   PC_i      /    2\eta    */
    fprintf(code,"  U10 = Data->U[2][0]*Data->oo2z + Data->U[3][0]*Data->oo2n;\n");
    fprintf(code,"  U11 = Data->U[2][1]*Data->oo2z + Data->U[3][1]*Data->oo2n;\n");
    fprintf(code,"  U12 = Data->U[2][2]*Data->oo2z + Data->U[3][2]*Data->oo2n;\n");
  }

  return;
}
