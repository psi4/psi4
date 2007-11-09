/*! \file 
    \ingroup (INT)
    \brief Enter brief description of file here 
*/
#include <math.h>
#include <stdio.h>
#include "build_libint.h"
#include <libint/constants.h>

extern FILE *outfile, *hrr_header;
extern LibintParams_t Params;

extern void punt(char *);

int emit_hrr_build()
{
  int new_am = Params.new_am;
  int max_class_size = Params.max_class_size;
  int am_to_inline = Params.max_am_to_inline_hrr_worker;

  FILE *code;
  int p,q,r,s;
  int ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz;
  int t0, t1, t2, t3, t4;
  int i,j,nj,i_i0,i_i1;
  int k,l,nl,k_i0,k_i1;
  int i0_step,i1_step;
  int a, b;
  int flag;
  int am_in[2];
  int am[2][3];
  int current_highest_am, to_inline;
  int xyz;
  int class_size;
  int split;
  int la, lb;
  int ld, lc, ld_max;
  int curr_count,curr_subfunction;
  int num_subfunctions, subbatch_length;
  int f;
  char code_name[20];
  char function_name[18];
  char **subfunction_name;
  

  for(lc=0;lc<=new_am;lc++) {
    ld_max = lc/2 + 1;
    if (ld_max > lc)
      ld_max = lc;
    for(ld=1;ld<=ld_max;ld++) {

      /*-----------------------
	HRR on centers C and D
       -----------------------*/

      am_in[0] = lc-ld;
      am_in[1] = ld;

      /* Is this function to be made inline */
      current_highest_am = (am_in[0] > am_in[1]) ? am_in[0] : am_in[1];
      to_inline = (current_highest_am <= am_to_inline) ? 1 : 0;
      
      class_size = ((am_in[0]+1)*(am_in[0]+2)*(am_in[1]+1)*(am_in[1]+2))/4;

      if (to_inline)
	fprintf(hrr_header, "#ifndef INLINE_HRR_WORKER\n");
      fprintf(hrr_header,"void hrr3_build_%c%c(const REALTYPE *, REALTYPE *, REALTYPE *, REALTYPE *, int);\n",
	      am_letter[am_in[0]],am_letter[am_in[1]]);
      if (to_inline)
	fprintf(hrr_header, "#endif\n");
      /* Decide if the routine has to be split into several routines producing "subbatches" */
      if (class_size > max_class_size) {
	split = 1;
	num_subfunctions = ceil((double)class_size/max_class_size);
	subbatch_length = 1 + class_size/num_subfunctions;
      }
      else {
	split = 0;
      }
      
      sprintf(function_name,"hrr3_build_%c%c",am_letter[am_in[0]],am_letter[am_in[1]]);
      if (split) {
	subfunction_name = (char **) malloc (num_subfunctions*sizeof(char *));
	for(i=0;i<num_subfunctions;i++) {
	  subfunction_name[i] = (char *) malloc(20*sizeof(char));
	  sprintf(subfunction_name[i],"_%s_%d",
		  function_name,i);
	}
      }
      
      sprintf(code_name,"%s.cc",function_name);
      code = fopen(code_name,"w");

      fprintf(code,"  /* These machine-generated functions compute a quartet of |%c%c) integrals */\n\n",
	      am_letter[am_in[0]],am_letter[am_in[1]]);
      fprintf(code,"#include \"libint.h\"\n\n");
      if (split) {
	for(i=0;i<num_subfunctions;i++) {
	  fprintf(code,"REALTYPE *%s(const REALTYPE *, REALTYPE *, REALTYPE *, REALTYPE *);\n",
		  subfunction_name[i]);
	}
	fprintf(code,"\n");
      }
      fprintf(code,"void %s(const REALTYPE *CD, REALTYPE *vp, REALTYPE *I0, REALTYPE *I1, int ab_num)\n{\n",
	      function_name);
      if (split == 1) {
	curr_subfunction = 0;
	curr_count = 0;
      }
      else {
	fprintf(code,"  const REALTYPE CD0 = CD[0];\n");
	fprintf(code,"  const REALTYPE CD1 = CD[1];\n");
	fprintf(code,"  const REALTYPE CD2 = CD[2];\n");
      }
      fprintf(code,"  int ab;\n\n");

      nl = (am_in[1]*(am_in[1]+1))/2;
      i0_step = (am_in[0]+2)*(am_in[0]+3)*nl/2;
      i1_step = (am_in[0]+1)*(am_in[0]+2)*nl/2;
      fprintf(code,"  for(ab=0;ab<ab_num;ab++) {\n");
      if (split == 1) {
	for(f=0;f<num_subfunctions;f++)
	  fprintf(code,"    vp = %s(CD, vp, I0, I1);\n",
		subfunction_name[f]);
	fprintf(code,"    I0 += %d;\n    I1 += %d;\n",i0_step,i1_step);
	fprintf(code,"  }\n}\n\n");

	fprintf(code,"REALTYPE *%s(const REALTYPE *CD, REALTYPE *vp, REALTYPE *I0, REALTYPE *I1)\n{\n",
	      subfunction_name[0]);
	fprintf(code,"  const REALTYPE CD0 = CD[0];\n");
	fprintf(code,"  const REALTYPE CD1 = CD[1];\n");
	fprintf(code,"  const REALTYPE CD2 = CD[2];\n\n");
      }

      for(p = 0; p <= am_in[0]; p++){
	am[0][0] = am_in[0] - p;
	for(q = 0; q <= p; q++){
	  am[0][1] = p - q;
	  am[0][2] = q;
	  
	  for(r = 0; r <= am_in[1]; r++){
	    am[1][0] = am_in[1] - r;
	    for(s = 0; s <= r; s++){
	      am[1][1] = r - s;
	      am[1][2] = s;

	      if (am[1][0]) /* build along x */
		xyz = 0;
	      else if (am[1][1]) /* build along y */
		xyz = 1;
	      else /*build along z */
		xyz = 2;

	      am[0][xyz] += 1;
	      am_in[0] += 1;
	      am[1][xyz] -= 1;
	      am_in[1] -= 1;
	      t0 = hash(am,am_in);
	      am[0][xyz] -= 1;
	      am_in[0] -= 1;
	      t1 = hash(am,am_in);
	      am[1][xyz] += 1;
	      am_in[1] += 1;
	      
	      fprintf(code, "    *(vp++) = I0[%d] + CD%d*I1[%d];\n",t0,xyz,t1);

	      curr_count++;
	      if (curr_count == subbatch_length && split == 1) {
		curr_count = 0;
		curr_subfunction++;
		fprintf(code,"  return vp;\n}\n\n");
		fprintf(code,"REALTYPE *%s(const REALTYPE *CD, REALTYPE *vp, REALTYPE *I0, REALTYPE *I1)\n{\n",
			subfunction_name[curr_subfunction]);
		fprintf(code,"  const REALTYPE CD0 = CD[0];\n");
		fprintf(code,"  const REALTYPE CD1 = CD[1];\n");
		fprintf(code,"  const REALTYPE CD2 = CD[2];\n\n");
	      }
	    }
	  }
	}
      }
      if (split == 0) {
	fprintf(code,"    I0 += %d;\n    I1 += %d;\n",i0_step,i1_step);
	fprintf(code,"  }\n}\n");
      }
      else {
	fprintf(code,"  return vp;\n}\n");
      }
      fclose(code);
      printf("Done with %s\n",code_name);
      
      
      /*-----------------------
	HRR on centers A and B
       -----------------------*/

      la = lc-ld;  lb = ld;
      am_in[0] = la;
      am_in[1] = lb;

      class_size = ((am_in[0]+1)*(am_in[0]+2)*(am_in[1]+1)*(am_in[1]+2))/4;

      if (to_inline)
	fprintf(hrr_header,"#ifndef INLINE_HRR_WORKER\n");
      fprintf(hrr_header,"void hrr1_build_%c%c(const REALTYPE *, REALTYPE *, REALTYPE *, REALTYPE *, int);\n",
	      am_letter[am_in[0]],am_letter[am_in[1]]);
      if (to_inline)
	fprintf(hrr_header,"#endif\n");
      /* Decide if the routine has to be split into several routines producing "subbatches" */
      if (class_size > max_class_size) {
	split = 1;
	num_subfunctions = ceil((double)class_size/max_class_size);
	subbatch_length = 1 + class_size/num_subfunctions;
      }
      else {
	split = 0;
      }
      
      sprintf(function_name,"hrr1_build_%c%c",am_letter[am_in[0]],am_letter[am_in[1]]);
      if (split) {
	subfunction_name = (char **) malloc (num_subfunctions*sizeof(char *));
	for(i=0;i<num_subfunctions;i++) {
	  subfunction_name[i] = (char *) malloc(20*sizeof(char));
	  sprintf(subfunction_name[i],"_%s_%d",
		  function_name,i);
	}
      }
      
      sprintf(code_name,"%s.cc",function_name);
      code = fopen(code_name,"w");
      fprintf(code,"  /* This machine-generated function computes a quartet of (%c%c| integrals */\n\n",
	      am_letter[am_in[0]],am_letter[am_in[1]]);
      fprintf(code,"#include \"libint.h\"\n\n");
      if (split) {
	for(i=0;i<num_subfunctions;i++) {
	  fprintf(code,"REALTYPE *%s(const REALTYPE *, REALTYPE *, REALTYPE *, REALTYPE *, int);\n",
		  subfunction_name[i]);
	}
	fprintf(code,"\n");
      }
      fprintf(code,"void %s(const REALTYPE *AB, REALTYPE *vp, REALTYPE *I0, REALTYPE *I1, int cd_num)\n{\n",
	      function_name);
      if (split == 1) {
	curr_subfunction = 0;
	curr_count = 0;
      }
      else {
	fprintf(code,"  const REALTYPE AB0 = AB[0];\n");
	fprintf(code,"  const REALTYPE AB1 = AB[1];\n");
	fprintf(code,"  const REALTYPE AB2 = AB[2];\n");
	fprintf(code,"  REALTYPE *i0, *i1;\n");
	fprintf(code,"  int cd;\n\n");
      }
      fprintf(code,"\n");

      if (split == 1) {
	for(f=0;f<num_subfunctions;f++)
	  fprintf(code,"  vp = %s(AB, vp, I0, I1, cd_num);\n",
		subfunction_name[f]);
	fprintf(code,"}\n\n");

	fprintf(code,"REALTYPE *%s(const REALTYPE *AB, REALTYPE *vp, REALTYPE *I0, REALTYPE *I1, int cd_num)\n{\n",
	      subfunction_name[0]);
	fprintf(code,"  const REALTYPE AB0 = AB[0];\n");
	fprintf(code,"  const REALTYPE AB1 = AB[1];\n");
	fprintf(code,"  const REALTYPE AB2 = AB[2];\n");
	fprintf(code,"  REALTYPE *i0, *i1;\n");
	fprintf(code,"  int cd;\n\n");
      }
      
      nj = (lb*(lb+1))/2;

      for(p = 0; p <= am_in[0]; p++){
	am[0][0] = am_in[0] - p;
	for(q = 0; q <= p; q++){
	  am[0][1] = p - q;
	  am[0][2] = q;
	  
	  for(r = 0; r <= am_in[1]; r++){
	    am[1][0] = am_in[1] - r;
	    for(s = 0; s <= r; s++){
	      am[1][1] = r - s;
	      am[1][2] = s;

	      if (am[1][0]) /* build along x */
		xyz = 0;
	      else if (am[1][1]) /* build along y */
		xyz = 1;
	      else /* build along z */
		xyz = 2;

	      am[0][xyz] += 1;
	      am_in[0] += 1;
	      am[1][xyz] -= 1;
	      am_in[1] -= 1;
	      t0 = hash(am,am_in);
	      am[0][xyz] -= 1;
	      am_in[0] -= 1;
	      t1 = hash(am,am_in);
	      am[1][xyz] += 1;
	      am_in[1] += 1;
	      
	      if (t0)
		fprintf(code,"  i0 = I0 + %d*cd_num;\n",t0);
	      else
		fprintf(code,"  i0 = I0;\n");
	      if (t1)
		fprintf(code,"  i1 = I1 + %d*cd_num;\n",t1);
	      else
		fprintf(code,"  i1 = I1;\n");

	      fprintf(code,"  for(cd=0;cd<cd_num;cd++)\n");
	      fprintf(code,"    *(vp++) = *(i0++) + AB%d*(*(i1++));\n",xyz);

	      curr_count++;
	      if (curr_count == subbatch_length && split == 1) {
		curr_count = 0;
		curr_subfunction++;
		fprintf(code,"  return vp;\n}\n\n");
		fprintf(code,"REALTYPE *%s(const REALTYPE *AB, REALTYPE *vp, REALTYPE *I0, REALTYPE *I1, int cd_num)\n{\n",
			subfunction_name[curr_subfunction]);
		fprintf(code,"  const REALTYPE AB0 = AB[0];\n");
		fprintf(code,"  const REALTYPE AB1 = AB[1];\n");
		fprintf(code,"  const REALTYPE AB2 = AB[2];\n");
		fprintf(code,"  REALTYPE *i0, *i1;\n");
		fprintf(code,"  int cd;\n\n");
	      }
	    }
	  }
	}
      }
      if (split == 1) {
	fprintf(code,"  return vp;\n");
      }
      fprintf(code,"}\n");
      fclose(code);
      printf("Done with %s\n",code_name);
    }
  }
}


