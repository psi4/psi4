/*! \file
    \ingroup R12
    \brief Enter brief description of file here 
*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "build_libr12.h"
#include <libint/constants.h>

extern FILE *hrr_header;
extern Libr12Params_t Params;

extern void punt(char *);

int emit_hrr_t_build()
{
  int new_am = Params.new_am;
  int max_class_size = Params.max_class_size;

  FILE *code;
  int p,q,r,s;
  int ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz;
  int t0, t1, t2, t3, t4;
  int i,j,nj,i_i0,i_i1;
  int k,l,nl,k_i0,k_i1;
  int cp1dm1_num,cdm1_num;
  int a, b;
  int flag;
  int am_in[2];
  int am[2][3];
  int xyz;
  int class_size;
  int split;
  int la, lb;
  int ld, lc, ld_max;
  int curr_count,curr_subfunction;
  int num_subfunctions, subbatch_length;
  int f;
  char code_name[21];
  char function_name[18];
  char **subfunction_name;
  
  for(lc=0;lc<=new_am;lc++) {
    ld_max = (lc+1)/2;
    for(ld=1;ld<=ld_max;ld++) {

      /*-----------------------
	HRR on centers C and D
       -----------------------*/

      am_in[0] = lc-ld;
      am_in[1] = ld;

      class_size = ((am_in[0]+1)*(am_in[0]+2)*(am_in[1]+1)*(am_in[1]+2))/4;
      nl = (am_in[1]*(am_in[1]+1))/2;
      cp1dm1_num = (am_in[0]+2)*(am_in[0]+3)*nl/2;
      cdm1_num = (am_in[0]+1)*(am_in[0]+2)*nl/2;

      /* Decide if the routine has to be split into several routines producing "subbatches" */
      if (class_size > max_class_size) {
	split = 1;
	num_subfunctions = ceil((double)class_size/max_class_size);
	subbatch_length = 1 + class_size/num_subfunctions;
      }
      else {
	split = 0;
      }
      
      sprintf(function_name,"t2hrr3_build_%c%c",am_letter[am_in[0]],am_letter[am_in[1]]);
      if (split) {
	subfunction_name = (char **) malloc (num_subfunctions*sizeof(char *));
	for(i=0;i<num_subfunctions;i++) {
	  subfunction_name[i] = (char *) malloc(22*sizeof(char));
	  sprintf(subfunction_name[i],"_%s_%d",
		  function_name,i);
	}
      }
      sprintf(code_name,"%s.cc",function_name);
      code = fopen(code_name,"w");

      /* include the function into the hrr_header.h */
      fprintf(hrr_header,"void %s(REALTYPE *, REALTYPE *, REALTYPE *, const REALTYPE *, const REALTYPE *, ",function_name);
      fprintf(hrr_header,"const REALTYPE *, const REALTYPE *, const REALTYPE *, int, int);\n");

      fprintf(code,"  /* This machine-generated function computes a quartet of [r12,T2]|%c%c) integrals */\n\n",
	      am_letter[am_in[0]],am_letter[am_in[1]]);
      fprintf(code,"#include<libint/libint.h>\n\n");
      if (split) {
	for(i=0;i<num_subfunctions;i++) {
	  fprintf(code,"REALTYPE *%s(REALTYPE *, REALTYPE *, REALTYPE *, const REALTYPE *, const REALTYPE *, ",
		  subfunction_name[i]);
	  fprintf(code,"const REALTYPE *, const REALTYPE *, const REALTYPE *, int, int);\n");
	}
	fprintf(code,"\n");
      }
      fprintf(code,"void %s(REALTYPE *CD, REALTYPE *AC, REALTYPE *vp, const REALTYPE *I0, const REALTYPE *I1, ",function_name);
      fprintf(code,"const REALTYPE *I2, const REALTYPE *I3, const REALTYPE *I4, int la, int lb)\n{\n");
      if (split == 1) {
	curr_subfunction = 0;
	curr_count = 0;
      }
      else {
	fprintf(code,"  const REALTYPE CD0 = CD[0];\n");
	fprintf(code,"  const REALTYPE CD1 = CD[1];\n");
	fprintf(code,"  const REALTYPE CD2 = CD[2];\n");
	fprintf(code,"  const REALTYPE AC0 = AC[0];\n");
	fprintf(code,"  const REALTYPE AC1 = AC[1];\n");
	fprintf(code,"  const REALTYPE AC2 = AC[2];\n");
	fprintf(code,"  static int io[] = { 1");
	for(i=1;i<=new_am+1;i++)
	  fprintf(code,", %d",(i+1)*(i+2)/2);
	fprintf(code,"};\n");
	fprintf(code,"  int bcdm1_num = io[lb]*%d;\n\n",cdm1_num);
      }
      fprintf(code,"  int pa, qa, b;\n");

      fprintf(code,"  for(pa=0;pa<=la;pa++)\n");
      fprintf(code,"    for(qa=0;qa<=pa;qa++)\n");
      fprintf(code,"      for(b=0;b<((lb+1)*(lb+2)/2);b++) {\n");

      if (split == 1) {
	for(f=0;f<num_subfunctions;f++)
	  fprintf(code,"    vp = %s(CD, AC, vp, I0, I1, I2, I3, I4, pa, lb);\n",
		subfunction_name[f]);
	fprintf(code,"        I0 += %d;\n",cp1dm1_num);
	fprintf(code,"        I1 += %d;\n",cdm1_num);
	fprintf(code,"        I2 += %d;\n",cp1dm1_num);
	fprintf(code,"        I3 += %d;\n",cdm1_num);
	fprintf(code,"        I4 += %d;\n",cdm1_num);
	fprintf(code,"  }\n}\n\n");

	fprintf(code,"REALTYPE *%s(REALTYPE *CD, REALTYPE *AC, REALTYPE *vp, const REALTYPE *I0, const REALTYPE *I1, ",
	      subfunction_name[0]);
	fprintf(code,"const REALTYPE *I2, const REALTYPE *I3, const REALTYPE *I4, int pa, int lb)\n{\n");
	fprintf(code,"  const REALTYPE CD0 = CD[0];\n");
	fprintf(code,"  const REALTYPE CD1 = CD[1];\n");
	fprintf(code,"  const REALTYPE CD2 = CD[2];\n");
	fprintf(code,"  const REALTYPE AC0 = AC[0];\n");
	fprintf(code,"  const REALTYPE AC1 = AC[1];\n");
	fprintf(code,"  const REALTYPE AC2 = AC[2];\n");
	fprintf(code,"  static int io[] = { 1");
	for(i=1;i<=new_am+1;i++)
	  fprintf(code,", %d",(i+1)*(i+2)/2);
	fprintf(code,"};\n");
	fprintf(code,"  int bcdm1_num = io[lb]*%d;\n\n",cdm1_num);
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
	      t2 = t0;
	      t3 = t1;
	      t4 = t1;
	      
	      fprintf(code,"        *(vp++) = I0[%d] + CD%d*I1[%d] + I2[%d] - AC%d*I4[%d] - ",t0,xyz,t1,t2,xyz,t4);
	      if (xyz == 0)
		fprintf(code,"I3[%d];\n",t3);
	      else
		fprintf(code,"I3[(pa+%d)*bcdm1_num+%d];\n",xyz,t3);

	      curr_count++;
	      if (curr_count == subbatch_length && split == 1) {
		curr_count = 0;
		curr_subfunction++;
		fprintf(code,"  return vp;\n}\n\n");
		fprintf(code,"REALTYPE *%s(REALTYPE *CD, REALTYPE *AC, REALTYPE *vp, const REALTYPE *I0, const REALTYPE *I1, ",
			subfunction_name[curr_subfunction]);
		fprintf(code,"const REALTYPE *I2, const REALTYPE *I3, const REALTYPE *I4, int pa, int lb)\n{\n");
		fprintf(code,"  const REALTYPE CD0 = CD[0];\n");
		fprintf(code,"  const REALTYPE CD1 = CD[1];\n");
		fprintf(code,"  const REALTYPE CD2 = CD[2];\n");
		fprintf(code,"  const REALTYPE AC0 = AC[0];\n");
		fprintf(code,"  const REALTYPE AC1 = AC[1];\n");
		fprintf(code,"  const REALTYPE AC2 = AC[2];\n");
		fprintf(code,"  static int io[] = { 1");
		for(i=1;i<=new_am+1;i++)
		  fprintf(code,", %d",(i+1)*(i+2)/2);
		fprintf(code,"};\n");
		fprintf(code,"  int bcdm1_num = io[lb]*%d;\n\n",cdm1_num);
	      }
	    }
	  }
	}
      }
      if (split == 0) {
	fprintf(code,"        I0 += %d;\n",cp1dm1_num);
	fprintf(code,"        I1 += %d;\n",cdm1_num);
	fprintf(code,"        I2 += %d;\n",cp1dm1_num);
	fprintf(code,"        I3 += %d;\n",cdm1_num);
	fprintf(code,"        I4 += %d;\n",cdm1_num);
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

      /*--- Whether to split has been decided already ---*/
      
      sprintf(function_name,"t1hrr1_build_%c%c",am_letter[am_in[0]],am_letter[am_in[1]]);
      if (split) {
	subfunction_name = (char **) malloc (num_subfunctions*sizeof(char *));
	for(i=0;i<num_subfunctions;i++) {
	  subfunction_name[i] = (char *) malloc(22*sizeof(char));
	  sprintf(subfunction_name[i],"_%s_%d",
		  function_name,i);
	}
      }
      sprintf(code_name,"t1hrr1_build_%c%c.cc",am_letter[am_in[0]],am_letter[am_in[1]]);
      code = fopen(code_name,"w");

      /* include the function into the hrr_header.h */
      fprintf(hrr_header,"void %s(REALTYPE *, REALTYPE *, REALTYPE *, const REALTYPE *, const REALTYPE *, ",function_name);
      fprintf(hrr_header,"const REALTYPE *, const REALTYPE *, const REALTYPE *, int, int);\n");
      
      fprintf(code,"  /* This machine-generated function computes a quartet of (%c%c|[r12,T1] integrals */\n\n",
	      am_letter[am_in[0]],am_letter[am_in[1]]);
      fprintf(code,"#include<libint/libint.h>\n\n");
      
      if (split) {
	for(i=0;i<num_subfunctions;i++) {
	  fprintf(code,"REALTYPE *%s(REALTYPE *, REALTYPE *, REALTYPE *, const REALTYPE *, const REALTYPE *, ",
		  subfunction_name[i]);
	  fprintf(code,"const REALTYPE *, const REALTYPE *, const REALTYPE *, int, int);\n");
	}
	fprintf(code,"\n");
      }
      fprintf(code,"void %s(REALTYPE *AB, REALTYPE *AC, REALTYPE *vp, const REALTYPE *I0, const REALTYPE *I1, ",
	      function_name);
      fprintf(code,"const REALTYPE *I2, const REALTYPE *I3, const REALTYPE *I4, int lc, int ld)\n{\n");
      if (split == 1) {
	curr_subfunction = 0;
	curr_count = 0;
      }
      else {
	fprintf(code,"  int cd, cd_num, c_num, cp1_num, d_num;\n");
	fprintf(code,"  int pc, qc, d, ind_c, ind_cp1d;\n");
	fprintf(code,"  const REALTYPE *i0, *i1, *i2, *i3, *i4;\n");
	fprintf(code,"  const REALTYPE AB0 = AB[0];\n");
	fprintf(code,"  const REALTYPE AB1 = AB[1];\n");
	fprintf(code,"  const REALTYPE AB2 = AB[2];\n");
	fprintf(code,"  const REALTYPE AC0 = AC[0];\n");
	fprintf(code,"  const REALTYPE AC1 = AC[1];\n");
	fprintf(code,"  const REALTYPE AC2 = AC[2];\n");
	fprintf(code,"  static int io[] = { 1");
	for(i=1;i<=new_am+1;i++)
	  fprintf(code,", %d",(i+1)*(i+2)/2);
	fprintf(code,"};\n\n");
	fprintf(code,"  c_num = io[lc];\n");
	fprintf(code,"  cp1_num = io[lc+1];\n");
	fprintf(code,"  d_num = io[ld];\n\n");
      }

      if (split == 1) {
	for(f=0;f<num_subfunctions;f++)
	  fprintf(code,"  vp = %s(AB, AC, vp, I0, I1, I2, I3, I4, lc, ld);\n",
		subfunction_name[f]);
	fprintf(code,"}\n\n");

	fprintf(code,"REALTYPE *%s(REALTYPE *AB, REALTYPE *AC, REALTYPE *vp, const REALTYPE *I0, const REALTYPE *I1, ",
		subfunction_name[0]);
	fprintf(code,"const REALTYPE *I2, const REALTYPE *I3, const REALTYPE *I4, int lc, int ld)\n{\n");
	fprintf(code,"  int cd, cd_num, c_num, cp1_num, d_num;\n");
	fprintf(code,"  int pc, qc, d, ind_c, ind_cp1d;\n");
	fprintf(code,"  const REALTYPE *i0, *i1, *i2, *i3, *i4;\n");
	fprintf(code,"  const REALTYPE AB0 = AB[0];\n");
	fprintf(code,"  const REALTYPE AB1 = AB[1];\n");
	fprintf(code,"  const REALTYPE AB2 = AB[2];\n");
	fprintf(code,"  const REALTYPE AC0 = AC[0];\n");
	fprintf(code,"  const REALTYPE AC1 = AC[1];\n");
	fprintf(code,"  const REALTYPE AC2 = AC[2];\n");
	fprintf(code,"  static int io[] = { 1");
	for(i=1;i<=new_am+1;i++)
	  fprintf(code,", %d",(i+1)*(i+2)/2);
	fprintf(code,"};\n\n");
	fprintf(code,"  c_num = io[lc];\n");
	fprintf(code,"  cp1_num = io[lc+1];\n");
	fprintf(code,"  d_num = io[ld];\n\n");
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
	      t2 = t0;
	      t3 = t1;
	      t4 = t1;
	      
	      if (t0) {
		fprintf(code,"  i0 = I0 + %d*(c_num*d_num);\n",t0);
		fprintf(code,"  i2 = I2 + %d*(c_num*d_num);\n",t2);
	      }
	      else {
		fprintf(code,"  i0 = I0;\n");
		fprintf(code,"  i2 = I2;\n");
	      }
	      if (t1) {
		fprintf(code,"  i1 = I1 + %d*(c_num*d_num);\n",t1);
		fprintf(code,"  i4 = I4 + %d*(c_num*d_num);\n",t4);
		fprintf(code,"  i3 = I3 + %d*(cp1_num*d_num);\n",t3);
	      }
	      else {
		fprintf(code,"  i1 = I1;\n");
		fprintf(code,"  i3 = I3;\n");
		fprintf(code,"  i4 = I4;\n");
	      }

	      if (xyz == 0) {
		fprintf(code,"  for(cd=0;cd<c_num*d_num;cd++)\n");
		fprintf(code,"    *(vp++) = *(i0++) + AB0*(*(i1++)) + *(i2++) + AC0*(*(i4++)) - *(i3++);\n\n");
	      }
	      else {
		fprintf(code,"  ind_c = %d;\n",xyz);
		fprintf(code,"  for(pc=0;pc<=lc;pc++) {\n");
		fprintf(code,"    for(qc=0;qc<=pc;qc++) {\n");
		fprintf(code,"      ind_cp1d = (ind_c + pc)*d_num;\n");
		fprintf(code,"      for(d=0;d<d_num;d++) {\n");
		fprintf(code,"        *(vp++) = *(i0++) + AB%d*(*(i1++)) + *(i2++) + AC%d*(*(i4++)) - i3[ind_cp1d];\n",
			xyz,xyz);
		fprintf(code,"        ind_cp1d++;\n");
		fprintf(code,"      }\n");
		fprintf(code,"      ind_c++;\n");
		fprintf(code,"    }\n");
		fprintf(code,"  }\n\n");
	      }

	      curr_count++;
	      if (curr_count == subbatch_length && split == 1) {
		curr_count = 0;
		curr_subfunction++;
		fprintf(code,"  return vp;\n}\n\n");
		fprintf(code,"REALTYPE *%s(REALTYPE *AB, REALTYPE *AC, REALTYPE *vp, const REALTYPE *I0, const REALTYPE *I1, ",
			subfunction_name[curr_subfunction]);
		fprintf(code,"const REALTYPE *I2, const REALTYPE *I3, const REALTYPE *I4, int lc, int ld)\n{\n");
		fprintf(code,"  int cd, cd_num, c_num, cp1_num, d_num;\n");
		fprintf(code,"  int pc, qc, d, ind_c, ind_cp1d;\n");
		fprintf(code,"  const REALTYPE *i0, *i1, *i2, *i3, *i4;\n");
		fprintf(code,"  const REALTYPE AB0 = AB[0];\n");
		fprintf(code,"  const REALTYPE AB1 = AB[1];\n");
		fprintf(code,"  const REALTYPE AB2 = AB[2];\n");
		fprintf(code,"  const REALTYPE AC0 = AC[0];\n");
		fprintf(code,"  const REALTYPE AC1 = AC[1];\n");
		fprintf(code,"  const REALTYPE AC2 = AC[2];\n");
		fprintf(code,"  static int io[] = { 1");
		for(i=1;i<=new_am+1;i++)
		  fprintf(code,", %d",(i+1)*(i+2)/2);
		fprintf(code,"};\n\n");
		fprintf(code,"  c_num = io[lc];\n");
		fprintf(code,"  cp1_num = io[lc+1];\n");
		fprintf(code,"  d_num = io[ld];\n\n");
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


