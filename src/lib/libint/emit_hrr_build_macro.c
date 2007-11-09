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

int emit_hrr_build_macro()
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
      if (!to_inline)
	continue;
      
      class_size = ((am_in[0]+1)*(am_in[0]+2)*(am_in[1]+1)*(am_in[1]+2))/4;

      /* If the routine has to be split into several - user probably doesn't know what he/she is doing */
      if (class_size > max_class_size)
	punt("MAX_CLASS_SIZE is too small for the given inlining threshold");
      else {
	split = 0;
      }
      
      sprintf(function_name,"hrr3_build_%c%c",am_letter[am_in[0]],am_letter[am_in[1]]);
      sprintf(code_name,"%s.h",function_name);
      code = fopen(code_name,"w");

      fprintf(code,"#ifndef _libint_%s\n",function_name);
      fprintf(code,"#define _libint_%s\n",function_name);
      fprintf(code,"  /* These machine-generated functions compute a quartet of |%c%c) integrals */\n\n",
	      am_letter[am_in[0]],am_letter[am_in[1]]);
      fprintf(code,"#define %s(CD, vp, I0, I1, ab_num)\\\n{\\\n",
	      function_name);
      fprintf(code,"  const REALTYPE CD0 = CD[0];\\\n");
      fprintf(code,"  const REALTYPE CD1 = CD[1];\\\n");
      fprintf(code,"  const REALTYPE CD2 = CD[2];\\\n");
      fprintf(code,"  int ab;\\\n");
      fprintf(code,"  REALTYPE *target = (vp);\\\n");
      fprintf(code,"  REALTYPE *i0 = (I0);\\\n");
      fprintf(code,"  REALTYPE *i1 = (I1);\\\n\\\n");

      nl = (am_in[1]*(am_in[1]+1))/2;
      i0_step = (am_in[0]+2)*(am_in[0]+3)*nl/2;
      i1_step = (am_in[0]+1)*(am_in[0]+2)*nl/2;
      fprintf(code,"  for(ab=0;ab<ab_num;ab++) {\\\n");

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
	      
	      fprintf(code, "    *(target++) = i0[%d] + CD%d*i1[%d];\\\n",t0,xyz,t1);

	      curr_count++;
	    }
	  }
	}
      }
      fprintf(code,"    i0 += %d;\\\n    i1 += %d;\\\n",i0_step,i1_step);
      fprintf(code,"  }\\\n}\n");
      fprintf(code,"\n#endif\n"); /* end of #ifndef _libint_.... */
      fclose(code);
      printf("Done with %s\n",code_name);
      
      
      /*-----------------------
	HRR on centers A and B
       -----------------------*/

      la = lc-ld;  lb = ld;
      am_in[0] = la;
      am_in[1] = lb;

      class_size = ((am_in[0]+1)*(am_in[0]+2)*(am_in[1]+1)*(am_in[1]+2))/4;

      sprintf(function_name,"hrr1_build_%c%c",am_letter[am_in[0]],am_letter[am_in[1]]);
      sprintf(code_name,"%s.h",function_name);
      code = fopen(code_name,"w");
      fprintf(code,"#ifndef _libint_%s\n",function_name);
      fprintf(code,"#define _libint_%s\n",function_name);
      fprintf(code,"  /* This machine-generated function computes a quartet of (%c%c| integrals */\n\n",
	      am_letter[am_in[0]],am_letter[am_in[1]]);
      fprintf(code,"#define %s(AB, vp, I0, I1, cd_num)\\\n{\\\n",
	      function_name);
      fprintf(code,"  const REALTYPE AB0 = AB[0];\\\n");
      fprintf(code,"  const REALTYPE AB1 = AB[1];\\\n");
      fprintf(code,"  const REALTYPE AB2 = AB[2];\\\n");
      fprintf(code,"  REALTYPE *i0, *i1;\\\n");
      fprintf(code,"  int cd;\\\n");
      fprintf(code,"  REALTYPE *target = (vp);\\\n\\\n");

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
		fprintf(code,"  i0 = (I0) + %d*cd_num;\\\n",t0);
	      else
		fprintf(code,"  i0 = (I0);\\\n");
	      if (t1)
		fprintf(code,"  i1 = (I1) + %d*cd_num;\\\n",t1);
	      else
		fprintf(code,"  i1 = (I1);\\\n");

	      fprintf(code,"  for(cd=0;cd<cd_num;cd++)\\\n");
	      fprintf(code,"    *(target++) = *(i0++) + AB%d*(*(i1++));\\\n",xyz);

	      curr_count++;
	    }
	  }
	}
      }
      fprintf(code,"}\n");
      fprintf(code,"\n#endif\n"); /* end of #ifndef _libint_.... */
      fclose(code);
      printf("Done with %s\n",code_name);
    }
  }
}


