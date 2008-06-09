/*! \file
    \ingroup DERIV
    \brief Enter brief description of file here 
*/
#include <math.h>
#include <stdio.h>
#include "build_libderiv.h"

extern FILE *outfile, *libint_src, *d1hrr_header;
extern LibderivParams_t Params;

extern void punt(char *);

int emit_d1hrr_build()
{
  int new_am = Params.new_am;
  int old_am = Params.old_am;
  int am_to_inline = Params.max_am_to_inline_d1hrr_worker;

  FILE *code;
  int p,q,r,s;
  int ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz;
  int i,j,nj,i_i0,i_i1;
  int k,l,nl,k_i0,k_i1;
  int i0_step,i1_step;
  int a, b;
  int flag;
  int am_in[2];
  int class_size;
  int la, lb;
  int ld, lc, ld_max;
  int curr_count,curr_subfunction;
  int current_highest_am, to_inline;
  int errcod;
  static int io[] = {1,3,6,10,15,21,28,36,45,55,66,78,91,105,120,136,153};
  static const char am_letter[] = "0pdfghiklmnoqrtuvwxyz";
  char code_name[19];
  char function_name[17];
  

  for(lc=0;lc<=new_am;lc++) {
    ld_max = (lc+1)/2;
    for(ld=1;ld<=ld_max;ld++) {

      /*-----------------------
	HRR on centers C and D
       -----------------------*/

      am_in[0] = lc-ld;
      am_in[1] = ld;

      /* Is this function to be made inline */
      current_highest_am = (am_in[0] > am_in[1]) ? am_in[0] : am_in[1];
      to_inline = (current_highest_am <= am_to_inline) ? 1 : 0;

      sprintf(function_name,"d1hrr3_build_%c%c",am_letter[am_in[0]],am_letter[am_in[1]]);
      sprintf(code_name,"d1hrr3_build_%c%c.cc",am_letter[am_in[0]],am_letter[am_in[1]]);
      code = fopen(code_name,"w");

      /* include the function into the d1hrr_header.h */
      if (to_inline)
	fprintf(d1hrr_header, "#ifndef INLINE_D1HRR_WORKER\n");
      fprintf(d1hrr_header,"void %s(const double *, double *, const double *, const double *, \n",function_name);
      fprintf(d1hrr_header,"        double, const double *, double, const double *, double, const double *,\n");
      fprintf(d1hrr_header,"        double, const double *, double, const double *, double, const double *, int);\n");
      if (to_inline)
	fprintf(d1hrr_header, "#endif\n");

      fprintf(code,"  /* This machine-generated function computes a quartet of |%c%c) first derivative ERIs */\n\n",
	      am_letter[am_in[0]],am_letter[am_in[1]]);
      fprintf(code,"void %s(const double *CD, double *vp, const double *I0, const double *I1,\n",function_name);
      fprintf(code,"        double c2, const double *I2, double c3, const double *I3, double c4, const double *I4,\n");
      fprintf(code,"        double c5, const double *I5, double c6, const double *I6, double c7, const double *I7, int ab_num)\n{\n");
      fprintf(code,"  int ab;\n");
      fprintf(code,"  const double CD0 = CD[0];\n");
      fprintf(code,"  const double CD1 = CD[1];\n");
      fprintf(code,"  const double CD2 = CD[2];\n");

      nl = (am_in[1]*(am_in[1]+1))/2;
      i0_step = (am_in[0]+2)*(am_in[0]+3)*nl/2;
      i1_step = (am_in[0]+1)*(am_in[0]+2)*nl/2;
      fprintf(code,"  for(ab=0;ab<ab_num;ab++) {\n");

      for(p = 0; p <= am_in[0]; p++){
	cx = am_in[0] - p;
	for(q = 0; q <= p; q++){
	  cy = p - q;
	  cz = q;
	  k_i1 = io[p]-p+q-1;
	  
	  for(r = 0; r <= am_in[1]; r++){
	    dx = am_in[1] - r;
	    for(s = 0; s <= r; s++){
	      dy = r - s;
	      dz = s;

	      if (dx) { /* build along x */
		k_i0 = k_i1;
		l = io[r]-r+s-1;
		fprintf(code, "    *(vp++) = I0[%d] + CD0*I1[%d] + c2*I2[%d] - c5*I5[%d];\n",
			k_i0*nl+l,k_i1*nl+l,k_i1*nl+l,k_i1*nl+l);
	      }
	      else if (dy) { /* build along y */
		k_i0 = io[p+1]-p+q-2;
		l = io[r-1]-r+s;
		fprintf(code, "    *(vp++) = I0[%d] + CD1*I1[%d] + c3*I3[%d] - c6*I6[%d];\n",
			k_i0*nl+l,k_i1*nl+l,k_i1*nl+l,k_i1*nl+l);
	      }
	      else { /* build along z */
		k_i0 = io[p+1]-p+q-1;
		l = io[r-1]-r+s-1;
		fprintf(code, "    *(vp++) = I0[%d] + CD2*I1[%d] + c4*I4[%d] - c7*I7[%d];\n",
			k_i0*nl+l,k_i1*nl+l,k_i1*nl+l,k_i1*nl+l);
	      }
	    }
	  }
	}
      }
      fprintf(code,"    I0 += %d;\n    I1 += %d;\n",i0_step,i1_step);
      fprintf(code,"    I2 += %d;\n    I3 += %d;\n    I4 += %d;\n",i1_step,i1_step,i1_step);
      fprintf(code,"    I5 += %d;\n    I6 += %d;\n    I7 += %d;\n",i1_step,i1_step,i1_step);
      fprintf(code,"  }\n}\n");
      fclose(code);
      printf("Done with %s\n",code_name);
      
      
      /*-----------------------
	HRR on centers A and B
       -----------------------*/

      la = lc-ld;  lb = ld;
      am_in[0] = la;
      am_in[1] = lb;

      sprintf(function_name,"d1hrr1_build_%c%c",am_letter[am_in[0]],am_letter[am_in[1]]);
      sprintf(code_name,"d1hrr1_build_%c%c.cc",am_letter[am_in[0]],am_letter[am_in[1]]);
      code = fopen(code_name,"w");

      /* include the function into the d1hrr_header.h */
      if (to_inline)
	fprintf(d1hrr_header, "#ifndef INLINE_D1HRR_WORKER\n");
      fprintf(d1hrr_header,"void %s(const double *, double *, const double *, const double *,",function_name);
      fprintf(d1hrr_header,"        double, const double *, double, const double *, double, const double *,\n");
      fprintf(d1hrr_header,"        double, const double *, double, const double *, double, const double *, int);\n");
      if (to_inline)
	fprintf(d1hrr_header, "#endif\n");
      
      fprintf(code,"  /* This machine-generated function computes a quartet of (%c%c| first derivative integrals */\n\n",
	      am_letter[am_in[0]],am_letter[am_in[1]]);
      fprintf(code,"void %s(const double *AB, double *vp, const double *I0, const double *I1,\n",function_name);
      fprintf(code,"        double c2, const double *I2, double c3, const double *I3, double c4, const double *I4,\n");
      fprintf(code,"        double c5, const double *I5, double c6, const double *I6, double c7, const double *I7, int cd_num)\n{\n");
      fprintf(code,"  int cd;\n");
      fprintf(code,"  const double *i0, *i1, *i2, *i3, *i4, *i5, *i6, *i7;\n");
      fprintf(code,"  const double AB0 = AB[0];\n");
      fprintf(code,"  const double AB1 = AB[1];\n");
      fprintf(code,"  const double AB2 = AB[2];\n");

      nj = (lb*(lb+1))/2;

      for(p = 0; p <= am_in[0]; p++){
	ax = am_in[0] - p;
	for(q = 0; q <= p; q++){
	  ay = p - q;
	  az = q;
	  i_i1 = io[p]-p+q-1;
	  
	  for(r = 0; r <= am_in[1]; r++){
	    bx = am_in[1] - r;
	    for(s = 0; s <= r; s++){
	      by = r - s;
	      bz = s;

	      if (bx) { /* build along x */
		i_i0 = i_i1;
		j = io[r]-r+s-1;
		if (i_i0*nj+j)
		  fprintf(code,"  i0 = I0 + %d*cd_num;\n",i_i0*nj+j);
		else
		  fprintf(code,"  i0 = I0;\n");
		if (i_i1*nj+j)
		  fprintf(code,"  i1 = I1 + %d*cd_num;\n",i_i1*nj+j);
		else
		  fprintf(code,"  i1 = I1;\n");
		if (i_i1*nj+j)
		  fprintf(code,"  i2 = I2 + %d*cd_num;\n",i_i1*nj+j);
		else
		  fprintf(code,"  i2 = I2;\n");
		if (i_i1*nj+j)
		  fprintf(code,"  i5 = I5 + %d*cd_num;\n",i_i1*nj+j);
		else
		  fprintf(code,"  i5 = I5;\n");
		fprintf(code,"  for(cd=0;cd<cd_num;cd++)\n");
		fprintf(code,"    *(vp++) = *(i0++) + AB0*(*(i1++)) + c2*(*(i2++)) - c5*(*(i5++));\n");
	      }
	      else if (by) { /* build along y */
		i_i0 = io[p+1]-p+q-2;
		j = io[r-1]-r+s;
		if (i_i0*nj+j)
		  fprintf(code,"  i0 = I0 + %d*cd_num;\n",i_i0*nj+j);
		else
		  fprintf(code,"  i0 = I0;\n");
		if (i_i1*nj+j)
		  fprintf(code,"  i1 = I1 + %d*cd_num;\n",i_i1*nj+j);
		else
		  fprintf(code,"  i1 = I1;\n");
		if (i_i1*nj+j)
		  fprintf(code,"  i3 = I3 + %d*cd_num;\n",i_i1*nj+j);
		else
		  fprintf(code,"  i3 = I3;\n");
		if (i_i1*nj+j)
		  fprintf(code,"  i6 = I6 + %d*cd_num;\n",i_i1*nj+j);
		else
		  fprintf(code,"  i6 = I6;\n");
		fprintf(code,"  for(cd=0;cd<cd_num;cd++)\n");
		fprintf(code,"    *(vp++) = *(i0++) + AB1*(*(i1++)) + c3*(*(i3++)) - c6*(*(i6++));\n");
	      }
	      else { /* build along z */
		i_i0 = io[p+1]-p+q-1;
		j = io[r-1]-r+s-1;
		if (i_i0*nj+j)
		  fprintf(code,"  i0 = I0 + %d*cd_num;\n",i_i0*nj+j);
		else
		  fprintf(code,"  i0 = I0;\n");
		if (i_i1*nj+j)
		  fprintf(code,"  i1 = I1 + %d*cd_num;\n",i_i1*nj+j);
		else
		  fprintf(code,"  i1 = I1;\n");
		if (i_i1*nj+j)
		  fprintf(code,"  i4 = I4 + %d*cd_num;\n",i_i1*nj+j);
		else
		  fprintf(code,"  i4 = I4;\n");
		if (i_i1*nj+j)
		  fprintf(code,"  i7 = I7 + %d*cd_num;\n",i_i1*nj+j);
		else
		  fprintf(code,"  i7 = I7;\n");
		fprintf(code,"  for(cd=0;cd<cd_num;cd++)\n");
		fprintf(code,"    *(vp++) = *(i0++) + AB2*(*(i1++)) + c4*(*(i4++)) - c7*(*(i7++));\n");
	      }
	    }
	  }
	}
      }
      fprintf(code,"}\n");
      fclose(code);
      printf("Done with %s\n",code_name);
    }
  }
}


