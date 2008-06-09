/*! \file
    \ingroup INT
    \brief Enter brief description of file here 
*/

/*! \defgroup INT libint: The Integral Library */

/*------------------------------------------------------------------------------------------------------
                                           BUILD_LIBINT
  		       Written by Dr. Justin T. Fermann and Edward F. Valeev
  This program generates files necessary for compiling the LIBINT library. LIBINT is a library of
  streamlined highly efficient routines for recursive computation of ERI integrals of the form (a0|b0).
  The library is used by the integral code CINTS. The following files must be present in the directory
  where the current program is intended to run: libint_config.h and Makefile.libint (all included with the source).
  libint_config.h contains set of defines necessary for building the library.
  libint_config.h must define the following macro variables:
  LIBINT_NEW_AM - twice the desired maximum angular momentum (default is 8);
  LIBINT_OPT_AM - twice the angular momentum up to which VRR Level 0 routines are machine generated.
  A generic VRR Level 0 function is used past this value. OPT_AM=8 should be enough for almost anyone.
  LIBINT_MAX_CLASS_SIZE - maximum length of _build_a0b0 (if (a0|b0) class is longer than MAX_CLASS_SIZE, the routine is
  going to be split into several smaller ones. This is done to prevent compiler from exhausting system resources.
  Defaults to 785).

  *EXAMPLE*: if one wants LIBINT work for up to (gg|gg) integrals (the current PSI limit), NEW_AM has to be set
  to 8 (g corresponds to l=4, (gg|gg) class will require at most the (l0|l0) class). The intended angular
  momentum limit for PSI 3 is i (l=6), therefore up to (q0|q0) classes are required. NEW_AM must be set to 12.

  Accessing functions in LIBINT is very simple - the program has to call init_libint() just once before it
  starts computing integrals. After that all top_build_... functions (****BUT**** top_build_0000, which should
  never be neccessary, since (ss|ss) class is easily computed from the auxiliary function) will be arranged
  in a matrix of pointers. E.g., to call top_build_i0f0(...) one has to invoke top_build_a0b0[6][3](...).
 ------------------------------------------------------------------------------------------------------*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "build_libint.h"
#include <libint/constants.h>

#include "libint_config.h"

/* Global data */
FILE *outfile, *vrr_header, *hrr_header, *libint_header, *init_code;
char *real_type;   /*--- C type for real numbers ---*/
int libint_stack_size[MAX_AM/2+1];
LibintParams_t Params;

void punt();
int emit_vrr_build();
int emit_vrr_build_macro();
int emit_order();
int emit_hrr_build();
int emit_hrr_build_macro();

int main()
{
  int i,j,k,l,f;
  int j_min, j_max, k_min, k_max, l_min, l_max;
  int errcod;
  int old_am = 0;
  int new_am, opt_am;
  int max_class_size = DEFAULT_MAX_CLASS_SIZE;
  int long_double = 0;     /*--- Whether to use long doubles ---*/
  int stack_size;

  /*-------------------------------
    Initialize files and libraries
   -------------------------------*/
  outfile = fopen("./output.dat", "w");
  vrr_header = fopen("./vrr_header.h","w");
  hrr_header = fopen("./hrr_header.h","w");
  libint_header = fopen("./libint.h","w");
  init_code = fopen("./init_libint.cc","w");

  new_am = LIBINT_NEW_AM;
  if (new_am <= 0)
    punt("  MAX_AM must be positive.");
  if (new_am > MAX_AM)
    punt("  Maximum MAX_AM exceeded. Contact the program author.");

  opt_am = LIBINT_OPT_AM;
  if (opt_am < 2)
    opt_am = DEFAULT_OPT_AM;
  if (opt_am > new_am) opt_am = new_am;

  max_class_size = LIBINT_MAX_CLASS_SIZE;
  if (max_class_size < 10)
    punt("  MAX_CLASS_SIZE cannot be smaller than 10.");

  long_double = LIBINT_LONG_DOUBLE;
  if (long_double)
    real_type = strdup("long double");
  else
    real_type = strdup("double");

  /*-------------
    Init globals
   -------------*/
  for(l=0;l<=new_am/2;l++)
    libint_stack_size[l] = 1;  /* 1 is a safe value and ensures that int_stack is allocated */
  Params.new_am = new_am;
  Params.old_am = 0;
  Params.opt_am = opt_am;
  Params.max_class_size = max_class_size;
  Params.max_am_to_inline_vrr_worker = -1;
  Params.max_am_manager_to_inline_vrr_worker = -1;
  Params.max_am_to_inline_hrr_worker = -1;
  Params.max_am_manager_to_inline_hrr_worker = -1;
  Params.max_am_to_inline_vrr_manager = -1;

  /* Setting up init_libint.c, header.h */
  fprintf(init_code,"#include <stdio.h>\n");
  fprintf(init_code,"#include <stdlib.h>\n");
  fprintf(init_code,"#include <psi4.h>\n");
  fprintf(init_code,"#include \"libint.h\"\n");
  fprintf(init_code,"#include \"hrr_header.h\"\n\n");
  fprintf(init_code,"extern \"C\" {\n");
  fprintf(init_code,"REALTYPE *(*build_eri[%d][%d][%d][%d])(Libint_t *, int);\n",new_am/2+1,new_am/2+1,new_am/2+1,new_am/2+1);
  fprintf(init_code,"int libint_stack_size[%d];\n\n",new_am/2+1);
  fprintf(init_code,"/* This function initializes a matrix of pointers to routines */\n");
  fprintf(init_code,"/* for computing ERI classes up to (%cs|%cs) - the base of the library */\n\n",
	  am_letter[new_am],am_letter[new_am]);
  fprintf(init_code,"void init_libint_base()\n{\n");

  /* Declare generic build routines */
  fprintf(vrr_header,"extern \"C\" REALTYPE *vrr_build_xxxx(int am[2], prim_data *, REALTYPE *, const REALTYPE *,");
  fprintf(vrr_header,"const REALTYPE *, const REALTYPE *, const REALTYPE *, const REALTYPE *);\n\n");

  emit_order();
  emit_vrr_build();
  emit_vrr_build_macro();
  emit_hrr_build();
  emit_hrr_build_macro();

  /* put computed stack sizes for each angular momentum level into init_libint_base() */
  for(l=0;l<=new_am/2;l++)
    fprintf(init_code,"\n  libint_stack_size[%d] = %d;",l,libint_stack_size[l]);
  
  fprintf(init_code,"\n}\n");
  fprintf(init_code,"/* These functions initialize library objects */\n");
  fprintf(init_code,"/* Library objects operate independently of each other */\n");
  fprintf(init_code,"int init_libint(Libint_t *libint, int max_am, int max_num_prim_quartets)\n{\n");
  fprintf(init_code,"  int memory = 0;\n\n");
  fprintf(init_code,"  if (max_am >= LIBINT_MAX_AM) return -1;\n");
  fprintf(init_code,"  libint->int_stack = (REALTYPE *) malloc(libint_stack_size[max_am]*sizeof(REALTYPE));\n");
  fprintf(init_code,"  memory += libint_stack_size[max_am];\n");
  fprintf(init_code,"  libint->PrimQuartet = (prim_data *) malloc(max_num_prim_quartets*sizeof(prim_data));\n");
  fprintf(init_code,"  memory += max_num_prim_quartets*sizeof(prim_data)/sizeof(REALTYPE);\n");
  fprintf(init_code,"  return memory;\n}\n\n");
  fprintf(init_code,"void free_libint(Libint_t *libint)\n{\n");
  fprintf(init_code,"  if (libint->int_stack != NULL) {\n");
  fprintf(init_code,"    free(libint->int_stack);\n");
  fprintf(init_code,"    libint->int_stack = NULL;\n");
  fprintf(init_code,"  }\n");
  fprintf(init_code,"  if (libint->PrimQuartet != NULL) {\n");
  fprintf(init_code,"    free(libint->PrimQuartet);\n");
  fprintf(init_code,"    libint->PrimQuartet = NULL;\n");
  fprintf(init_code,"  }\n\n");
  fprintf(init_code,"  return;\n}\n\n");
  fprintf(init_code,"int libint_storage_required(int max_am, int max_num_prim_quartets)\n{\n");
  fprintf(init_code,"  int memory = 0;\n\n");
  fprintf(init_code,"  if (max_am >= LIBINT_MAX_AM) return -1;\n");
  fprintf(init_code,"  memory += libint_stack_size[max_am];\n");
  fprintf(init_code,"  memory += max_num_prim_quartets*sizeof(prim_data)/sizeof(REALTYPE);\n");
  fprintf(init_code,"  return memory;\n}\n\n");
  fprintf(init_code,"}\n"); /* end of extern "C" */
  fclose(init_code);
  fclose(vrr_header);
  fclose(hrr_header);
  
    /* Setting up libint.h */
  fprintf(libint_header,"#ifndef _psi3_libint_h\n");
  fprintf(libint_header,"#define _psi3_libint_h\n\n");
  fprintf(libint_header,"/* Maximum angular momentum of functions in a basis set plus 1 */\n");
  fprintf(libint_header,"#define REALTYPE %s\n",real_type);
  if (long_double)
    fprintf(libint_header,"#define NONDOUBLE_INTS\n");
  fprintf(libint_header,"#define LIBINT_MAX_AM %d\n",1+new_am/2);
  fprintf(libint_header,"#define LIBINT_OPT_AM %d\n",1+opt_am/2);
  fprintf(libint_header,"typedef struct pdata{\n");
  fprintf(libint_header,"  REALTYPE F[%d];\n",2*new_am+1);
  fprintf(libint_header,"  REALTYPE U[6][3];\n");
  fprintf(libint_header,"  REALTYPE twozeta_a;\n"); 
  fprintf(libint_header,"  REALTYPE twozeta_b;\n"); 
  fprintf(libint_header,"  REALTYPE twozeta_c;\n");
  fprintf(libint_header,"  REALTYPE twozeta_d;\n");
  fprintf(libint_header,"  REALTYPE oo2z;\n");
  fprintf(libint_header,"  REALTYPE oo2n;\n");
  fprintf(libint_header,"  REALTYPE oo2zn;\n");
  fprintf(libint_header,"  REALTYPE poz;\n");
  fprintf(libint_header,"  REALTYPE pon;\n");
  fprintf(libint_header,"  REALTYPE oo2p;\n");
  fprintf(libint_header,"  REALTYPE ss_r12_ss;\n");
  fprintf(libint_header,"  } prim_data;\n\n");
  fprintf(libint_header,"typedef struct {\n");
  fprintf(libint_header,"  REALTYPE *int_stack;\n"); 
  fprintf(libint_header,"  prim_data *PrimQuartet;\n"); 
  fprintf(libint_header,"  REALTYPE AB[3];\n");
  fprintf(libint_header,"  REALTYPE CD[3];\n");
  fprintf(libint_header,"  REALTYPE *vrr_classes[%d][%d];\n",1+new_am,1+new_am);
  fprintf(libint_header,"  REALTYPE *vrr_stack;\n");
  fprintf(libint_header,"  } Libint_t;\n\n");
  fprintf(libint_header,"#ifdef __cplusplus\n");
  fprintf(libint_header,"extern \"C\" {\n");
  fprintf(libint_header,"#endif\n");
  fprintf(libint_header,"extern REALTYPE *(*build_eri[%d][%d][%d][%d])(Libint_t *, int);\n",
	  new_am/2+1,new_am/2+1,new_am/2+1,new_am/2+1);
  fprintf(libint_header,"void init_libint_base();\n");
  fprintf(libint_header,"int  init_libint(Libint_t *, int max_am, int max_num_prim_comb);\n");
  fprintf(libint_header,"void free_libint(Libint_t *);\n");
  fprintf(libint_header,"int  libint_storage_required(int max_am, int max_num_prim_comb);\n");
  fprintf(libint_header,"#ifdef __cplusplus\n");
  fprintf(libint_header,"}\n");
  fprintf(libint_header,"#endif\n\n");
  fprintf(libint_header,"#endif\n");
  fclose(libint_header);
  fclose(outfile);
  exit(0);
}


void punt(char* str)
{
  printf(str);
  exit(1);
}


