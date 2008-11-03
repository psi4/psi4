/*!
** \file
** \brief Initialize input, output, file prefix, etc.
** \ingroup CIOMR
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <psifiles.h>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>

#if defined HAVE_DECL_SETENV && !HAVE_DECL_SETENV
  extern int setenv(const char *, const char *, int);
#endif

namespace psi {

static char *ifname = NULL;
static char *ofname = NULL;
static char *fprefix = NULL;

/*!
** psi_start():
** This function initializes the input, output files, file prefix, etc.,
** by checking command line arguments and environmental variables. It also 
** initializes the Input Parsing library.
**
** \param argc       = number of command-line arguments passed
** \param argv       = command-line arguments
** \param overwrite  = whether to overwrite output file (1) or append to it (0).
**                     Most PSI modules will want to append.
**
** Returns: one of standard PSI error codes
** \ingroup CIOMR
*/
int psi_start(FILE** infile, FILE** outfile, char** psi_file_prefix, 
  int argc, char *argv[], int overwrite_output)
{
  int i, errcod;
                                /* state flags */
  int found_if_np = 0;          /* found input file name without -f */
  int found_of_np = 0;          /* found output file name without -o */
  int found_fp_np = 0;          /* found file prefix name without -p */
  int found_if_p = 0;           /* found input file name with -f */
  int found_of_p = 0;           /* found output file name with -o */
  int found_fp_p = 0;           /* found file prefix name with -p */
  char *cfname = NULL;
  char *userhome;
  FILE *psirc;
  char *arg;
  char *tmpstr1;

  /* process command-line arguments in sequence */
  for(i=0; i<argc; i++) {
    arg = argv[i];
    if (!strcmp(arg,"-f") && !found_if_p) {
      ifname = argv[++i];
      found_if_p = 1;
    }
    else if (!strcmp(arg,"-o") && !found_of_p) {
      ofname = argv[++i];
      found_of_p = 1;
    }
    else if (!strcmp(arg,"-p") && !found_fp_p) {
      fprefix = argv[++i];
      found_fp_p = 1;
    }
    else if (arg[0] == '-') {
      fprintf(stderr, "Error: unrecognized command-line argument %s\n", arg);
      return(PSI_RETURN_FAILURE);
    }
    else if (!found_if_np) {
      ifname = arg;
      found_if_np = 1;
    }
    else if (!found_of_np) {
      ofname = arg;
      found_of_np = 1;
    }
    else if (!found_fp_np) {
      fprefix = arg;
      found_fp_np = 1;
    }
    else {
      fprintf(stderr, "Error: too many command-line arguments given\n");
      return(PSI_RETURN_FAILURE);
    }
  }
  

  /* check if some args were specified in both prefixed and nonprefixed form */
  if (found_if_p && found_if_np) {
    fprintf(stderr, 
      "Error: input file name specified both with and without -f\n");
    fprintf(stderr, 
      "Usage: (module) [options] -f input -o output [-p prefix]  OR\n");
    fprintf(stderr, "       (module) [options] input output [prefix]\n");
    return(PSI_RETURN_FAILURE);
  }
  if (found_of_p && found_of_np) {
    fprintf(stderr, 
      "Error: output file name specified both with and without -o\n");
    fprintf(stderr, 
      "Usage: (module) [options] -f input -o output [-p prefix]  OR\n");
    fprintf(stderr, "       (module) [options] input output [prefix]\n");
    return(PSI_RETURN_FAILURE);
  }
  if (found_fp_p && found_fp_np) {
    fprintf(stderr, 
      "Error: file prefix specified both with and without -p\n");
    fprintf(stderr, 
      "Usage: (module) [options] -f input -o output -p prefix  OR\n");
    fprintf(stderr, "       (module) [options] input output prefix\n");
    return(PSI_RETURN_FAILURE);
  }

  
  /* if some args were not specified on command-line - check the environment */
  if (ifname == NULL)
    ifname = getenv("PSI_INPUT");
  if (ofname == NULL)
    ofname = getenv("PSI_OUTPUT");
  if (fprefix == NULL)
    fprefix = getenv("PSI_PREFIX");

  /* if some arguments still not defined - assign default values */
  if (ifname == NULL)
    ifname = strdup("input.dat");
  if (ofname == NULL)
    ofname = strdup("output.dat");
  /* default prefix is not assigned here yet because need to check 
     input file first */

  /* open input and output files */
  if(ifname[0]=='-' && ifname[1]=='\x0')
    *infile=stdin;
  else
    *infile = fopen(ifname, "r");
  if (*infile == NULL) {
    fprintf(stderr, "Error: could not open input file %s\n",ifname);
    return(PSI_RETURN_FAILURE);
  }
  if (overwrite_output)
  {
    if(ofname[0]=='-' && ofname[1]=='\x0')
      *outfile=stdout;
    else
      *outfile = fopen(ofname, "w");
  }
  else
  {
    if(ofname[0]=='-' && ofname[1]=='\x0')
      *outfile=stdout;
    else
      *outfile = fopen(ofname, "a");
  }
  if (*outfile == NULL) {
    fprintf(stderr, "Error: could not open output file %s\n",ofname);
    return(PSI_RETURN_FAILURE);
  }

  /* initialize libipv1 */
  ip_set_uppercase(1);
  ip_initialize(*infile, *outfile);
  ip_cwk_clear();

  /* open user's PSI configuration file (default, $HOME/.psirc) */
  cfname = getenv("PSI_RC");
  if (cfname == NULL) {
    userhome = getenv("HOME");
    cfname = (char *) malloc((10+strlen(userhome))*sizeof(char));
    sprintf(cfname, "%s%s", userhome, "/.psirc");
    psirc = fopen(cfname, "r");
    free(cfname);
  }
  else psirc = fopen(cfname, "r");
  if(psirc != NULL) {
    ip_append(psirc, stderr);
    fclose(psirc);
  }  

  /* lastly, everybody needs DEFAULT and PSI sections */
  ip_cwk_add(const_cast<char*>(":DEFAULT"));
  ip_cwk_add(const_cast<char*>(":PSI"));

  /* if prefix still NULL - check input file */
  if (fprefix == NULL)
    errcod = ip_string(const_cast<char*>(":DEFAULT:FILES:DEFAULT:NAME"),&fprefix,0);
  if (fprefix == NULL)
    errcod = ip_string(const_cast<char*>(":DEFAULT:NAME"),&fprefix,0);
  if (fprefix == NULL)
    errcod = ip_string(const_cast<char*>(":PSI:FILES:DEFAULT:NAME"),&fprefix,0);
  if (fprefix == NULL)
    errcod = ip_string(const_cast<char*>(":PSI:NAME"),&fprefix,0);

  /* copy over file prefix, etc. into their appropriate variables */
  if (fprefix == NULL) {
    fprefix = strdup(PSI_DEFAULT_FILE_PREFIX);
  }
  *psi_file_prefix = strdup(fprefix);

  /* other Psi modules called by this module should read from the same input 
     file set the value of PSI_INPUT for the duration of this run */
#if HAVE_PUTENV
  tmpstr1 = (char *) malloc(11+strlen(ifname));
  sprintf(tmpstr1, "PSI_INPUT=%s", ifname);
  putenv(tmpstr1);  /* note potential memory leak */
#elif HAVE_SETENV
  setenv("PSI_OUTPUT",ifname,1);
#else
#error "Have neither putenv nor setenv. Something must be very broken on this system."
#endif

  /* By default, other Psi modules called by this module should write to 
     the same output file set the value of PSI_OUTPUT for the duration of 
     this run */
#if HAVE_PUTENV
  tmpstr1 = (char *) malloc(12+strlen(ofname));
  sprintf(tmpstr1, "PSI_OUTPUT=%s", ofname);
  putenv(tmpstr1); /* note potential memory leak */
#elif HAVE_SETENV
  setenv("PSI_OUTPUT",ofname,1);
#else
#error "Have neither putenv nor setenv. Something must be very broken on this system."
#endif

  /* By default, other Psi modules called by this module should use the same
     prefix too set the value of PSI_PREFIX for the duration of this run */
#if HAVE_PUTENV
  tmpstr1 = (char *) malloc(12+strlen(fprefix));
  sprintf(tmpstr1, "PSI_PREFIX=%s", fprefix);
  putenv(tmpstr1); /* note potential memory leak */
#elif HAVE_SETENV
  setenv("PSI_PREFIX",fprefix,1);
#else
#error "Have neither putenv nor setenv. Something must be very broken on this system."
#endif

  return(PSI_RETURN_SUCCESS);
}


/*!
** psi_ifname()
**
** This function returns the input file name
**
** Arguments: none
**
** Returns: the pointer to the string containing the input
**          file name if it has been determined, NULL otherwise
** \ingroup CIOMR
*/

char* psi_ifname()
{
  return ifname;
}

/*!
** psi_ofname()
**
** This function returns the output file name
**
** Arguments: none
**
** Returns: the pointer to the string containing the output
**          file name if it has been determined, NULL otherwise
** \ingroup CIOMR
*/

char* psi_ofname()
{
  return ofname;
}

/*!
** psi_fprefix()
**
** This function returns the PSI file prefix
**
** Arguments: none
**
** Returns: the pointer to the string containing the PSI
**          file prefix if it has been determined, NULL otherwise
** \ingroup CIOMR
*/

char* psi_fprefix()
{
  return fprefix;
}

}

