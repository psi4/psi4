/*! 
** \file
** \ingroup GEOM
** \brief Enter brief description of file here 
*/

/*
** READ_GEOM() : Function reads geometry (bohr) from file (usu. geom.dat)
**
** David Sherrill  -  June 1993
**
** Last modified 9/1/94 to remove dependency on ds_io library and to be
** a little more flexible parsing the parentheses.
**
*/

/* include's */
#include <cstdio>
#include <cstdlib>
#include <cstring>

/* define's */
#define MAX_LINE 132
#define MAX_ATOM_LABEL 10

namespace psi { namespace geom {

extern void malloc_ck(void *array, const char *mesg);
extern void prf_abort(FILE *file, const char *mesg);

/* global variables */
static char line1[MAX_LINE+1];

extern int label2an(char *label);


/*
** READ_GEOM(): This function reads geometry from file and stores it in an
**      array.
**
** Arguments: 
**      maxlines  = maximum number of lines to read
**      natom     = pointer to variable to hold number of atoms
**      X, Y, Z   = pointers to arrays of cartesian coordinates (assume bohr)
**                  (these are currently allocated in THIS function)
**      fname     = filename for input file 
**
** Returns: success (1) or failure (0) 
**
** NOTE : READS FIRST (MOST RECENT) ENTRY OF INPUT FILE
**        Also, according to the format of geom.dat, we must strip all leading
**        comment lines, and then read until the next comment line (or EOF)
**        is encountered.  Actually, also beware that the geometry data can
**        be enclosed in "geometry = ( ... )"
*/
int read_geom(int maxlines, int* natom, double** X, double** Y, double** Z, 
  char* fname) 
{
int i ;                      /* loop variable */
int datalines = 0 ;          /* number of data-containing lines read */
int commentlines = 0 ;       /* number of comment lines */
char *parenpos ;             /* point at position of parens in data lines */
FILE *fpi ;                  /* file pointer for input */

/* open input file */
   if ( (fpi = fopen(fname, "r")) == NULL) {
      fprintf(stderr, "(read_geom): Trouble opening file\n") ;
      return(0) ;
      }

/* make sure input file is at beginning */
   rewind(fpi) ;

/* read through all initial comment lines */
   while (1) {
      if (fgets(line1, MAX_LINE, fpi) == NULL) { /* error reading line */
         return(0) ;
         }
      else if ( (strchr(line1, '%') == NULL) && 
         (strstr(line1, "geometry =") == NULL) ) break ; /* if non-comment */
      else commentlines++ ;
      }

/* now read until we reach EOF or the next comment line */
   datalines = 1 ;
   while (datalines <= maxlines) {
      if (fgets(line1, MAX_LINE, fpi) == NULL) { /* error reading line */
         if (feof(fpi)) break ;
         else return(0) ;
         }
      else if (strchr(line1, '%') != NULL) break;       /* brk on comnt ln */

      i=0;
      while (line1[i] == ' ') i++;

      if (line1[i] == ')') break;                       /* brk if only paren */
      else datalines++ ;
      }

   if (datalines > maxlines) {
      fprintf(stderr, "(read_geom): Read too many lines\n") ;
      return(0) ;
      }

   *natom = datalines ;

/* make room for the Cartesian coordinates */
   *X = (double *) malloc (*natom * sizeof(double)) ;
   *Y = (double *) malloc (*natom * sizeof(double)) ;
   *Z = (double *) malloc (*natom * sizeof(double)) ;
   if (*X == NULL || *Y == NULL || *Z == NULL) {
      fprintf(stderr, "(read_geom): Trouble allocating Cartesian array\n");
      return(0);
      }

/* rewind file, read through comment lines, and then input data lines */
   rewind(fpi) ;
   for (i=0; i<commentlines; i++) {
      if (fgets(line1, MAX_LINE, fpi) == NULL) { /* error reading line */
         fprintf(stderr,"(read_geom): Read error on second comment pass\n") ;
         return(0) ;
         }
      }

/* read Cartesians */
      for (i=0; i<(*natom); i++) {
         if (fgets(line1, MAX_LINE, fpi) == NULL) { /* error reading line */
            fprintf(stderr,"(read_geom): Read error on second data pass\n");
            return(0) ;
            }
         /* strip out open and close parentheses (if any) */
         parenpos = strchr(line1, '(') ;
         if (parenpos != NULL) (*parenpos) = ' ' ;
         parenpos = strchr(line1, ')') ;
         if (parenpos != NULL) (*parenpos) = ' ' ;
       
         if (sscanf(line1, "%lf %lf %lf", (*X+i), (*Y+i), (*Z+i)) != 3) {
            fprintf(stderr,"(read_geom): Trouble reading Cartesian coords\n");
            return(0) ;
            }
      }

/* close up and return successfully */
   fclose(fpi) ;
   return(1) ;
}


/*
** PRINT_GEOM(): Function prints out the the geometry read in from
**      geom.dat
**
**
** Arguments: 
**      natom   = number of atoms
**      X, Y, Z = arrays of cartesian coordinates (assume bohr)
**      fpo     = file pointer for output
*/
void print_geom(int natom, double *X, double *Y, double *Z, FILE *fpo) 
{
   int i;

   for (i=0; i<natom; i++) {
      fprintf(fpo, "     %12.7lf    %12.7lf    %12.7lf\n",
            X[i], Y[i], Z[i]) ;
      }
   fprintf(fpo, "\n") ;
   fflush(fpo) ;
}





/*
** READ_ACES_GEOM()
**
** This function reads geometry data in the form of output from ACES II
**
** Arguments: 
**      fname   = filename to read
**      natom   = pointer to variable to hold number of atoms
**      X, Y, Z = pointers to arrays of cartesian coordinates (assume bohr)
**                (these are currently allocated in THIS function)
**      AN      = pointer to atomic number array (allocated here)
**      fpo     = file pointer for output
**
** Returns: success (1) or failure (0) 
**
*/
int read_aces_geom(char *fname, int *natom, double **X, double **Y, 
		   double **Z, double **AN, FILE *fpo) 
{
  int i;                          /* loop variable */
  int datalines = 0;              /* number of lines of input file */
  FILE *fpi;                      /* pointer for input file */
  char alabel[MAX_ATOM_LABEL];
  
  /* open file for reading */
  if ((fpi = fopen(fname, "r")) == NULL) {
    printf("(read_aces_geom): Error opening %s\n", fname);
    exit(0);
  }

  /* get the number of data lines */
  while (fgets(line1, MAX_LINE, fpi) != NULL) 
    datalines++;

  *natom = datalines;

  /* now make room for the Cartesian coordinates */
  *X = (double *) malloc (*natom * sizeof(double));
  malloc_ck((void*)*X, "(read_file11): Trouble allocating Cartesian array\n");
  *Y = (double *) malloc (*natom * sizeof(double));
  malloc_ck((void*)*Y, "(read_file11): Trouble allocating Cartesian array\n");
  *Z = (double *) malloc (*natom * sizeof(double));
  malloc_ck((void*)*Z, "(read_file11): Trouble allocating Cartesian array\n");
  *AN = (double *) malloc (*natom * sizeof(double));
  malloc_ck((void*)*AN, "(read_file11): Trouble allocating atomic num array\n");
  
  /* read symbol, atomic number, and Cartesians */
  rewind(fpi);
  i=0;
  while (fgets(line1, MAX_LINE, fpi) != NULL) {
    if (sscanf(line1, "%s %lf %lf %lf %lf", alabel, (*AN+i), (*X+i),
	       (*Y+i), (*Z+i)) != 5) 
      prf_abort(fpo, "(read_aces_geom): Trouble reading data\n");
    i++; 
  }

  *natom = i;
  fclose(fpi);
  return(1);
}



/*
** PRINT_ACES_GEOM()
**
** Function prints out the geometry information obtained in ACES format
** (or XYZ)
**
** Arguments:
**      natom   = number of atoms
**      X, Y, Z = arrays of cartesian coordinates (assume bohr)
**      AN      = atomic number array (allocated here)
**      labels  = atomic labels
**      fpo     = file pointer for output
*/
void print_aces_geom(int natom, double *X, double *Y, double *Z, double *AN, 
		     const char **labels, FILE *fpo) 
{
  int i;
  
  /* fprintf(fpo, "DATA FROM INPUT\n"); */
  /* fprintf(fpo, "Number of atoms = %d\n", natom); */
  /* fprintf(fpo, "\nCartesian coordinates (bohr) :\n"); */
  for (i=0; i<natom; i++) {
    fprintf(fpo, "     %3s    %12.7lf    %12.7lf    %12.7lf\n",
            labels[(int) AN[i]], X[i], Y[i], Z[i]);
    /*
    fprintf(fpo, "     %11.6lf    %12.7lf    %12.7lf    %12.7lf\n",
            AN[i], X[i], Y[i], Z[i]);
    */
  }
  fprintf(fpo, "\n");
  fflush(fpo);
}



/*
** READ_QCHEM_GEOM()
**
** This function reads geometry data in the form of output from QCHEM
**
** Arguments: 
**      fname   = filename to read
**      natom   = pointer to variable to hold number of atoms
**      X, Y, Z = pointers to arrays of cartesian coordinates (assume bohr)
**                (these are currently allocated in THIS function)
**      AN      = pointer to array of atomic numbers (alloc'd here)
**      fpo     = file pointer for output
**
** Returns: success (1) or failure (0) 
**
*/
int read_qchem_geom(char *fname, int *natom, double **X, double **Y, 
                    double **Z, double **AN, FILE *fpo) 
{
  int i;                          /* loop variable */
  int datalines = 0;              /* number of lines of input file */
  FILE *fpi;                      /* pointer for input file */
  char alabel[MAX_ATOM_LABEL];
  int junk;                       /* for holding the atom numbers */
  
  /* open file for reading */
  if ((fpi = fopen(fname, "r")) == NULL) {
    printf("(read_qchem_geom): Error opening %s\n", fname);
    exit(0);
  }

  /* get the number of data lines */
  while (fgets(line1, MAX_LINE, fpi) != NULL) 
    datalines++;

  *natom = datalines;

  /* now make room for the Cartesian coordinates */
  *X = (double *) malloc (*natom * sizeof(double));
  malloc_ck((void*)*X, "(read_file11): Trouble allocating Cartesian array\n");
  *Y = (double *) malloc (*natom * sizeof(double));
  malloc_ck((void*)*Y, "(read_file11): Trouble allocating Cartesian array\n");
  *Z = (double *) malloc (*natom * sizeof(double));
  malloc_ck((void*)*Z, "(read_file11): Trouble allocating Cartesian array\n");
  *AN = (double *) malloc (*natom * sizeof(double));
  malloc_ck((void*)*AN, "(read_file11): Trouble allocating atomic num array\n");
  
  /* read atom number and Cartesians */
  rewind(fpi);
  i=0;
  while (fgets(line1, MAX_LINE, fpi) != NULL) {
    if (sscanf(line1, "%d %s %lf %lf %lf", &junk, alabel, (*X+i),
	       (*Y+i), (*Z+i)) != 5) 
      prf_abort(fpo, "(read_qchem_geom): Trouble reading data\n");
    junk = label2an(alabel);
    if (junk == -1) 
      prf_abort(fpo, "(read_qchem_geom): Trouble parsing atom label\n");
    (*AN)[i] = junk;
    i++; 
  }

  *natom = i;
  fclose(fpi);
  return(1);
}


/*
** READ_XYZ_GEOM()
**
** This function reads geometry data in the form of a standard XYZ
** file, which has the number of atoms on the first line, a commnent
** on the second line, and following lines contain an atomic symbol
** followed by X Y and Z coordinates.
**
** Arguments: 
**      fname   = filename to read
**      natom   = pointer to variable to hold number of atoms
**      X, Y, Z = pointers to arrays of cartesian coordinates (assume bohr)
**                (these are currently allocated in THIS function)
**      AN      = pointer to array of atomic numbers (alloc'd here)
**      fpo     = file pointer for output
**
** Returns: success (1) or failure (0) 
**
*/
int read_xyz_geom(char *fname, int *natom, double **X, double **Y, 
                  double **Z, double **AN, FILE *fpo) 
{
  int i;                          /* loop variable */
  FILE *fpi;                      /* pointer for input file */
  char alabel[MAX_ATOM_LABEL];
  int junk;                       /* for holding the atom numbers */
  
  /* open file for reading */
  if ((fpi = fopen(fname, "r")) == NULL) {
    printf("(read_xyz_geom): Error opening %s\n", fname);
    exit(0);
  }

  /* get the number of data lines */
  fgets(line1, MAX_LINE, fpi);
  if (sscanf(line1, "%d", natom) != 1) {
    printf("(read_xyz_geom): Trouble reading number of atoms\n");
    exit(0);
  }

  /* now make room for the Cartesian coordinates */
  *X = (double *) malloc (*natom * sizeof(double));
  malloc_ck((void*)*X, "(read_xyz_geom): Trouble allocating Cartesian array\n");
  *Y = (double *) malloc (*natom * sizeof(double));
  malloc_ck((void*)*Y, "(read_xyz_geom): Trouble allocating Cartesian array\n");
  *Z = (double *) malloc (*natom * sizeof(double));
  malloc_ck((void*)*Z, "(read_xyz_geom): Trouble allocating Cartesian array\n");
  *AN = (double *) malloc (*natom * sizeof(double));
  malloc_ck((void*)*AN, "(read_xyz_geom): Trouble allocating atomic num array\n");
  
  /* read the comment line and ignore it */
  fgets(line1, MAX_LINE, fpi);
  
  /* read atom number and Cartesians */
  i=0;
  while (fgets(line1, MAX_LINE, fpi) != NULL) {
    if (sscanf(line1, "%s %lf %lf %lf", alabel, (*X+i),
	       (*Y+i), (*Z+i)) != 4) 
      prf_abort(fpo, "(read_xyz_geom): Trouble reading data\n");
    junk = label2an(alabel);
    if (junk == -1) 
      prf_abort(fpo, "(read_xyz_geom): Trouble parsing atom label\n");
    (*AN)[i] = junk;
    i++; 
  }

  if (*natom != i) {
    fprintf(fpo, "(read_xyz_geom): Should be %d atoms but read %d!\n",
            *natom, i);
    fprintf(fpo, "Extra atoms will be ignored.\n");
  }
  fclose(fpi);
  return(1);
}

}} // namespace psi::geom
