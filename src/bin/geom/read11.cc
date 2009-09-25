/*! \file
    \ingroup GEOM
    \brief Enter brief description of file here 
*/
/*
** NEW VERSION OF READ_FILE11() TO READ PAST ALL THE OLD ENTRIES ON TOP
**
** David Sherrill  -  June 1993
**
*/

/* include's */
#include <cstdio>
#include <cstdlib>

/* define's */
#define MAX_LINE 132

namespace psi { namespace geom {

extern void malloc_ck(void *array, const char *mesg);
extern void prf_abort(FILE *file, const char *mesg);

/* global variables */
static char line1[MAX_LINE] ;

/*
** READ_FILE11(): This function reads 'file11.dat' and prints out again
**      some of the info to an output file.
**
** Arguments: 
**      label   = character string to hold file11 label field
**      natom   = pointer to variable to hold number of atoms
**      energy  = pointer to var to hold energy
**      X, Y, Z = pointers to arrays of cartesian coordinates (assume bohr)
**                (these are currently allocated in THIS function)
**      AN      = pointer to atomic number array (allocated here)
**      fpo     = file pointer for output
**
** Returns: success (1) or failure (0) 
**
*/
int read_file11(char* label, int* natom, double* energy, double** X, double** Y, double** Z, double** AN, FILE* fpo)
{
int i ;                  /* loop variable */
FILE *fpi ;              /* pointer for input file file11.dat */
double g1, g2, g3;       /* junk variables for gradients */

/* open file 11 for reading */
   if ((fpi = fopen("file11.dat", "r")) == NULL) {
      printf("(read_file11): Error opening file11.dat\n");
      exit(0);
      }

/* get first file11 label */
   if (fgets(label, MAX_LINE, fpi) == NULL)
      prf_abort(fpo, "(read_file11): Trouble reading label\n") ;
   fflush(fpo) ;

/* now read the number of atoms and the energy */
   fgets(line1, MAX_LINE, fpi) ;
   if (sscanf(line1, "%d %lf", natom, energy) != 2) 
      prf_abort(fpo, "(read_file11): Trouble reading natoms and energy\n") ;
   
/* now make room for the Cartesian coordinates */
   *X = (double *) malloc (*natom * sizeof(double)) ;
   malloc_ck((void*)*X, "(read_file11): Trouble allocating Cartesian array\n") ;
   *Y = (double *) malloc (*natom * sizeof(double)) ;
   malloc_ck((void*)*Y, "(read_file11): Trouble allocating Cartesian array\n") ;
   *Z = (double *) malloc (*natom * sizeof(double)) ;
   malloc_ck((void*)*Z, "(read_file11): Trouble allocating Cartesian array\n") ;
   *AN = (double *) malloc (*natom * sizeof(double)) ;
   malloc_ck((void*)*AN, "(read_file11): Trouble allocating atomic num array\n") ;

/* now rewind file11 and read it one chunk at a time */
   rewind(fpi) ;

/* loop over each chunk */
/* read number of atoms and energy */
   while (fgets(label, MAX_LINE, fpi) != NULL) {
      fgets(line1, MAX_LINE, fpi) ;
      if (sscanf(line1, "%d %lf", natom, energy) != 2) 
         prf_abort(fpo, "(read_file11): Trouble reading natoms and energy\n") ;

/* read Cartesians */
      for (i=0; i<(*natom); i++) {
         if (fscanf(fpi, "%lf %lf %lf %lf", (*AN+i), (*X+i), (*Y+i), 
            (*Z+i)) != 4)
            prf_abort(fpo, "(read_file11): Trouble reading Cartesian coords\n");
            }

/* read Gradients */
      for (i=0; i<(*natom); i++) {
         if (fscanf(fpi, "%lf %lf %lf", &g1, &g2, &g3) != 3)
            prf_abort(fpo, "(read_file11): Trouble reading gradients\n");
            }
      fgets(line1, MAX_LINE, fpi) ;   /* go to end of line */
      } /* end loop over file11 chunk */

   fclose(fpi) ;
   return(1) ;
}



/*
** PRINT_FILE11(): Function prints out the information from file11.dat
** obtained from the read_file11() function.
**
** Arguments: 
**      label   = character string to hold file11 label field
**      natom   = number of atoms
**      energy  = energy from file11
**      X, Y, Z = arrays of cartesian coordinates (assume bohr)
**      AN      = atomic number array (allocated here)
**      fpo     = file pointer for output
*/
void print_file11(char* label, int natom, double energy, double* X, double* Y, double* Z, double* AN, FILE* fpo) 
{
int i ;

   fprintf(fpo, "DATA FROM FILE11.DAT\n") ;
   fprintf(fpo, "Label :\n%s\n", label) ;
   fprintf(fpo, "Number of atoms = %d\n", natom) ;
   fprintf(fpo, "Energy = %.10lf\n", energy) ;
   fprintf(fpo, "\nCartesian coordinates (bohr) :\n") ;
   for (i=0; i<natom; i++) {
      fprintf(fpo, "     %11.6lf    %12.7lf    %12.7lf    %12.7lf\n",
            AN[i], X[i], Y[i], Z[i]) ;
      }
   fprintf(fpo, "\n") ;
   fflush(fpo) ;
}

}} // namespace psi::geom
