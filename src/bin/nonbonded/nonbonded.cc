/*
** NONBONDED
**
** Evaluate nonbonded interactions such as van der Waals terms and
** point-charge electrostatics
**
** C. David Sherrill
** Edward G. Hohenstein
**
** January 2008: Initial version
** November 2008: Hessian of -D term
** 
*/


/*! 
** \file
** \ingroup NONBONDED
** \brief Evaluate empirical non-bonded interactions
*/

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cctype>               // for toupper()
#include <libciomr/libciomr.h>
#include <libipv1/ip_lib.h>
#include <libchkpt/chkpt.h>
#include <libpsio/psio.h>
#include <psifiles.h>
#include <cstring>
#include <physconst.h>
#include <psifiles.h>
#include <masses.h>
#include <string>

#include "globals.h"
#include "nonbonded.h"

namespace psi { namespace nonbonded {

/* function prototypes this module */
void start_io(int argc, char *argv[]);
void stop_io(void);
void print_intro(void);
void compute_R(int natom, double **geom, double *R);
void compute_dXdYdZ(int natom, double **geom, double **dX, double **dY,
  double **dZ);
double compute_ddisp(int natom, double *R, double *AN, double s6, double d,
  char *damp_type);
void compute_ddisp_gradient(int natom, double *AN, double **geom,
  double s6, double d, char *damp_type, double *gradient);
void compute_disp_hess(int natom, double *AN, double **geom,
  double s6, double **hessian);
void compute_ddisp_hess(int natom, double *AN, double **geom,
  double s6, double d, double **hessian);
int label2an(char *label);
double compute_estatic(int natom, double *R, double *AC);

}} // close namespace decl

using namespace psi::nonbonded;

int main(int argc, char *argv[]) {

  double **geom;                  /* geometry matrix (cols are x,y,z)   */
  double *AC;                     /* atomic charges                     */
  double *AN;                     /* atomic numbers                     */
  double *R;                      /* distance matrix (lwr triangle)     */
  double **dX, **dY, **dZ;        /* x, y, z differences (lwr triangle) */
  int natom;                      /* number of atoms                    */
  int have_partial_charges=0;     /* flag for partial charges available */
  double s6=1.0;                  /* global scaling parameter (Grimme),
                                     0.75 PBE, 1.2 BLYP, 1.05 B-P86,
                                     1.0 TPSS, 1.05 B3LYP               */
  double d=20.0;                  /* exponent for damping term (Grimme) */
  double energy_dd;               /* damped dispersion energy (J/mol)   */
  double energy_dd_hartree;       /* in hartree                         */
  double e_scf;                   /* Hartree-Fock energy (hartree)      */
  double etot;                    /* prev tot energy+energy_dd_hartree  */
  int do_dispersion;              /* compute Grimme dispersion terms?   */
  int do_estatic;                 /* compute electrostatic terms?       */
  char *damp_type;                /* which damping equation to be used? */
  char *disp_method;              /* which C6 function to be used       */
  int num_array;                  /* count elements in parsed array     */
  int errcod;                     /* input parsing error code           */
  double tval;                    /* temp var                           */
  char tmpstr[100];               /* temp string for parsing input      */
  char *tmpstr2;                  /* another temporary string           */
  char *dertype;                  /* derivative type from input         */
  double *gradient;               /* gradient vector                    */
  double **hessian;               /* hessian matrix                     */
  int derlvl;                     /* do energy(0), deriv(1), or Hess(2)?*/
  int i,j;

  start_io(argc,argv);
  print_intro();
  geom = chkpt_rd_geom();
  natom = chkpt_rd_natom();
  AN = chkpt_rd_zvals();

  do_estatic = 0;
  if (ip_exist("PARTIAL_CHARGES",0)) {
    AC = init_array(natom);
    errcod = ip_double_array("PARTIAL_CHARGES",AC,natom);
    have_partial_charges = 1;
    do_estatic = 1;
  }

  errcod = ip_data("S6","%lf",&s6,0);
  errcod = ip_data("D_DMP","%lf",&d,0);

  do_dispersion = 1;
  damp_type = (char *) malloc(80);
  if (ip_exist("DISPERSION",0)) {
    errcod = ip_string("DISPERSION",&disp_method,0);   
    if (strcmp(disp_method, "GRIMME")==0) {
      strcpy(damp_type, "GRIMME");
    }
    else if (strcmp(disp_method, "NONE")==0 ||
             strcmp(disp_method, "FALSE")==0) {
      do_dispersion = 0;
      strcpy(damp_type, "NONE");
    }
    else {
      fprintf(outfile, "Error: unrecognized dispersion type %s\n", disp_method);
      do_dispersion = 0;
    }
  }
  else { // no DISPERSION keyword in input, so set default
    disp_method = (char *) malloc(80);
    strcpy(disp_method, "GRIMME");
    strcpy(damp_type, "GRIMME");
  }

  if (ip_exist("DAMP", 0)) {
    errcod = ip_string("DAMP",&tmpstr2,0);
    strcpy(damp_type, tmpstr2);
    free(tmpstr2);
    if (strcmp(damp_type, "TRUE")==0)
      strcpy(damp_type, "GRIMME");
    else if (strcmp(damp_type, "FALSE")==0)
      strcpy(damp_type, "NONE");

    if (strcmp(damp_type, "GRIMME")!=0 &&
        strcmp(damp_type, "NONE")!=0)
      fprintf(outfile, "Error: unrecognized option for DAMP\n");
      exit(PSI_RETURN_FAILURE);
  }

  errcod = ip_boolean("ELECTROSTATICS",&(do_estatic),0);

  /* parse overridden vdW radii */
  if (ip_exist("VDW_RADII",0)) {
    errcod = ip_count("VDW_RADII", &num_array, 0);
    if (errcod != IPE_OK) {
      fprintf(stderr,"ERROR: %s\n", ip_error_message(errcod));
    }  
    else {
      for (i=0; i<num_array; i+=2) {
        errcod = ip_data("VDW_RADII", "%s", tmpstr, 1, i); 
        errcod = ip_data("VDW_RADII", "%lf", &tval, 1, i+1); 
        if ((j = label2an(tmpstr)) != -1) {
          vdw_radii_grimme[j] = tval;
          fprintf(outfile, 
            "   Using custom value of %6.3lf for vdW radius of element %s\n",
            tval, tmpstr);
        }
      }
    }
  } /* end parsing overridden vdW radii */

  /* parse overridden C6 coefficients */
  if (ip_exist("C6",0)) {
    errcod = ip_count("C6", &num_array, 0);
    if (errcod != IPE_OK) {
      fprintf(stderr,"ERROR: %s\n", ip_error_message(errcod));
    }  
    else {
      for (i=0; i<num_array; i+=2) {
        errcod = ip_data("C6", "%s", tmpstr, 1, i); 
        errcod = ip_data("C6", "%lf", &tval, 1, i+1); 
        if ((j = label2an(tmpstr)) != -1) {
          vdw_C6_grimme[j] = tval;
          fprintf(outfile, 
            "   Using custom value of %6.3lf for C6 for element %s\n",
            tval, tmpstr);
        }
      }
    }
  } /* end parsing overridden C6 coefficients */
 
  fprintf(outfile, "\n");
  if (do_dispersion) {
    fprintf(outfile, "   Universal scaling coefficient s6 = %6.4lf\n", s6);
    fprintf(outfile, "   Universal damping exponent d     = %6.4lf\n", d);
    fprintf(outfile, "   Dispersion method                = %s\n", 
      disp_method);
    fprintf(outfile, "   Damping function                 = %s\n", 
      damp_type);
  }

  /* get checkpoint energy and print it */
  if (chkpt_exist_add_prefix("Total energy")) {
    etot = chkpt_rd_etot();
    fprintf(outfile, "\n");
    fprintf(outfile, "   Previous energy           = %14.9lf hartree\n", etot);
  }
  else etot = 0.0;

  R = init_array((natom*(natom+1))/2);
  compute_R(natom, geom, R);

  if (do_dispersion) {

    energy_dd = compute_ddisp(natom, R, AN, s6, d, damp_type);
    energy_dd_hartree = energy_dd / (_na * _hartree2J);

    fprintf(outfile, "\n");
    fprintf(outfile, "   Damped dispersion energy  = %14.9lf hartree ", 
      energy_dd_hartree);
    fprintf(outfile, "(%10.4lf kcal/mol)\n", (energy_dd / 4184.0));
    fflush(outfile);
    etot += energy_dd_hartree;
  }

  if (do_estatic) {
    double estatic = compute_estatic(natom, R, AC);
    double estatic_hartree = estatic / (_na * _hartree2J);
    fprintf(outfile, "   Electrostatic energy      = %14.9lf hartree ", 
      estatic_hartree);
    fprintf(outfile, "(%10.4lf kcal/mol)\n", (estatic / 4184.0));
    etot += estatic_hartree;
  }

  fprintf(outfile, " * Total energy + empirical  = %14.9lf hartree\n", etot);
  fprintf(outfile, "\n");
  chkpt_wt_etot(etot);

  /* calculate the gradient if requested */
  errcod = ip_string("DERTYPE", &dertype, 0);
  if (errcod == IPE_KEY_NOT_FOUND) {
     dertype = (char *) malloc(sizeof(char)*5);
     strcpy(dertype, "NONE");
   }
  if (strcmp(dertype, "NONE")==0) derlvl = 0;
  else if (strcmp(dertype, "FIRST")==0) derlvl = 1;
  else if (strcmp(dertype, "SECOND")==0) derlvl = 2;
  else derlvl = 0; 

  if (derlvl==1) { // do a gradient

    /* this should never happen if called before $deriv
    if (chkpt_exist_add_prefix("Energy Gradient")) {
      gradient = chkpt_rd_grad();
      fprintf(outfile, "  Energy gradient from checkpoint (hartree/bohr)\n");
      for (i=0; i<natom; i++) {
        fprintf(outfile, "  %12.9lf  %12.9lf  %12.9lf\n",
          gradient[i*3], gradient[i*3+1], gradient[i*3+2]);
      }
      fprintf(outfile, "\n");
    }
    else {
      gradient = init_array(3*natom);
    }
    */

    gradient = init_array(3*natom);
    compute_ddisp_gradient(natom, AN, geom, s6, d, damp_type, gradient);

    fprintf(outfile, "  Gradient of empirical contribution (hartree/bohr)\n");
    for (i=0; i<natom; i++) {
      fprintf(outfile, "  %12.9lf  %12.9lf  %12.9lf\n",
        gradient[i*3], gradient[i*3+1], gradient[i*3+2]);
    }
    fprintf(outfile, "\n");
    chkpt_wt_grad(gradient);
  }

  if (derlvl==2) { // Hessian

    hessian = block_matrix(3*natom,3*natom);

    psio_open(PSIF_DERINFO, PSIO_OPEN_OLD);
    psio_read_entry(PSIF_DERINFO, "Skeleton Hessian",
               (char *) hessian[0], natom*3*natom*3*sizeof(double));
    psio_close(PSIF_DERINFO, 1);

    if (strcmp(damp_type, "GRIMME")==0) {
      compute_ddisp_hess(natom, AN, geom, s6, d, hessian);
    }
    else if (strcmp(damp_type, "NONE")==0) {
      compute_disp_hess(natom, AN, geom, s6, hessian);
    }
    else {
      fprintf(outfile, "Error: unrecognized damp type for Hessian!\n");
    }

  }

  if (do_estatic && derlvl > 0) {
    fprintf(outfile, "   Warning: Electrostatic interactions not ");
    fprintf(outfile, "included in derivatives!!\n");
  }

  /* clean up */
  free(AN);
  free(R);
  if (have_partial_charges) free(AC);
  if (derlvl==1) free(gradient);
  if (derlvl==2) free_block(hessian);
  free_block(geom);
  stop_io();
  return(PSI_RETURN_SUCCESS);
}

extern "C" { const char *gprgid(void) { const char *prgid = "NONBONDED"; return (prgid); } }

namespace psi { namespace nonbonded {

/*!
** start_io(): Initiate PSI IO routines
**
** \param argc = number of command-line arguments
** \param argv = the command-line arguments
**
** \ingroup NONBONDED
*/
void start_io(int argc, char *argv[])
{
  int errcod;
  
  errcod = psi_start(&infile,&outfile,&psi_file_prefix,argc-1,argv+1,0);
  if (errcod != PSI_RETURN_SUCCESS)
    abort();
  ip_cwk_add(":NONBONDED");
  tstart(outfile);
  psio_init(); psio_ipv1_config();
  chkpt_init(PSIO_OPEN_OLD);
   
  return;
}


/*!
** stop_io(): Shut down PSI IO
** 
** Returns: none
**
** \ingroup NONBONDED
*/
void stop_io(void)
{
  tstop(outfile);
  psio_done();
  psi_stop(infile,outfile,psi_file_prefix);
}
 


/*!
** print_intro(): Print an intro for the module
**
** Returns: none
**
** \ingroup NONBONDED
*/
void print_intro(void)
{
 fprintf(outfile,"             ---------------------------------------\n");
 fprintf(outfile,"                            NONBONDED                  \n");
 fprintf(outfile,"               Evaluate empirical non-bonded terms     \n");
 fprintf(outfile,"                        C. David Sherrill              \n");
 fprintf(outfile,"             ---------------------------------------\n");
 fprintf(outfile,"\n");
}


/*!
** compute_R(): Compute the distances between each pair of atoms.  
** Convert to Angstrom because those are the units Grimme's parameters 
** are in.
**
** \param natom = number of atoms
** \param geom  = geometry matrix (cols are x, y, z; rows are atoms) 
** \param R     = matrix of distances (packed lower triangle, ij indexing)
**
** Returns: none
**
** \ingroup NONBONDED
*/
void compute_R(int natom, double **geom, double *R)
{
  int i, j, ij;
  double dx, dy, dz, dist;

  for (i=0,ij=0; i<natom; i++) {
    for (j=0; j<i; j++,ij++) {
      dx = geom[i][0] - geom[j][0];
      dy = geom[i][1] - geom[j][1];
      dz = geom[i][2] - geom[j][2];
      dist = sqrt(dx * dx + dy * dy + dz * dz);
      R[ij] = dist * _bohr2angstroms;
    }
  }
 
}


/*!
** compute_dXdYdZ(): Compute the X, Y, and Z differences between each pair 
** of atoms.  Need this for the gradients.
**
** \param natom = number of atoms
** \param geom  = geometry matrix (cols are x, y, z; rows are atoms)
** \param dX    = matrix of X distances
** \param dY    = matrix of Y distances
** \param dZ    = matrix of Z distances
**
** Returns: none
**
** \ingroup NONBONDED
*/
void compute_dXdYdZ(int natom, double **geom, double **dX, double **dY,
  double **dZ)
{
  int i, j;

  for (i=0; i<natom; i++) {
    for (j=0; j<natom; j++) {
      dX[i][j] = (geom[i][0] - geom[j][0]) * _bohr2angstroms;
      dY[i][j] = (geom[i][1] - geom[j][1]) * _bohr2angstroms;
      dZ[i][j] = (geom[i][2] - geom[j][2]) * _bohr2angstroms;
    }
  }

}


/*!
** compute_ddisp(): Compute damped dispersion terms 
** (ala S. Grimme, J. Comput. Chem. 27, 1787, 2006)
**
** \param natom      = number of atoms
** \param R          = matrix of interatomic distances (packed lower triangle)
** \param AN         = atomic number for each atom
** \param s6         = s6 universal dispersion scaling parameter
** \param d          = d universal dispersion damping exponent
** \param damp_type  = string giving damping function to use
**
** Returns: the dispersion energy (J/mol), damped if requested
**
** \ingroup NONBONDED
*/
double compute_ddisp(int natom, double *R, double *AN, double s6, double d,
  char *damp_type)
{
  int i, j, ij, Zi, Zj;
  double energy=0.0, tval, r, r_i, r_j, r_vdw;
  double C6i, C6j, fdmp; 

  /* loop over unique pairs of atoms and evaluate the C6 dispersion terms */
  for (i=0,ij=0; i<natom; i++) {
    Zi = (int) AN[i];
    C6i = vdw_C6_grimme[Zi];
    r_i = vdw_radii_grimme[Zi];
    for (j=0; j<i; j++,ij++) {
      Zj = (int) AN[j];
      r = R[ij];
      C6j = vdw_C6_grimme[Zj];
      r_j = vdw_radii_grimme[Zj];
      r_vdw = r_i + r_j;
      tval = sqrt(C6i * C6j) * 1000000.0; /* nm^6 -> Ang^6 */
      tval = tval/pow(r, 6.0);
      // printf("Undamped Dispersion term: %14.9lf\n", tval / 4184);
      if (strcmp(damp_type, "GRIMME") == 0) {
        fdmp = 1.0 / (1.0 + exp(-d * (r / r_vdw - 1.0))); 
        // printf("Damping factor          : %14.9lf\n", fdmp);
        tval *= fdmp;
      }
      energy += tval;
    }
  } 

  energy = -s6 * energy; /* in J mol^-1 */
  return energy;
}


/*!
** compute_ddisp_gradient(): Compute the gradient of the damped dispersion 
** terms from S. Grimme, J. Comput. Chem. 27, 1787, 2006
**
** \param natom      = number of atoms
** \param AN         = atomic number for each atom
** \param geom       = geometry matrix (in bohr)
** \param s6         = s6 universal dispersion scaling parameter
** \param d          = d universal dispersion damping exponent
** \param damp_type  = string giving damping type
** \param gradient   = gradient vector to add the results to (x1, y1, z1, ..)
**
** Returns: None
**
** Don't forget to convert to hartree/bohr
**
** \ingroup NONBONDED
*/
void compute_ddisp_gradient(int natom, double *AN, double **geom,
  double s6, double d, char *damp_type, double *gradient)
{
  int i, j, ij, Zi, Zj;
  double energy=0.0, tval, r, r_i, r_j, r_vdw;
  double C6i, C6j, C6, fdmp, qdmp, prefact, unitconv;
  double dx, dy, dz;

  /* loop over unique pairs of atoms and evaluate gradient */
  for (i=0; i<natom; i++) {
    Zi = (int) AN[i];
    C6i = vdw_C6_grimme[Zi];
    r_i = vdw_radii_grimme[Zi];
    for (j=0; j<natom; j++) {
      if (j==i) continue;

      dx = (geom[i][0] - geom[j][0]) * _bohr2angstroms;
      dy = (geom[i][1] - geom[j][1]) * _bohr2angstroms;
      dz = (geom[i][2] - geom[j][2]) * _bohr2angstroms;
      r =  sqrt(dx * dx + dy * dy + dz * dz); /* already ang now */

      Zj = (int) AN[j];
      C6j = vdw_C6_grimme[Zj];
      r_j = vdw_radii_grimme[Zj];
      r_vdw = r_i + r_j;
      prefact = s6 * sqrt(C6i * C6j) * 1000000.0; /* nm^6 -> Ang^6 */
      tval = 6.0/pow(r, 8.0);

      if (strcmp(damp_type, "GRIMME") == 0) {
        qdmp = exp(-d * (r / r_vdw - 1.0));
        fdmp = 1.0 / (1.0 + qdmp);
        tval *= fdmp;

        tval -= (d / (r_vdw * pow(r, 7.0))) * (qdmp * fdmp * fdmp);
      }

      tval *= prefact;

      /* right now we're in J / (mol * Angstrom) */
      unitconv = _bohr2angstroms / (_na * _hartree2J);
      tval *= unitconv;

      /* x */
      gradient[i*3] += tval * dx;

      /* y */
      gradient[i*3+1] += tval * dy;
  
      /* z */
      gradient[i*3+2] += tval * dz;

    } /* end loop over j */
  } /* end loop over i */

}


/*!
** compute_disp_hess(): Compute the Hessian of the undamped dispersion
** terms from S. Grimme, J. Comput. Chem. 27, 1787, 2006
**
** \param natom      = number of atoms
** \param AN         = atomic number for each atom
** \param geom       = geometry matrix (in bohr)
** \param s6         = s6 universal dispersion scaling parameter
** \param hessian    = Hessian to add the results to 
**
** Returns: none
**
** \ingroup NONBONDED
*/
void compute_disp_hess(int natom, double *AN, double **geom,
  double s6, double **hessian)
{
  int i,j,k;
  int Zi,Zj,Zk;
  double r;
  double C6i,C6j,C6k;
  double dx, dy, dz;
  double tval_10, tval_8;
  double prefact;
  double unitconv = _bohr2angstroms * _bohr2angstroms / (_na * _hartree2J);

  for(i=0; i < natom; i++) {
    for(j=0; j <= i; j++) { /*Loop over all unique atom pairs (including i=j)*/

      Zi = (int) AN[i];
      Zj = (int) AN[j];

      C6i = vdw_C6_grimme[Zi];
      C6j = vdw_C6_grimme[Zj];

      if (i==j) { /* Handle cases where i=j first */

        for(k=0; k < natom; ++k){
          if (i == k) continue;

          Zk = (int) AN[k];
          C6k = vdw_C6_grimme[Zk];

          prefact = s6 * sqrt(C6i * C6k) * 1000000.0 * unitconv;

          dx = (geom[i][0] - geom[k][0]) * _bohr2angstroms;
          dy = (geom[i][1] - geom[k][1]) * _bohr2angstroms;
          dz = (geom[i][2] - geom[k][2]) * _bohr2angstroms;
          r =  sqrt(dx * dx + dy * dy + dz * dz); /* in ang */

          tval_8  = 6.0/pow(r, 8.0);
          tval_10 = -48.0/pow(r, 10.0);

          hessian[3*i][3*i] += prefact * ( tval_10 * dx * dx + tval_8);
          hessian[3*i+1][3*i+1] += prefact * ( tval_10 * dy * dy + tval_8);
          hessian[3*i+2][3*i+2] += prefact * ( tval_10 * dz * dz + tval_8);
          hessian[3*i+1][3*i] += prefact * tval_10 * dx * dy;
          hessian[3*i+2][3*i] += prefact * tval_10 * dx * dz;
          hessian[3*i+2][3*i+1] += prefact * tval_10 * dy * dz;

        }
      }
      else { /* Now i != j */

        prefact = s6 * sqrt(C6i * C6j) * 1000000.0 * unitconv;

        dx = (geom[i][0] - geom[j][0]) * _bohr2angstroms;
        dy = (geom[i][1] - geom[j][1]) * _bohr2angstroms;
        dz = (geom[i][2] - geom[j][2]) * _bohr2angstroms;
        r =  sqrt(dx * dx + dy * dy + dz * dz); /* in ang */

        tval_8  = -6.0/pow(r, 8.0);
        tval_10 = 48.0/pow(r, 10.0); /* Sign change from xij = - xji */

        hessian[3*i][3*j] += prefact * ( tval_10 * dx * dx + tval_8);
        hessian[3*i+1][3*j+1] += prefact * ( tval_10 * dy * dy + tval_8);
        hessian[3*i+2][3*j+2] += prefact * ( tval_10 * dz * dz + tval_8);
        hessian[3*i+1][3*j] += prefact * tval_10 * dx * dy;
        hessian[3*i+2][3*j] += prefact * tval_10 * dx * dz;
        hessian[3*i+2][3*j+1] += prefact * tval_10 * dy * dz;
        hessian[3*i][3*j+1] += prefact * tval_10 * dx * dy;
        hessian[3*i][3*j+2] += prefact * tval_10 * dx * dz;
        hessian[3*i+1][3*j+2] += prefact * tval_10 * dy * dz;

      } // end i!=j
    } // end loop over atoms j
  } // end loop over atoms i

  for(i=0; i < 3*natom; i++) { // symmetrize Hessian
    for(j=0; j < i; j++){
      hessian[j][i] = hessian[i][j];
    }
  }

  fprintf(outfile,
    "   Hessian contribution from undamped empirical correction\n");
  print_mat(hessian,3*natom,3*natom,outfile);

  psio_open(PSIF_DERINFO, PSIO_OPEN_NEW);
  psio_write_entry(PSIF_DERINFO, "Skeleton Hessian", (char *) hessian[0], 
    natom*3*natom*3*sizeof(double));
  psio_close(PSIF_DERINFO, 1);
}


/*!
** compute_ddisp_hess(): Compute the Hessian of the damped dispersion 
** terms from S. Grimme, J. Comput. Chem. 27, 1787, 2006
**
** \param natom      = number of atoms
** \param AN         = atomic number for each atom
** \param geom       = geometry matrix (in bohr)
** \param s6         = s6 universal dispersion scaling parameter
** \param d          = d universal dispersion damping exponent
** \param hessian    = Hessian to add the results to
**
** Returns: none
**
** \ingroup NONBONDED
*/
void compute_ddisp_hess(int natom, double *AN, double **geom,
  double s6, double d, double **hessian)
{
  int i,j,k;
  int Zi,Zj,Zk;
  double r;
  double C6i,C6j,C6k;
  double dx, dy, dz;
  double dmp;
  double r_i, r_j, r_k, rvdw;
  double term1, term2, term3, term4, term5, term6;
  double clump1, clump2;
  double prefact;
  double unitconv = _bohr2angstroms * _bohr2angstroms / (_na * _hartree2J);

  for(i=0; i < natom; i++){
    for(j=0; j <= i; j++){ /*Loop over all unique atom pairs (including i=j)*/

      Zi = (int) AN[i];
      Zj = (int) AN[j];

      C6i = vdw_C6_grimme[Zi];
      C6j = vdw_C6_grimme[Zj];

      r_i = vdw_radii_grimme[Zi];
      r_j = vdw_radii_grimme[Zj];

      if (i==j) { /* Handle cases where i=j first */

        for(k=0; k < natom; ++k){
          if (i == k) continue;

          Zk = (int) AN[k];
          C6k = vdw_C6_grimme[Zk];
          r_k = vdw_radii_grimme[Zk];

          prefact = s6 * sqrt(C6i * C6k) * 1000000.0 * unitconv;
          rvdw = r_i + r_k;

          dx = (geom[i][0] - geom[k][0]) * _bohr2angstroms;
          dy = (geom[i][1] - geom[k][1]) * _bohr2angstroms;
          dz = (geom[i][2] - geom[k][2]) * _bohr2angstroms;
          r =  sqrt(dx * dx + dy * dy + dz * dz); /* in ang */

          dmp = exp(-d * (r / rvdw - 1.0));

          term1 = -48.0 / (pow(r, 10.0) * ( 1 + dmp ));
          term2 = 13.0 * d * dmp / (pow(r, 9.0) * pow(1 + dmp, 2.0) * rvdw);
          term3 = 6.0 / (pow(r, 8.0) * ( 1 + dmp ));
          term4 = -2.0*d*d*dmp*dmp / (pow(r, 8.0) * pow(1 + dmp, 3.0) * rvdw * rvdw );
          term5 = -d*dmp/(pow(r, 7.0) * pow(1 + dmp, 2.0) * rvdw );
          term6 = d*d*dmp/(pow(r, 8.0) * pow(1 + dmp, 2.0) * rvdw * rvdw );

          clump1 = term1 + term2 + term4 + term6;
          clump2 = term3 + term5;

          hessian[3*i][3*i] += prefact * ( dx * dx * clump1 + clump2 );
          hessian[3*i+1][3*i+1] += prefact * ( dy * dy * clump1 + clump2 );
          hessian[3*i+2][3*i+2] += prefact * ( dz * dz * clump1 + clump2 );
          hessian[3*i+1][3*i] += prefact * clump1 * dx * dy;
          hessian[3*i+2][3*i] += prefact * clump1 * dx * dz;
          hessian[3*i+2][3*i+1] += prefact * clump1 * dy * dz;

        }
      } // end i==j

      else { // i!=j

        prefact = s6 * sqrt(C6i * C6j) * 1000000.0 * unitconv;
        rvdw = r_i + r_j;

        dx = (geom[i][0] - geom[j][0]) * _bohr2angstroms;
        dy = (geom[i][1] - geom[j][1]) * _bohr2angstroms;
        dz = (geom[i][2] - geom[j][2]) * _bohr2angstroms;
        r =  sqrt(dx * dx + dy * dy + dz * dz); /* in ang */

        dmp = exp(-d * (r / rvdw - 1.0));

        term1 = 48.0 / (pow(r, 10.0) * ( 1 + dmp ));
        term2 = -13.0 * d * dmp / (pow(r, 9.0) * pow(1 + dmp, 2.0) * rvdw);
        term3 = -6.0 / (pow(r, 8.0) * ( 1 + dmp ));
        term4 = 2.0*d*d*dmp*dmp / (pow(r, 8.0) * pow(1 + dmp, 3.0) * 
          rvdw * rvdw );
        term5 = d*dmp/(pow(r, 7.0) * pow(1 + dmp, 2.0) * rvdw );
        term6 = -d*d*dmp/(pow(r, 8.0) * pow(1 + dmp, 2.0) * rvdw * rvdw );

        clump1 = term1 + term2 + term4 + term6;
        clump2 = term3 + term5;

        hessian[3*i][3*j] += prefact * ( dx * dx * clump1 + clump2 );
        hessian[3*i+1][3*j+1] += prefact * ( dy * dy * clump1 + clump2 );
        hessian[3*i+2][3*j+2] += prefact * ( dz * dz * clump1 + clump2 );
        hessian[3*i+1][3*j] += prefact * clump1 * dx * dy;
        hessian[3*i+2][3*j] += prefact * clump1 * dx * dz;
        hessian[3*i+2][3*j+1] += prefact * clump1 * dy * dz;
        hessian[3*i][3*j+1] += prefact * clump1 * dx * dy;
        hessian[3*i][3*j+2] += prefact * clump1 * dx * dz;
        hessian[3*i+1][3*j+2] += prefact * clump1 * dy * dz;
      }

    } // end loop over j
  } // end loop over i

  for(i=0; i < 3*natom; ++i){
    for(j=0; j < i; ++j){
      hessian[j][i] = hessian[i][j];
    }
  }

  fprintf(outfile,"   Hessian contribution from damped empirical correction\n");
  print_mat(hessian,3*natom,3*natom,outfile);
  fflush(outfile);

  psio_open(PSIF_DERINFO, PSIO_OPEN_NEW);
  psio_write_entry(PSIF_DERINFO, "Skeleton Hessian", (char *) hessian[0], 
    natom*3*natom*3*sizeof(double));
  psio_close(PSIF_DERINFO, 1);
}


/*!
** compute_estatic(): This function computes the electrostatic contribution
** to the energy, assuming each atom has a partial charge centered on
** the atomic center.
** 
** \param natom = number of atoms
** \param R     = lower triangle array of distances
** \param AC    = array of atomic charges
**
** Returns: electrostatic energy (J/mol)
** Based on code by Michael Marshall, 2008
**
** \ingroup NONBONDED
*/
double compute_estatic(int natom, double *R, double *AC) 
{

  int i, j, ij;
  double q1, q2;
  double estatic = 0.0;
  double au_to_coulomb = 1.60219E-19;
  double convfact;

  convfact = au_to_coulomb * au_to_coulomb * _na / (4.0 * _pi * _e0);
  convfact *= 1.0E10; /* Angstroms to meters in denominator */

  for(i=0,ij=0; i<natom; i++){
    q1=AC[i];
    for(j=0; j<i; j++,ij++){
       q2=AC[j];
       estatic += q1*q2/R[ij];
    }
  }
  // fprintf(outfile, "estatic before conversion = %lf\n", estatic);
  estatic *= convfact;
  // fprintf(outfile, "estatic after conversion  = %lf\n", estatic);
  return estatic;
}


/*!
** label2an(): This function returns the atomic number corresponding to a 
** given mass label, or -1 if not found.
**
** \param label = Atom label we're trying to match (case insensitive)
**
** Returns: atomic number of matched label (or -1 if not found)
**
** Depends on atomic_labels in masses.h
**
** C. David Sherrill
** July 1999
**
** \ingroup NONBONDED
*/
int label2an(char *label)
{
  int i, j, k;
  char p, q;

  for (i=0; i<LAST_ATOMIC_INDEX; i++) {
    k = strlen(label);
    for (j=0; j < k; j++) {
       p = label[j];
       q = atomic_labels[i][j];
       if (toupper(p) != toupper(q)) break;
    }
    if (j == k) return(i);
  }

  /* couldn't find the label */
  return(-1);

} // end label2an()


}} // close off psi::nonbonded
