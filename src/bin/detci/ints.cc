/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

/*! \file
    \ingroup DETCI
    \brief Enter brief description of file here
*/


/*
** INTS.C
**
** Return values of one and two-electron integrals
**
** C. David Sherrill
** Center for Computational Quantum Chemistry
** University of Georgia
**
** Updated 3/18/95 to exclude frozen virtual orbitals.
** Updated 3/28/95 to exclude frozen core orbitals.
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <libiwl/iwl.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <psifiles.h>
#include <libtrans/integraltransform.h>
#include <libmints/mints.h>
#include <psi4-dec.h>
#include "structs.h"
#include "ciwave.h"
#define EXTERN
#include "globals.h"
#include "globaldefs.h"

namespace psi { namespace detci {

// #define MIN0(a,b) (((a)<(b)) ? (a) : (b))
// #define MAX0(a,b) (((a)>(b)) ? (a) : (b))

void CIWavefunction::transform_ci_integrals()
{
   extern double check_energy(double *H, double *twoel_ints, int *docc,
      int *dropped_docc, int drc_flag, double escf, double enuc, double edrc,
      int nirreps, int *reorder, int *opi, int print_lvl, std::string out);

  // Grab orbitals
  SharedMatrix drc = get_orbitals("DRC");
  SharedMatrix act = get_orbitals("ACT");
  SharedMatrix vir = get_orbitals("VIR");
  SharedMatrix fzv = get_orbitals("FZV");

  // Build up active space
  std::vector<boost::shared_ptr<MOSpace> > spaces;

  std::vector<int> orbitals(CalcInfo.num_ci_orbs, 0);
  std::vector<int> indices(CalcInfo.num_ci_orbs, 0);
  for (int h = 0, cinum = 0, orbnum = 0; h < CalcInfo.nirreps; h++){
    orbnum += CalcInfo.dropped_docc[h];
    for (int i = 0; i < CalcInfo.ci_orbs[h]; i++){
      orbitals[cinum] = orbnum;
      indices[cinum] = CalcInfo.reorder[orbnum] - CalcInfo.num_drc_orbs;
      orbnum++;
      cinum++;
    }
    orbnum += CalcInfo.dropped_uocc[h];
  }

  boost::shared_ptr<MOSpace> custom_space(new MOSpace('X', orbitals, indices));
  spaces.push_back(custom_space);

  IntegralTransform *ints = new IntegralTransform(drc, act, vir, fzv, spaces,
                            IntegralTransform::Restricted,
                            IntegralTransform::IWLOnly,
                            IntegralTransform::PitzerOrder,
                            IntegralTransform::OccAndVir,
                            true);

  // Incase we do two ci runs
  ints->set_keep_iwl_so_ints(true);
  ints->transform_tei(custom_space, custom_space, custom_space, custom_space);

  // Build useful information
  int nmotri_full = (CalcInfo.nmo * (CalcInfo.nmo + 1)) / 2 ;
  int nmotri = (CalcInfo.num_ci_orbs * (CalcInfo.num_ci_orbs + 1)) / 2 ;

  // Build desired arrays
  double*  tmp_onel_ints1 = (double *) init_array(nmotri_full);
  double*  tmp_onel_ints2 = (double *) init_array(nmotri_full);
  CalcInfo.onel_ints = (double *) init_array(nmotri);
  CalcInfo.twoel_ints = (double *) init_array(nmotri * (nmotri + 1) / 2);
  CalcInfo.maxK = (double *) init_array(CalcInfo.num_ci_orbs);


  // Read one electron integrals
  iwl_rdone(Parameters.oei_file, PSIF_MO_FZC, tmp_onel_ints1, nmotri_full,
               0, (Parameters.print_lvl>4), "outfile");

  // IntegralTransform does not properly order one electron integrals for whatever reason
  for (int i=0, cnt=0; i<CalcInfo.nmo; i++){
    for (int j=i; j<CalcInfo.nmo; j++){
      int reorder_idx = INDEX(CalcInfo.reorder[i],CalcInfo.reorder[j]);
      tmp_onel_ints2[reorder_idx] = tmp_onel_ints1[INDEX(i,j)];
    }
  }

  filter(tmp_onel_ints2, CalcInfo.onel_ints, ioff, CalcInfo.nmo,
  CalcInfo.num_drc_orbs, CalcInfo.num_drv_orbs);
  free(tmp_onel_ints1);
  free(tmp_onel_ints2);

  CalcInfo.edrc = ints->get_frozen_core_energy();
  iwl_rdtwo(PSIF_MO_TEI, CalcInfo.twoel_ints, ioff, CalcInfo.num_ci_orbs,
            0, 0, (Parameters.print_lvl>4), "outfile");

   /* Determine maximum K integral for use in averaging the diagonal */
   /* Hamiltonian matrix elements over spin-coupling set */
   if (Parameters.hd_ave) {
     for(int i=0; i<CalcInfo.num_ci_orbs; i++)
        for(int j=0; j<CalcInfo.num_ci_orbs; j++) {
           /* if (i==j) continue; */
           int ij = ioff[MAX0(i,j)] + MIN0(i,j);
           int ijij = ioff[ij] + ij;
           double value = CalcInfo.twoel_ints[ijij];
           if(value > CalcInfo.maxK[i]) CalcInfo.maxK[i] = value;
           }
      for(int i=0; i<CalcInfo.num_ci_orbs; i++) {
        if(CalcInfo.maxK[i] > CalcInfo.maxKlist)
          CalcInfo.maxKlist = CalcInfo.maxK[i];
        if (Parameters.print_lvl > 4)
          outfile->Printf("maxK[%d] = %lf\n",i, CalcInfo.maxK[i]);
        }
    }
 CalcInfo.eref = check_energy(CalcInfo.onel_ints, CalcInfo.twoel_ints,
    CalcInfo.docc, CalcInfo.dropped_docc, 1, CalcInfo.escf,
    CalcInfo.enuc, CalcInfo.edrc, CalcInfo.nirreps, CalcInfo.reorder,
    CalcInfo.orbs_per_irr, Parameters.print_lvl, "outfile");

  delete ints;
}
// void CIWavefunction::transform_dfmcscf_ints(){
//   boost::shared_ptr<BasisSet> primary = BasisSet::pyconstruct_orbital(
//     Process::environment.molecule(), "BASIS", options.get_str("BASIS"));
//   boost::shared_ptr<BasisSet> auxiliary = BasisSet::pyconstruct_auxiliary(primary->molecule(),
//       "DF_BASIS_SCF", options.get_str("DF_BASIS_SCF"), "JKFIT",
//       options.get_str("BASIS"), primary->has_puream());

//   /// Build JK, DFERI, and SOMCSCF objects
//   boost::shared_ptr<JK> jk = JK::build_JK();
//   jk->set_do_J(true);
//   jk->set_do_K(true);
//   jk->initialize();
//   jk->print_header();

//   boost::shared_ptr<DFERI> df = DFERI::build(primary,auxiliary,options);


// }


double get_onel(int i, int j)
{
   int ij ;
   double value ;

   if (i > j) {
      ij = ioff[i] + j;
      value = CalcInfo.onel_ints[ij] ;
      return(value) ;
      }
   else {
      ij = ioff[j] + i ;
      value = CalcInfo.onel_ints[ij] ;
      return(value) ;
      }
   return(CalcInfo.onel_ints[ij]) ;
}

double get_twoel(int i, int j, int k, int l)
{
   int ij, kl, ijkl ;

   ij = ioff[MAX0(i,j)] ;
   ij += MIN0(i,j) ;
   kl = ioff[MAX0(k,l)] ;
   kl += MIN0(k,l) ;
   ijkl = ioff[MAX0(ij,kl)] ;
   ijkl += MIN0(ij,kl) ;


   return(CalcInfo.twoel_ints[ijkl]) ;
}

/*
** tf_onel_ints(): Function lumps together one-electron contributions
**    so that h'_{ij} = h_{ij} - 1/2 SUM_k (ik|kj)
**    The term h' arises in the calculation of sigma1 and sigma2 via
**    equation (20) of Olsen, Roos, et. al. JCP 1988
**
*/
void tf_onel_ints(int printflag, std::string out)
{
   int i, j, k, ij, ik, kj, ikkj ;
   int nbf ;
   double *tei, *teptr ;
   double tval ;
   int ntri;

   /* set up some shorthand notation (speed up access) */
   nbf = CalcInfo.num_ci_orbs ;
   tei = CalcInfo.twoel_ints ;
   ntri = (nbf * (nbf + 1)) / 2;

   /* ok, new special thing for CASSCF...if there are *no* excitations
      into restricted orbitals, and if Parameters.fci=TRUE, then we
      do *not* want to sum over the restricted virts in h' or else
      we would need to account for RAS-out-of-space contributions
      (requiring fci=false).
    */
   if (Parameters.fci && (nbf > Parameters.ras3_lvl) &&
       Parameters.ras34_max == 0)
      nbf = Parameters.ras3_lvl;

   /* allocate space for the new array */
   CalcInfo.tf_onel_ints = init_array(ntri) ;

   /* fill up the new array */
   for (i=0,ij=0; i<nbf; i++)
      for (j=0; j<=i; j++) {
         tval = CalcInfo.onel_ints[ij] ;

         for (k=0; k<nbf; k++) {
            ik = ioff[MAX0(i,k)] + MIN0(i,k) ;
            kj = ioff[MAX0(k,j)] + MIN0(k,j) ;
            ikkj = ioff[ik] + kj ;
            teptr = tei + ikkj ;
            tval -= 0.5 * (*teptr) ;
            }

         CalcInfo.tf_onel_ints[ij++] = tval ;
         }
   /* print if necessary */
   if (printflag) {
      outfile->Printf( "\nh' matrix\n") ;
      print_array(CalcInfo.tf_onel_ints, nbf, "outfile") ;
      outfile->Printf( "\n") ;
      }
}



/*
** form_gmat(): Form the g matrix necessary for restriction to the RAS
**    subspaces (i.e. to eliminate contributions of out-of-space terms).
**    See equations (28-29) in Olsen, Roos, et. al. JCP 1988
**
*/
void form_gmat(int printflag, std::string out)
{
   int nbf ;
   double *tei, *oei ;
   double tval ;
   int i, j, k, ij, ii, ik, kj, ikkj, iiij ;


   /* set up some shorthand notation (speed up access) */
   nbf = CalcInfo.num_ci_orbs ;
   oei = CalcInfo.onel_ints ;
   tei = CalcInfo.twoel_ints ;

   /* allocate space for the new array */
   /* CalcInfo.gmat = init_matrix(nbf, nbf) ; */
   /* why not use init_blockmatix here? */
   CalcInfo.gmat = (double **) malloc (nbf * sizeof(double *));
   CalcInfo.gmat[0] = (double *) malloc (nbf * nbf * sizeof(double));
   for (i=1; i<nbf; i++) {
      CalcInfo.gmat[i] = CalcInfo.gmat[i-1] + nbf;
      }

   /* fill up the new array */
   for (i=0; i<nbf; i++) {
      for (j=i+1; j<nbf; j++) {
         ij = ioff[j] + i ;
         tval = oei[ij] ;
         for (k=0; k<i; k++) {
            ik = ioff[i] + k ;
            kj = ioff[j] + k ;
            ikkj = ioff[kj] + ik ;
            tval -= tei[ikkj] ;
            }
         CalcInfo.gmat[i][j] = tval ;
         }
      }

   for (i=0, ij=0; i<nbf; i++) {
      for (j=0; j<=i; j++,ij++) {
         tval = oei[ij] ;
         for (k=0; k<i; k++) {
            ik = ioff[i] + k ;
            kj = ioff[MAX0(k,j)] + MIN0(k,j) ;
            ikkj = ioff[ik] + kj ;
            tval -= tei[ikkj] ;
            }
         ii = ioff[i] + i ;
         iiij = ioff[ii] + ij ;
         if (i==j) tval -= 0.5 * tei[iiij] ;
         else tval -= tei[iiij] ;
         CalcInfo.gmat[i][j] = tval ;
         }
      }

   if (printflag) {
      outfile->Printf( "\ng matrix\n") ;
      print_mat(CalcInfo.gmat, nbf, nbf, "outfile") ;
      outfile->Printf( "\n") ;
      }
}


}} // namespace psi::detci

