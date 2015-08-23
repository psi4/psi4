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
#include <libmints/mints.h>
#include <psi4-dec.h>

#include <libtrans/integraltransform.h>
#include <libthce/thce.h>
#include <libthce/lreri.h>
#include <libfock/jk.h>

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
  SharedMatrix Cdrc = get_orbitals("DRC");
  SharedMatrix Cact = get_orbitals("ACT");
  SharedMatrix Cvir = get_orbitals("VIR");
  SharedMatrix Cfzv = get_orbitals("FZV");

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

  IntegralTransform *ints = new IntegralTransform(Cdrc, Cact, Cvir, Cfzv, spaces,
                            IntegralTransform::Restricted,
                            IntegralTransform::IWLOnly,
                            IntegralTransform::PitzerOrder,
                            IntegralTransform::OccAndVir,
                            true);

  // Incase we do two ci runs
  ints->set_keep_iwl_so_ints(true);
  ints->transform_tei(custom_space, custom_space, custom_space, custom_space);

  // Build temporary desired arrays
  int nmotri_full = (CalcInfo.nmo * (CalcInfo.nmo + 1)) / 2 ;
  double*  tmp_onel_ints1 = (double *) init_array(nmotri_full);
  double*  tmp_onel_ints2 = (double *) init_array(nmotri_full);

  // Read one electron integrals
  iwl_rdone(PSIF_OEI, PSIF_MO_FZC, tmp_onel_ints1, nmotri_full,
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

  // Read drc energy
  CalcInfo.edrc = ints->get_frozen_core_energy();

  // Read two electron integrals
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

  tf_onel_ints();
  form_gmat();
  delete ints;
}
void CIWavefunction::setup_dfmcscf_ints(){

  /// Grab and build basis sets
  boost::shared_ptr<BasisSet> primary = BasisSet::pyconstruct_orbital(
    Process::environment.molecule(), "BASIS", options_.get_str("BASIS"));
  boost::shared_ptr<BasisSet> auxiliary = BasisSet::pyconstruct_auxiliary(primary->molecule(),
      "DF_BASIS_SCF", options_.get_str("DF_BASIS_SCF"), "JKFIT",
      options_.get_str("BASIS"), primary->has_puream());

  /// Build JK object
  DFJK* jk = new DFJK(primary,auxiliary);

  if (options_["INTS_TOLERANCE"].has_changed())
      jk->set_cutoff(options_.get_double("INTS_TOLERANCE"));
  if (options_["PRINT"].has_changed())
      jk->set_print(options_.get_int("PRINT"));
  if (options_["DEBUG"].has_changed())
      jk->set_debug(options_.get_int("DEBUG"));
  if (options_["BENCH"].has_changed())
      jk->set_bench(options_.get_int("BENCH"));
  if (options_["DF_INTS_IO"].has_changed())
      jk->set_df_ints_io(options_.get_str("DF_INTS_IO"));
  if (options_["DF_FITTING_CONDITION"].has_changed())
      jk->set_condition(options_.get_double("DF_FITTING_CONDITION"));
  if (options_["DF_INTS_NUM_THREADS"].has_changed())
      jk->set_df_ints_num_threads(options_.get_int("DF_INTS_NUM_THREADS"));

  jk_ = boost::shared_ptr<JK>(jk);

  jk_->set_do_J(true);
  jk_->set_do_K(true);
  jk_->initialize();

  /// Build DF object
  dferi_ = DFERI::build(primary,auxiliary,options_);

  ints_init_ = true;
}
void CIWavefunction::transform_dfmcscf_ints(){

  if (!ints_init_) setup_dfmcscf_ints();

  // => AO C matrices <= //
  // We want pitzer order C matrix with appended Cact

  SharedMatrix Cocc = get_orbitals("DOCC");
  SharedMatrix Cact = get_orbitals("ACT");
  SharedMatrix Cvir = get_orbitals("VIR");

  int nao = AO2SO_->rowspi()[0];
  int nrot = Cocc->ncol() + Cact->ncol() + Cvir->ncol();
  int aoc_rowdim =  nrot + Cact->ncol();
  SharedMatrix AO_C = SharedMatrix(new Matrix("AO_C", nao, aoc_rowdim));

  double** AO_Cp = AO_C->pointer();
  for (int h=0, offset=0, offset_act=0; h < nirrep_; h++){
      int hnso = AO2SO_->colspi()[h];
      if (hnso == 0) continue;
      double** Up = AO2SO_->pointer(h);

      int noccpih = Cocc->colspi()[h];
      int nactpih = Cact->colspi()[h];
      int nvirpih = Cvir->colspi()[h];
      // occupied
      if (noccpih){
          double** CSOp = Cocc->pointer(h);
          C_DGEMM('N','N',nao,noccpih,hnso,1.0,Up[0],hnso,CSOp[0],noccpih,0.0,&AO_Cp[0][offset],aoc_rowdim);
          offset += noccpih;
      }
      // active
      if (nactpih){
          double** CSOp = Cact->pointer(h);
          C_DGEMM('N','N',nao,nactpih,hnso,1.0,Up[0],hnso,CSOp[0],nactpih,0.0,&AO_Cp[0][offset],aoc_rowdim);
          offset += nactpih;

          C_DGEMM('N','N',nao,nactpih,hnso,1.0,Up[0],hnso,CSOp[0],nactpih,0.0,&AO_Cp[0][offset_act + nrot],aoc_rowdim);
          offset_act += nactpih;
      }
      // virtual
      if (nvirpih){
          double** CSOp = Cvir->pointer(h);
          C_DGEMM('N','N',nao,nvirpih,hnso,1.0,Up[0],hnso,CSOp[0],nvirpih,0.0,&AO_Cp[0][offset],aoc_rowdim);
          offset += nvirpih;
      }
  }

  // => Compute DF ints <= //
  dferi_->clear();
  dferi_->set_C(AO_C);
  dferi_->add_space("N", 0, nrot);
  dferi_->add_space("a", nrot, aoc_rowdim);

  // Is it smart enough to order then untranspose?
  // In the future build this once then slice it
  dferi_->add_pair_space("aaQ", "a", "a");
  dferi_->add_pair_space("NaQ", "N", "a");
  dferi_->add_pair_space("NNQ", "N", "N");

  // dferi_->print_header(2);
  dferi_->compute();

  // => Compute onel ints <= //
  SharedMatrix Cdrc = get_orbitals("DRC");
  std::vector<SharedMatrix>& Cl = jk_->C_left();
  std::vector<SharedMatrix>& Cr = jk_->C_right();
  Cl.clear();
  Cr.clear();
  Cl.push_back(Cdrc);
  jk_->compute();
  Cl.clear();

  const std::vector<SharedMatrix>& J = jk_->J();
  const std::vector<SharedMatrix>& K = jk_->K();

  J[0]->scale(2.0);
  J[0]->subtract(K[0]);

  J[0]->add(H_);
  SharedMatrix onel_ints = Matrix::triplet(Cact, J[0], Cact, true, false, false);

  // Set onel ints
  for (int h=0, target=0, offset=0; h<nirrep_; h++){
    int nactpih = Cact->colspi()[h];
    if (!nactpih) continue;

    double** onep = onel_ints->pointer(h);
    for (int i=0; i<nactpih; i++){
      target += offset;
      for (int j=0; j<=i; j++){
        CalcInfo.onel_ints[target++] = onep[i][j];
      }
    }
    offset += nactpih;
  }

  // Compute Dropped core energy
  J[0]->add(H_);

  SharedMatrix D = Matrix::doublet(Cdrc, Cdrc, false, true);
  CalcInfo.edrc = J[0]->vector_dot(D);
  Cdrc.reset();
  D.reset();


  // Build two electron integrals
  int nQ = dferi_->size_Q();
  int nact = nact = Cact->ncol();

  std::map<std::string, boost::shared_ptr<Tensor> >& dfints = dferi_->ints();
  boost::shared_ptr<Tensor> aaQT = dfints["aaQ"];
  SharedMatrix aaQ(new Matrix("aaQ", nact * nact, nQ));

  double* aaQp = aaQ->pointer()[0];
  FILE* aaQF = aaQT->file_pointer();
  fseek(aaQF,0L,SEEK_SET);
  fread(aaQp, sizeof(double), nact * nact * nQ, aaQF);
  SharedMatrix actMO = Matrix::doublet(aaQ, aaQ, false, true);
  aaQ.reset();

  // Set twoel ints, ohboy
  int* myorder = new int[nact];
  for (int h=0, target=0, pos=0; h<nirrep_; h++){
    target += CalcInfo.dropped_docc[h];
    for (int i=0; i<CalcInfo.ci_orbs[h]; i++){
      myorder[pos++] = CalcInfo.reorder[target++] - CalcInfo.num_drc_orbs;
    }
    target += CalcInfo.dropped_uocc[h];
  }

  double** actMOp = actMO->pointer();
  int irel, jrel, krel, lrel, lmax;
  int target = 0;
  int ndrc = CalcInfo.num_drc_orbs;
  for (int i=0; i<nact; i++){
    irel = myorder[i];
    for (int j=0; j<=i; j++){
      jrel = myorder[j];
      for (int k=0; k<=i; k++){
        krel = myorder[k];
        lmax = (i==k) ? j+1 : k+1;
        for (int l=0; l<lmax; l++){
          lrel = myorder[l];
          // outfile->Printf("%d %d %d %d | %d %lf\n", i, j, k, l, target, val);
          CalcInfo.twoel_ints[target++] = actMOp[irel * nact + jrel][krel * nact + lrel];
  }}}}
  actMO.reset();
  delete[] myorder;

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
  tf_onel_ints();
  form_gmat();

}
void CIWavefunction::setup_mcscf_ints(){
  // We need to do a few weird things to make IntegralTransform work for us

  // Grab orbitals
  SharedMatrix Cdrc = get_orbitals("DRC");
  SharedMatrix Cact = get_orbitals("ACT");
  SharedMatrix Cvir = get_orbitals("VIR");
  SharedMatrix Cfzv = get_orbitals("FZV");

  // Need active and rot spaces
  std::vector<boost::shared_ptr<MOSpace> > spaces;

  std::vector<int> rot_orbitals(CalcInfo.num_ci_orbs, 0);
  std::vector<int> act_orbitals(CalcInfo.num_ci_orbs, 0);
  int act_orbnum = 0;
  int rot_orbnum = 0;
  for (int h = 0, rn = 0, an = 0; h < CalcInfo.nirreps; h++){
    act_orbnum += CalcInfo.dropped_docc[h];
    rot_orbnum += CalcInfo.frozen_docc[h];

    // Act space
    for (int i = 0; i < CalcInfo.ci_orbs[h]; i++){
      act_orbitals[an++] = act_orbnum;
    }
    act_orbnum += CalcInfo.dropped_uocc[h];

    int nrotorbs = CalcInfo.rstr_docc[h] + CalcInfo.ci_orbs[h] + CalcInfo.rstr_uocc[h];
    for (int i = 0; i < nrotorbs; i++){
      rot_orbitals[rn++] = rot_orbnum;
    }
    rot_orbnum += CalcInfo.frozen_uocc[h];
  }

  MOSpace* rot_space = new MOSpace('r', rot_orbitals);
  MOSpace* act_space = new MOSpace('a', act_orbitals);

  rot_space_ = boost::shared_ptr<MOSpace>(rot_space);
  act_space_ = boost::shared_ptr<MOSpace>(act_space);
  spaces.push_back(rot_space_);
  spaces.push_back(act_space_);


  // Now the occ space is active, the vir space is our rot space (FZC to FZV)
  IntegralTransform *ints = new IntegralTransform(Cdrc, Cact, Cvir, Cfzv, spaces,
                                                IntegralTransform::Restricted,
                                                IntegralTransform::DPDOnly,
                                                IntegralTransform::PitzerOrder,
                                                IntegralTransform::OccAndVir,
                                                true);
  ints_ = boost::shared_ptr<IntegralTransform>(ints);
  ints_->set_keep_iwl_so_ints(true);

  ints_init_ = true;


}
void CIWavefunction::transform_mcscf_ints(){

  if (!ints_init_) setup_mcscf_ints();

  // The orbital matrix need to be identical to the previous one
  ints_->set_orbitals(get_orbitals("ALL"));

  // We need (aa|aa), (aa|aN), (aa|NN), (aN|aN)
  ints_->transform_tei(act_space_, rot_space_, act_space_, rot_space_);

  // Half trans then work from there
  ints_->transform_tei_first_half(act_space_, act_space_);
  ints_->transform_tei_second_half(act_space_, act_space_, act_space_, rot_space_);
  ints_->transform_tei_second_half(act_space_, act_space_, rot_space_, rot_space_);
  ints_->transform_tei_second_half(act_space_, act_space_, act_space_, act_space_);


  CalcInfo.edrc = ints_->get_frozen_core_energy();

  // Build useful information
  int nmotri_full = (CalcInfo.nmo * (CalcInfo.nmo + 1)) / 2 ;
  int nmotri = (CalcInfo.num_ci_orbs * (CalcInfo.num_ci_orbs + 1)) / 2 ;

  // Build temporary arrays
  double*  tmp_onel_ints1 = (double *) init_array(nmotri_full);
  double*  tmp_onel_ints2 = (double *) init_array(nmotri_full);

  // Read one electron integrals
  iwl_rdone(PSIF_OEI, PSIF_MO_FZC, tmp_onel_ints1, nmotri_full,
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

  // Read two electron integrals
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

}


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
void CIWavefunction::tf_onel_ints()
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

   // /* allocate space for the new array */
   // CalcInfo.tf_onel_ints = init_array(ntri) ;

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
   // if (printflag) {
   //    outfile->Printf( "\nh' matrix\n") ;
   //    print_array(CalcInfo.tf_onel_ints, nbf, "outfile") ;
   //    outfile->Printf( "\n") ;
   //    }
}



/*
** form_gmat(): Form the g matrix necessary for restriction to the RAS
**    subspaces (i.e. to eliminate contributions of out-of-space terms).
**    See equations (28-29) in Olsen, Roos, et. al. JCP 1988
**
*/
void CIWavefunction::form_gmat()
{
   int nbf ;
   double *tei, *oei ;
   double tval ;
   int i, j, k, ij, ii, ik, kj, ikkj, iiij ;


   /* set up some shorthand notation (speed up access) */
   nbf = CalcInfo.num_ci_orbs ;
   oei = CalcInfo.onel_ints ;
   tei = CalcInfo.twoel_ints ;

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

   // if (printflag) {
   //    outfile->Printf( "\ng matrix\n") ;
   //    print_mat(CalcInfo.gmat, nbf, nbf, "outfile") ;
   //    outfile->Printf( "\n") ;
   //    }
}


}} // namespace psi::detci

