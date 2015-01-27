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

#ifndef _psi_src_bin_detci_mcscf_h_
#define _psi_src_bin_detci_mcscf_h_

#include <libmints/wavefunction.h>
#include "psi4-dec.h"
#include "libparallel/ParallelPrinter.h"
#include <libdiis/diismanager.h>
#include <libdiis/diisentry.h>
#include "MCSCF_indpairs.h"

namespace boost {
template<class T> class shared_ptr;
}

namespace psi {

namespace detci {

class MCSCF 
{
private:
    // Core functions
    Options& options_;

    // Independent pairs
    IndepPairs IndPairs;
    int num_indep_pairs_;

    // MCSCF.cc
    void calc_gradient(void);
    void bfgs_hessian(void);
    void ds_hessian(void);
    void calc_hessian(void);
    void scale_gradient(void);
    int check_conv(void);
    int take_step(void);
    void rotate_orbs(void);
    double** lagcalc(double **OPDM, double *TPDM, double *h, double *TwoElec,
                   int nmo, int npop, int print_lvl, int lag_file);
    void form_independent_pairs(void);

    // MCSCF_f_act
    void form_F_act(void);

    // get_mo_info
    void read_cur_orbs(void);

    //Hessians
    void form_appx_diag_mo_hess(int npairs, int *ppair, int *qpair,
                           double *F_core, double *tei, double **opdm,
                           double *tpdm, double *F_act, int firstact,
                           int lastact, double *hess);
    void form_diag_mo_hess(int npairs, int *ppair, int *qpair,
                           double *F_core, double *tei, double **opdm,
                           double *tpdm, double *F_act, int firstact,
                           int lastact, double *hess);
    void form_full_mo_hess(int npairs, int *ppair, int *qpair,
                      double *oei, double *tei, double **opdm,
                      double *tpdm, double **lag, double **hess);
    void form_diag_mo_hess_yy(int npairs, int *ppair, int *qpair,
                      double *oei, double *tei, double **opdm,
                      double *tpdm, double **lag, double *hess);

    // Rotate
    void calc_orb_step(int npairs, double *grad, double *hess_diag,
                       double *theta);
    void calc_orb_step_full(int npairs, double *grad, double **hess,
                       double *theta);
    void calc_orb_step_bfgs(int npairs, double *grad, double **hess,
                       double *theta);
    void print_step(int iter, int npairs, int steptype,
                    OutFile& IterSummaryOut);
    void postmult_by_U(int irrep, int dim, double **mo_coeffs,
                       int npairs, int *p_arr, int *q_arr,
                       double *theta_arr);
    void premult_by_U(int irrep, int dim, double **mo_coeffs,
                      int npairs, int *p_arr, int *q_arr,
                      double *theta_arr);
    void postmult_by_exp_R(int irrep, int dim, double **mat,
                           int npairs, int *p_arr, int *q_arr,
                           double *theta_arr);
    int  read_ref_orbs(void);
    int  write_ref_orbs(void);
    void read_thetas(int npairs);
    void write_thetas(int npairs);

    void calc_dE_dT(int n, double **dEU, int npairs, int *ppair,
                    int *qpair, double *theta, double *dET);

    // DIIS information
    boost::shared_ptr<DIISManager> diis_manager_;
    int diis_iter_;
    int ndiis_vec_;


public:

    /// Constructor
    MCSCF(Options& options);
    ~MCSCF();   
 
    void mcscf_cleanup(void);
    int mcscf_update(int iter, Options &options, OutFile& IterSummaryOut);
    void mcscf_title(void);
    void mcscf_get_mo_info(Options &options);


};
}} /* End Namespaces */

#endif // _psi_src_bin_detci_mcscf_h_
