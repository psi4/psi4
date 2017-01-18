/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
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
 * @END LICENSE
 */

#include <map>
#include "dcft.h"
#include <cmath>
#include "psi4/libiwl/iwl.hpp"
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libtrans/integraltransform.h"
#include "defines.h"


using namespace std;

namespace psi{ namespace dcft{

/**
 * Computes the initial Hartree-Fock orbitals by either reading them from the
 * checkpoint file or computing a core Hamiltonian guess
 * for RHF reference.
 */
void
DCFTSolver::scf_guess_RHF()
{
    dcft_timer_on("DCFTSolver::rhf_guess");
    SharedMatrix T = SharedMatrix(new Matrix("SO basis kinetic energy integrals", nirrep_, nsopi_, nsopi_));
    SharedMatrix V = SharedMatrix(new Matrix("SO basis potential energy integrals", nirrep_, nsopi_, nsopi_));
    double *ints = init_array(ntriso_);

    IWL::read_one(psio_.get(), PSIF_OEI, PSIF_SO_T, ints, ntriso_, 0, 0, "outfile");
    T->set(ints);
    IWL::read_one(psio_.get(), PSIF_OEI, PSIF_SO_V, ints, ntriso_, 0, 0, "outfile");
    V->set(ints);
    free(ints);

    so_h_->add(T);
    so_h_->add(V);

    std::string guess = options_.get_str("DCFT_GUESS"); // The default DCFT_GUESS is mp2

    epsilon_a_->copy(reference_wavefunction_->epsilon_a().get());
    epsilon_b_->copy(epsilon_a_.get());
    Ca_->copy(reference_wavefunction_->Ca());
    Cb_->copy(Ca_);
    moFa_->copy(reference_wavefunction_->Fa());
    moFa_->transform(Ca_);
    moFb_->copy(moFa_);
    update_scf_density_RHF();

    dcft_timer_off("DCFTSolver::rhf_guess");

}

/**
* Uses the MO coefficients to form the SCF density matrices in the SO basis.
* @param Whether to damp the update or not
* @return RMS density change
*/
double
DCFTSolver::update_scf_density_RHF(bool damp)
{
    dcft_timer_on("DCFTSolver::update_rhf_density");

    double dampingFactor = options_.get_double("DAMPING_PERCENTAGE"); // The default DAMPING_PERCENTAGE is 0.0
    double newFraction = damp ? 1.0 : 1.0 - dampingFactor/100.0;
    size_t nElements = 0;
    double sumOfSquares = 0.0;
    Matrix old(kappa_so_a_); // Zero matrix
    for (int h = 0; h < nirrep_; ++h) {
        for (int mu = 0; mu < nsopi_[h]; ++mu) {
            for (int nu = 0; nu < nsopi_[h]; ++nu) {
                double val = 0.0;
                for (int i = 0; i < naoccpi_[h]; ++i)
                    val += Ca_->get(h, mu, i) * Ca_->get(h, nu, i);
                kappa_so_a_->set(h, mu, nu, newFraction*val + (1.0-newFraction) * kappa_so_a_->get(h, mu, nu));
                ++nElements;
                sumOfSquares += pow(val - old.get(h, mu, nu), 2.0);
            }
        }
    }

    kappa_so_b_->copy(kappa_so_a_);

    // We're not converged until the RMS error vector *and* the RMS density
    // changes are below the threshold
    dcft_timer_off("DCFTSolver::update_rhf_density");

    return sqrt(sumOfSquares / nElements);
}

/**
* Builds the G matrix for Fock matrix, the external potential (tau contracted with the integrals)
* and other tensors, if requested, out-of-core using the SO integrals.
* Also builds the AO basis tensors for the <VV||VV>, <vv||vv>, and <Vv|Vv> terms in the G and X intermediates
* All quantities are built simultaneously to reduce I/O.
*/
void
DCFTSolver::process_so_ints_RHF()
{
        dcft_timer_on("DCFTSolver::process_so_ints");

        IWL *iwl = new IWL(psio_.get(), PSIF_SO_TEI, int_tolerance_, 1, 1);

        Label *lblptr = iwl->labels();
        Value *valptr = iwl->values();

        double *Da = init_array(ntriso_);
        double *Ta = init_array(ntriso_);
        double *Ga = init_array(ntriso_);
        double *Va = init_array(ntriso_);

        int soOffset = 0;
        for(int h = 0; h < nirrep_; ++h){
            for(int mu = 0; mu < nsopi_[h]; ++ mu){
                for(int nu = 0; nu <= mu; ++ nu){
                    int muNu = INDEX((nu+soOffset), (mu+soOffset));
                    Da[muNu] = kappa_so_a_->get(h, mu, nu);
                    Ta[muNu] = tau_so_a_->get(h,mu,nu);
                }
            }
            soOffset += nsopi_[h];
        }

      double value;
      int Gc, Gd;
      int pqArr, qpArr, rsArr, srArr, qrArr, rqArr;
      int qsArr, sqArr, psArr, spArr, prArr, rpArr;
      int offset, labelIndex, p, q, r, s, h, counter;
      int **pq_row_start, **CD_row_start, **Cd_row_start, **cd_row_start;
      dpdbuf4 tau_temp, lambda;
      dpdbuf4 tau1_AO_ab, tau2_AO_ab;

      bool buildTensors = (options_.get_str("AO_BASIS") == "DISK");

      if(buildTensors){

          counter = 0;

          //Build the offset arrays needed for the DGEMM in half_transform
          pq_row_start = init_int_matrix(nirrep_, nirrep_);
          Cd_row_start = init_int_matrix(nirrep_, nirrep_);
          for(h = 0; h < nirrep_; ++h){
              for(Gc = 0, offset = 0; Gc < nirrep_; ++Gc){
                  Gd = Gc ^ h;
                  pq_row_start[h][Gc] = offset;
                  offset += nsopi_[Gc] * nsopi_[Gd];
              }
              for(Gc = 0, offset = 0; Gc < nirrep_; ++Gc){
                  Gd = Gc ^ h;
                  Cd_row_start[h][Gc] = offset;
                  offset += navirpi_[Gc] * nbvirpi_[Gd];
              }
          }

          dpd_set_default(_ints->get_dpd_id());

          /********** AB ***********/
          global_dpd_->buf4_init(&lambda, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
              ID("[O,O]"), ID("[V,V]"), 0, "Lambda SF <OO|VV>"); // Lambda <Oo|Vv>
          global_dpd_->buf4_init(&tau1_AO_ab, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[n,n]"),
              ID("[O,O]"), ID("[n,n]"), 0, "tau1AO <OO|nn>"); // tau1AO <Oo|nn>
          global_dpd_->buf4_scm(&tau1_AO_ab, 0.0);
          half_transform(&tau1_AO_ab, &lambda, avir_c_, avir_c_, navirpi_, navirpi_,
                  pq_row_start, Cd_row_start, true, 1.0, 0.0);
          global_dpd_->buf4_close(&lambda);
          global_dpd_->buf4_close(&tau1_AO_ab);

          // Now sort for better memory access patterns
          global_dpd_->buf4_init(&tau1_AO_ab, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[n,n]"),
              ID("[O,O]"), ID("[n,n]"), 0, "tau1AO <OO|nn>"); // tau1AO <Oo|nn>
          global_dpd_->buf4_sort(&tau1_AO_ab, PSIF_DCFT_DPD, rspq, ID("[n,n]"), ID("[O,O]"), "tau1AO <nn|OO>"); // tau1AO <nn|Oo>
          global_dpd_->buf4_close(&tau1_AO_ab);

          // Reopen the two AO dpd_buf4's
          global_dpd_->buf4_init(&tau1_AO_ab, PSIF_DCFT_DPD, 0, ID("[n,n]"), ID("[O,O]"),
              ID("[n,n]"), ID("[O,O]"), 0, "tau1AO <nn|OO>"); // tau1AO <nn|Oo>
          global_dpd_->buf4_init(&tau2_AO_ab, PSIF_DCFT_DPD, 0, ID("[n,n]"), ID("[O,O]"),
              ID("[n,n]"), ID("[O,O]"), 0, "tau2AO <nn|OO>"); // tau2AO <nn|Oo>
          global_dpd_->buf4_scm(&tau2_AO_ab, 0.0);

          // Now put stuff in memory
          for(int h = 0; h < nirrep_; ++h){
              global_dpd_->buf4_mat_irrep_init(&tau1_AO_ab, h);
              global_dpd_->buf4_mat_irrep_rd(&tau1_AO_ab, h);
              global_dpd_->buf4_mat_irrep_init(&tau2_AO_ab, h);

          }

      }

      bool lastBuffer;
      do{
          lastBuffer = iwl->last_buffer();
          for(int index = 0; index < iwl->buffer_count(); ++index){
              labelIndex = 4*index;
              p = abs((int) lblptr[labelIndex++]);
              q = (int) lblptr[labelIndex++];
              r = (int) lblptr[labelIndex++];
              s = (int) lblptr[labelIndex++];
              value = (double) valptr[index];
              if(buildTensors){
                  AO_contribute(&tau1_AO_ab, &tau2_AO_ab, p, q, r, s, value);
                  ++counter;
              }

              qpArr = pqArr = INDEX(p, q);
              srArr = rsArr = INDEX(r, s);
              prArr = rpArr = INDEX(p, r);
              qsArr = sqArr = INDEX(q, s);
              spArr = psArr = INDEX(p, s);
              qrArr = rqArr = INDEX(q, r);

              /* (pq|rs) */
              Ga[rsArr] += (Da[pqArr] + Da[pqArr]) * value;
              Va[rsArr] += (Ta[pqArr] + Ta[pqArr]) * value;
              if(q >= r){
                  Ga[qrArr] -= Da[psArr] * value;
                  Va[qrArr] -= Ta[psArr] * value;
              }

              if(p!=q && r!=s && pqArr!=rsArr){
                  /* (pq|sr) */
                  if(s >= r){
                      Ga[srArr] += (Da[pqArr] + Da[pqArr]) * value;
                      Va[srArr] += (Ta[pqArr] + Ta[pqArr]) * value;
                  }
                  if(q >= s){
                      Ga[qsArr] -= Da[prArr] * value;
                      Va[qsArr] -= Ta[prArr] * value;
                  }

                  /* (qp|rs) */
                  if(r >= s){
                      Ga[rsArr] += (Da[qpArr] + Da[qpArr]) * value;
                      Va[rsArr] += (Ta[qpArr] + Ta[qpArr]) * value;
                  }
                  if(p >= r){
                      Ga[prArr] -= Da[qsArr] * value;
                      Va[prArr] -= Ta[qsArr] * value;
                  }

                  /* (qp|sr) */
                  if(s >= r){
                      Ga[srArr] += (Da[qpArr] + Da[qpArr]) * value;
                      Va[srArr] += (Ta[qpArr] + Ta[qpArr]) * value;
                  }
                  if(p >= s){
                      Ga[psArr] -= Da[qrArr] * value;
                      Va[psArr] -= Ta[qrArr] * value;
                  }

                  /* (rs|pq) */
                  if(p >= q){
                      Ga[pqArr] += (Da[rsArr] + Da[rsArr]) * value;
                      Va[pqArr] += (Ta[rsArr] + Ta[rsArr]) * value;
                  }
                  if(s >= p){
                      Ga[spArr] -= Da[rqArr] * value;
                      Va[spArr] -= Ta[rqArr] * value;
                  }

                  /* (sr|pq) */
                  if(p >= q){
                      Ga[pqArr] += (Da[srArr] + Da[srArr]) * value;
                      Va[pqArr] += (Ta[srArr] + Ta[srArr]) * value;
                  }
                  if(r >= p){
                      Ga[rpArr] -= Da[sqArr] * value;
                      Va[rpArr] -= Ta[sqArr] * value;
                  }

                  /* (rs|qp) */
                  if(q >= p){
                      Ga[qpArr] += (Da[rsArr] + Da[rsArr]) * value;
                      Va[qpArr] += (Ta[rsArr] + Ta[rsArr]) * value;
                  }
                  if(s >= q){
                      Ga[sqArr] -= Da[rpArr] * value;
                      Va[sqArr] -= Ta[rpArr] * value;
                  }

                  /* (sr|qp) */
                  if(q >= p){
                      Ga[qpArr] += (Da[srArr] + Da[srArr]) * value;
                      Va[qpArr] += (Ta[srArr] + Ta[srArr]) * value;
                  }
                  if(r >= q){
                      Ga[rqArr] -= Da[spArr] * value;
                      Va[rqArr] -= Ta[spArr] * value;
                  }
              }else if(p!=q && r!=s && pqArr==rsArr){
                  /* (pq|sr) */
                  if(s >= r){
                      Ga[srArr] += (Da[pqArr] + Da[pqArr]) * value;
                      Va[srArr] += (Ta[pqArr] + Ta[pqArr]) * value;
                  }
                  if(q >= s){
                      Ga[qsArr] -= Da[prArr] * value;
                      Va[qsArr] -= Ta[prArr] * value;
                  }
                  /* (qp|rs) */
                  if(r >= s){
                      Ga[rsArr] += (Da[qpArr] + Da[qpArr]) * value;
                      Va[rsArr] += (Ta[qpArr] + Ta[qpArr]) * value;
                  }
                  if(p >= r){
                      Ga[prArr] -= Da[qsArr] * value;
                      Va[prArr] -= Ta[qsArr] * value;
                  }

                  /* (qp|sr) */
                  if(s >= r){
                      Ga[srArr] += (Da[qpArr] + Da[qpArr]) * value;
                      Va[srArr] += (Ta[qpArr] + Ta[qpArr]) * value;
                  }
                  if(p >= s){
                      Ga[psArr] -= Da[qrArr] * value;
                      Va[psArr] -= Ta[qrArr] * value;
                  }
              }else if(p!=q && r==s){
                  /* (qp|rs) */
                  if(r >= s){
                      Ga[rsArr] += (Da[qpArr] + Da[qpArr]) * value;
                      Va[rsArr] += (Ta[qpArr] + Ta[qpArr]) * value;
                  }
                  if(p >= r){
                      Ga[prArr] -= Da[qsArr] * value;
                      Va[prArr] -= Ta[qsArr] * value;
                  }

                  /* (rs|pq) */
                  if(p >= q){
                      Ga[pqArr] += (Da[rsArr] + Da[rsArr]) * value;
                      Va[pqArr] += (Ta[rsArr] + Ta[rsArr]) * value;
                  }
                  if(s >= p){
                      Ga[spArr] -= Da[rqArr] * value;
                      Va[spArr] -= Ta[rqArr] * value;
                  }

                  /* (rs|qp) */
                  if(q >= p){
                      Ga[qpArr] += (Da[rsArr] + Da[rsArr]) * value;
                      Va[qpArr] += (Ta[rsArr] + Ta[rsArr]) * value;
                  }
                  if(s >= q){
                      Ga[sqArr] -= Da[rpArr] * value;
                      Va[sqArr] -= Ta[rpArr] * value;
                  }
              }else if(p==q && r!=s){
                  /* (pq|sr) */
                  if(s >= r){
                      Ga[srArr] += (Da[pqArr] + Da[pqArr]) * value;
                      Va[srArr] += (Ta[pqArr] + Ta[pqArr]) * value;
                  }
                  if(q >= s){
                      Ga[qsArr] -= Da[prArr] * value;
                      Va[qsArr] -= Ta[prArr] * value;
                  }

                  /* (rs|pq) */
                  if(p >= q){
                      Ga[pqArr] += (Da[rsArr] + Da[rsArr]) * value;
                      Va[pqArr] += (Ta[rsArr] + Ta[rsArr]) * value;
                  }
                  if(s >= p){
                      Ga[spArr] -= Da[rqArr] * value;
                      Va[spArr] -= Ta[rqArr] * value;
                  }

                  /* (sr|pq) */
                  if(p >= q){
                      Ga[pqArr] += (Da[srArr] + Da[srArr]) * value;
                      Va[pqArr] += (Ta[srArr] + Ta[srArr]) * value;
                  }
                  if(r >= p){
                      Ga[rpArr] -= Da[sqArr] * value;
                      Va[rpArr] -= Ta[sqArr] * value;
                  }
              }else if(p==q && r==s && pqArr!=rsArr){
                  /* (rs|pq) */
                  if(p >= q){
                      Ga[pqArr] += (Da[rsArr] + Da[rsArr]) * value;
                      Va[pqArr] += (Ta[rsArr] + Ta[rsArr]) * value;
                  }
                  if(s >= p){
                      Ga[spArr] -= Da[rqArr] * value;
                      Va[spArr] -= Ta[rqArr] * value;
                  }
              }
          } /* end loop through current buffer */
          if(!lastBuffer) iwl->fetch();
      }while(!lastBuffer);
      iwl->set_keep_flag(1);
      delete iwl;
      if(buildTensors){
          if(print_ > 1){
              outfile->Printf( "Processed %d SO integrals each for AB\n", counter);
          }
          for(int h = 0; h < nirrep_; ++h){
              global_dpd_->buf4_mat_irrep_wrt(&tau2_AO_ab, h);
              global_dpd_->buf4_mat_irrep_close(&tau1_AO_ab, h);
              global_dpd_->buf4_mat_irrep_close(&tau2_AO_ab, h);
          }

          global_dpd_->buf4_close(&tau1_AO_ab);
          global_dpd_->buf4_close(&tau2_AO_ab);

          /********** AB ***********/
          global_dpd_->buf4_init(&tau2_AO_ab, PSIF_DCFT_DPD, 0, ID("[n,n]"), ID("[O,O]"),
                  ID("[n,n]"), ID("[O,O]"), 0, "tau2AO <nn|OO>"); // tau2AO <nn|Oo>
          global_dpd_->buf4_sort(&tau2_AO_ab, PSIF_DCFT_DPD, rspq, ID("[O,O]"), ID("[n,n]"), "tau2AO <OO|nn>"); // tau2AO <Oo|nn>
          global_dpd_->buf4_close(&tau2_AO_ab);
          global_dpd_->buf4_init(&tau2_AO_ab, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[n,n]"),
                 ID("[O,O]"), ID("[n,n]"), 0, "tau2AO <OO|nn>"); // tau2AO <Oo|nn>
          global_dpd_->buf4_init(&tau_temp, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "tau(temp) SF <OO|VV>"); // tau(temp) <Oo|Vv>
          global_dpd_->buf4_scm(&tau_temp, 0.0);
          half_transform(&tau2_AO_ab, &tau_temp, avir_c_, avir_c_, navirpi_, navirpi_,
                  pq_row_start, Cd_row_start, false, 1.0, 0.0);
          global_dpd_->buf4_close(&tau2_AO_ab);
          global_dpd_->buf4_close(&tau_temp);

          free_int_matrix(pq_row_start);
          free_int_matrix(Cd_row_start);

      }

      // Build the Fock matrices from the H and G matrices
      soOffset = 0;
      for(int h = 0; h < nirrep_; ++h){
          for(int mu = 0; mu < nsopi_[h]; ++mu){
              for(int nu = 0; nu <= mu; ++nu){
                  int muNu = INDEX((nu+soOffset), (mu+soOffset));
                  double aVal   = Ga[muNu];
                  double aGTVal = Va[muNu];
                  Fa_->add(h, mu, nu, aVal);
                  g_tau_a_->set(h, mu, nu, aGTVal);
                  if(mu != nu){
                      Fa_->add(h, nu, mu, aVal);
                      g_tau_a_->set(h, nu, mu, aGTVal);
                  }
              }
          }
          soOffset += nsopi_[h];
      }

      Fb_->copy(Fa_);
      g_tau_b_->copy(g_tau_a_);

      delete [] Ta;
      delete [] Va;
      delete [] Da;
      delete [] Ga;

      dcft_timer_off("DCFTSolver::process_so_ints");
}


/**
* Computes the SCF energy from the latest Fock and density matrices.
*/
void
DCFTSolver::compute_scf_energy_RHF()
{
    dcft_timer_on("DCFTSolver::compute_scf_energy");

    // Escf = eNuc + 0.5 * (H + F) * (kappa + tau)
    scf_energy_ = enuc_;
    scf_energy_ += kappa_so_a_->vector_dot(so_h_);
    scf_energy_ += tau_so_a_->vector_dot(so_h_);

    if(options_.get_str("DCFT_TYPE") == "DF" && options_.get_str("AO_BASIS") == "NONE"){
        mo_gammaA_->add(kappa_mo_a_);

        scf_energy_ += mo_gammaA_->vector_dot(moFa_);
    }
    else{
        scf_energy_ += kappa_so_a_->vector_dot(Fa_);
        scf_energy_ += tau_so_a_->vector_dot(Fa_);
    }

    dcft_timer_off("DCFTSolver::compute_scf_energy");
}

double
DCFTSolver::compute_scf_error_vector_RHF()
{
    dcft_timer_on("DCFTSolver::compute_scf_error_vector");

    size_t nElements = 0;
    double sumOfSquares = 0.0;
    SharedMatrix tmp1(new Matrix("tmp1", nirrep_, nsopi_, nsopi_));
    SharedMatrix tmp2(new Matrix("tmp2", nirrep_, nsopi_, nsopi_));
    // form FDS
    tmp1->gemm(false, false, 1.0, kappa_so_a_, ao_s_, 0.0);
    scf_error_a_->gemm(false, false, 1.0, Fa_, tmp1, 0.0);
    // form SDF
    tmp1->gemm(false, false, 1.0, kappa_so_a_, Fa_, 0.0);
    tmp2->gemm(false, false, 1.0, ao_s_, tmp1, 0.0);
    scf_error_a_->subtract(tmp2);
    // Orthogonalize
    scf_error_a_->transform(s_half_inv_);
    scf_error_b_->copy(scf_error_a_);

    for(int h = 0; h < nirrep_; ++h){
        for(int p = 0; p < nsopi_[h]; ++p){
            for(int q = 0; q < nsopi_[h]; ++q){
                nElements += 2;
                sumOfSquares += pow(scf_error_a_->get(h, p, q), 2.0);
                sumOfSquares += pow(scf_error_b_->get(h, p, q), 2.0);
            }
        }
    }
    dcft_timer_off("DCFTSolver::compute_scf_error_vector");
    return sqrt(sumOfSquares / nElements);
}

}} // Namespace
