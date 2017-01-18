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
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libtrans/integraltransform.h"
#include "defines.h"


using namespace std;

namespace psi{ namespace dcft{


  /**
  * Checks to make sure that the phase, and ordering, of the MOs is consistent
  * with the previous set.
  *
  * @return Whether the phase correction was successful
  */
  bool
  DCFTSolver::correct_mo_phases(bool dieOnError)
  {
      dcft_timer_on("DCFTSolver::correct_mo_phases()");

#if 1
      Matrix temp("temp", nirrep_, nsopi_, nsopi_);
      Matrix overlap("Old - New Overlap", nirrep_, nsopi_, nsopi_);

      temp.gemm(true, false, 1.0, old_ca_, ao_s_, 0.0);
      overlap.gemm(false, false, 1.0, temp, Ca_, 0.0);
      temp.copy(Ca_);
      std::map<int, int> mosUsed;
      int offset = 0;
      for(int h = 0; h < nirrep_; ++h){
          for(int oldMO = 0; oldMO < nsopi_[h]; ++oldMO){
              int bestMO = 0;
              double maximumProjection = 0.0;
              double prefactor = 0.0;
              for(int newMO = 0; newMO < nsopi_[h]; ++newMO){
                  double val = overlap.get(h, oldMO, newMO);
                  if(fabs(val) > maximumProjection){
                      maximumProjection = fabs(val);
                      bestMO = newMO;
                      prefactor = val < 0.0 ? -1.0 : 1.0;
                  }
              }
              // Now we've found the MO to use, check it's not been used already then
              // copy it over.
              if(mosUsed[bestMO + offset]++){
                  if(dieOnError){
                      overlap.print();
                      old_ca_->print();
                      temp.print();
                      throw SanityCheckError("Duplicate MOs used in phase check", __FILE__, __LINE__);
                  }else{
                      // Copy Ca back from temp
                      Ca_->copy(temp);
                      dcft_timer_off("DCFTSolver::correct_mo_phases()");
                      return false;
                  }
              }
              for(int so = 0; so < nsopi_[h]; ++so){
                  Ca_->set(h, so, oldMO, prefactor * temp.get(h, so, bestMO));
              }
          }
          offset += nsopi_[h];
      }

      temp.gemm(true, false, 1.0, old_cb_, ao_s_, 0.0);
      overlap.gemm(false, false, 1.0, temp, Cb_, 0.0);
      temp.copy(Cb_);
      mosUsed.clear();
      offset = 0;
      for(int h = 0; h < nirrep_; ++h){
          for(int oldMO = 0; oldMO < nsopi_[h]; ++oldMO){
              int bestMO = 0;
              double bestOverlap = 0.0;
              double prefactor = 0.0;
              for(int newMO = 0; newMO < nsopi_[h]; ++newMO){
                  double val = overlap.get(h, oldMO, newMO);
                  if(fabs(val) > bestOverlap){
                      bestOverlap = fabs(val);
                      bestMO = newMO;
                      prefactor = val < 0.0 ? -1.0 : 1.0;
                  }
              }
              // Now we've found the MO to use, check it's not been used already then
              // copy it over.
              if(mosUsed[bestMO + offset]++){
                  if(dieOnError){
                      overlap.print();
                      old_cb_->print();
                      temp.print();
                      throw SanityCheckError("Duplicate MOs used in phase check", __FILE__, __LINE__);
                  }else{
                      // Copy Cb back from temp
                      Cb_->copy(temp);
                      dcft_timer_off("DCFTSolver::correct_mo_phases()");
                      return false;
                  }
              }
              for(int so = 0; so < nsopi_[h]; ++so){
                  Cb_->set(h, so, oldMO, prefactor * temp.get(h, so, bestMO));
              }
          }
          offset += nsopi_[h];
      }
#else
      Matrix temp("temp", nirrep_, nsopi_, nsopi_);
      Matrix overlap("Old - New Overlap", nirrep_, nsopi_, nsopi_);

      temp.gemm(true, false, 1.0, old_ca_, ao_s_, 0.0);
      overlap.gemm(false, false, 1.0, temp, Ca_, 0.0);
      temp.copy(Ca_);
      int offset = 0;
      std::vector<std::pair<double, int> > proj_orb;
      for(int h = 0; h < nirrep_; ++h){
          proj_orb.clear();
          for(int mo = 0; mo < nsopi_[h]; ++mo){
              double proj = 0.0;
              for(int occ = 0; occ < nalphapi_[h]; ++occ){
                  proj += overlap.get(h, occ, mo);
              }
              proj_orb.push_back(std::make_pair(proj, mo));
          }
              // Now we've found the MO to use, check it's not been used already then
              // copy it over.
              if(mosUsed[bestMO + offset]++){
                  if(dieOnError){
                      overlap.print();
                      old_ca_->print();
                      temp.print();
                      throw SanityCheckError("Duplicate MOs used in phase check", __FILE__, __LINE__);
                  }else{
                      // Copy Ca back from temp
                      Ca_->copy(temp);
                      dcft_timer_off("DCFTSolver::correct_mo_phases()");
                      return false;
                  }
              }
              for(int so = 0; so < nsopi_[h]; ++so){
                  Ca_->set(h, so, oldMO, prefactor * temp.get(h, so, bestMO));
              }
          }
          offset += nsopi_[h];
      }

      temp.gemm(true, false, 1.0, old_cb_, ao_s_, 0.0);
      overlap.gemm(false, false, 1.0, temp, Cb_, 0.0);
      temp.copy(Cb_);
      mosUsed.clear();
      offset = 0;
      for(int h = 0; h < nirrep_; ++h){
          for(int oldMO = 0; oldMO < nsopi_[h]; ++oldMO){
              int bestMO = 0;
              double bestOverlap = 0.0;
              double prefactor = 0.0;
              for(int newMO = 0; newMO < nsopi_[h]; ++newMO){
                  double val = overlap.get(h, oldMO, newMO);
                  if(fabs(val) > bestOverlap){
                      bestOverlap = fabs(val);
                      bestMO = newMO;
                      prefactor = val < 0.0 ? -1.0 : 1.0;
                  }
              }
              // Now we've found the MO to use, check it's not been used already then
              // copy it over.
              if(mosUsed[bestMO + offset]++){
                  if(dieOnError){
                      overlap.print();
                      old_cb_->print();
                      temp.print();
                      throw SanityCheckError("Duplicate MOs used in phase check", __FILE__, __LINE__);
                  }else{
                      // Copy Cb back from temp
                      Cb_->copy(temp);
                      dcft_timer_off("DCFTSolver::correct_mo_phases()");
                      return false;
                  }
              }
              for(int so = 0; so < nsopi_[h]; ++so){
                  Cb_->set(h, so, oldMO, prefactor * temp.get(h, so, bestMO));
              }
          }
          offset += nsopi_[h];
      }
#endif
      dcft_timer_off("DCFTSolver::correct_mo_phases()");
      return true;
  }

  /**
  * Figures out the orbital symmetries and prints them, ordered by energy and type
  */
  void
  DCFTSolver::print_orbital_energies()
  {
      std::vector<std::pair<double, int> > aPairs;
      std::vector<std::pair<double, int> > bPairs;
      for (int h = 0; h < nirrep_; ++h) {
          for (int i=0; i < nsopi_[h]; ++i){
              aPairs.push_back(make_pair(epsilon_a_->get(h, i), h));
              bPairs.push_back(make_pair(epsilon_b_->get(h, i), h));
          }
      }
      sort(aPairs.begin(), aPairs.end());
      sort(bPairs.begin(), bPairs.end());

      int *aIrrepCount = init_int_array(nirrep_);
      int *bIrrepCount = init_int_array(nirrep_);
      char **irrepLabels = molecule_->irrep_labels();

      outfile->Printf( "\n\tOrbital energies (a.u.):\n\t\tAlpha occupied orbitals\n\t\t");
      for (int i = 0, count = 0; i < nalpha_; ++i, ++count) {
          int irrep = aPairs[i].second;
          outfile->Printf( "%4d%-4s%11.6f  ", ++aIrrepCount[irrep], irrepLabels[irrep], aPairs[i].first);
          if (count % 4 == 3 && i != nalpha_)
              outfile->Printf( "\n\t\t");
      }
      outfile->Printf( "\n\n\t\tBeta occupied orbitals\n\t\t");
      for (int i = 0, count = 0; i < nbeta_; ++i, ++count) {
          int irrep = bPairs[i].second;
          outfile->Printf( "%4d%-4s%11.6f  ", ++bIrrepCount[irrep], irrepLabels[irrep], bPairs[i].first);
          if (count % 4 == 3 && i != nbeta_)
              outfile->Printf( "\n\t\t");
      }
      outfile->Printf( "\n\n\t\tAlpha virtual orbitals\n\t\t");
      for (int i = nalpha_, count = 0; i < nmo_; ++i, ++count) {
          int irrep = aPairs[i].second;
          outfile->Printf( "%4d%-4s%11.6f  ", ++aIrrepCount[irrep], irrepLabels[irrep], aPairs[i].first);
          if (count % 4 == 3 && i != nmo_)
              outfile->Printf( "\n\t\t");
      }
      outfile->Printf( "\n\n\t\tBeta virtual orbitals\n\t\t");
      for (int i = nbeta_, count = 0; i < nmo_; ++i, ++count) {
          int irrep = bPairs[i].second;
          outfile->Printf( "%4d%-4s%11.6f  ", ++bIrrepCount[irrep], irrepLabels[irrep], bPairs[i].first);
          if (count % 4 == 3 && i != nmo_)
              outfile->Printf( "\n\t\t");
      }
      outfile->Printf( "\n\n");

      outfile->Printf( "\n\tIrrep              ");
      for(int h = 0; h < nirrep_; ++h) outfile->Printf( "%4s ", irrepLabels[h]);
      outfile->Printf( "\n\t-------------------");
      for(int h = 0; h < nirrep_; ++h) outfile->Printf( "-----");
      outfile->Printf( "\n\t#Symmetry Orbitals ");
      for(int h = 0; h < nirrep_; ++h) outfile->Printf( "%4d ", nsopi_[h]);
      outfile->Printf( "\n\t#Alpha Occupied    ");
      for(int h = 0; h < nirrep_; ++h) outfile->Printf( "%4d ", naoccpi_[h]);
      outfile->Printf( "\n\t#Beta Occupied     ");
      for(int h = 0; h < nirrep_; ++h) outfile->Printf( "%4d ", nboccpi_[h]);
      outfile->Printf( "\n\t#Alpha Virtual     ");
      for(int h = 0; h < nirrep_; ++h) outfile->Printf( "%4d ", navirpi_[h]);
      outfile->Printf( "\n\t#Beta Virtual      ");
      for(int h = 0; h < nirrep_; ++h) outfile->Printf( "%4d ", nbvirpi_[h]);
      outfile->Printf( "\n\t-------------------");
      for(int h = 0; h < nirrep_; ++h) outfile->Printf( "-----");
      outfile->Printf( "\n\n");

      if(print_ > 2){
          Ca_->print();
          Cb_->print();
      }
      for (int h = 0; h < nirrep_; ++h)
          delete [] irrepLabels[h];
      delete[] irrepLabels;
      delete[] aIrrepCount;
      delete[] bIrrepCount;
  }

  /**
   * Computes the initial Hartree-Fock orbitals by either reading them from the
   * checkpoint file or computing a core Hamiltonian guess.
   */
  void
  DCFTSolver::scf_guess()
  {
      dcft_timer_on("DCFTSolver::scf_guess");
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

      std::string guess = options_.get_str("DCFT_GUESS");
      epsilon_a_->copy(reference_wavefunction_->epsilon_a().get());
      epsilon_b_->copy(reference_wavefunction_->epsilon_b().get());
      Ca_->copy(reference_wavefunction_->Ca());
      Cb_->copy(reference_wavefunction_->Cb());
      moFa_->copy(reference_wavefunction_->Fa());
      moFa_->transform(Ca_);
      moFb_->copy(reference_wavefunction_->Fb());
      moFb_->transform(Cb_);
      update_scf_density();

      dcft_timer_off("DCFTSolver::scf_guess");
  }


  /**
  * Computes the SCF energy from the latest Fock and density matrices.
  */
  void
  DCFTSolver::compute_scf_energy()
  {
      dcft_timer_on("DCFTSolver::compute_scf_energy");

      // Escf = eNuc + 0.5 * (H + F) * (kappa + tau)
      scf_energy_ = enuc_;
      scf_energy_ += 0.5 * kappa_so_a_->vector_dot(so_h_);
      scf_energy_ += 0.5 * kappa_so_b_->vector_dot(so_h_);
      scf_energy_ += 0.5 * tau_so_a_->vector_dot(so_h_);
      scf_energy_ += 0.5 * tau_so_b_->vector_dot(so_h_);

      if (options_.get_str("DCFT_TYPE") == "DF" && options_.get_str("AO_BASIS") == "NONE"){
          mo_gammaA_->add(kappa_mo_a_);
          mo_gammaB_->add(kappa_mo_b_);

          scf_energy_ += 0.5 * mo_gammaA_->vector_dot(moFa_);
          scf_energy_ += 0.5 * mo_gammaB_->vector_dot(moFb_);
      }
      else{
          scf_energy_ += 0.5 * kappa_so_a_->vector_dot(Fa_);
          scf_energy_ += 0.5 * kappa_so_b_->vector_dot(Fb_);
          scf_energy_ += 0.5 * tau_so_a_->vector_dot(Fa_);
          scf_energy_ += 0.5 * tau_so_b_->vector_dot(Fb_);
      }

      dcft_timer_off("DCFTSolver::compute_scf_energy");
  }

  /**
  * Computes the SCF error vector by transforming the Fock matrices to the
  * MO basis and computing [F, Kappa], and the RMS value of this quantity.
  * @return RMS error
  */
  double
  DCFTSolver::compute_scf_error_vector()
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

      // form FDS
      tmp1->gemm(false, false, 1.0, kappa_so_b_, ao_s_, 0.0);
      scf_error_b_->gemm(false, false, 1.0, Fb_, tmp1, 0.0);
      // form SDF
      tmp1->gemm(false, false, 1.0, kappa_so_b_, Fb_, 0.0);
      tmp2->gemm(false, false, 1.0, ao_s_, tmp1, 0.0);
      scf_error_b_->subtract(tmp2);
      // Orthogonalize
      scf_error_b_->transform(s_half_inv_);

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


  /**
  * Uses the MO coefficients to form the SCF density matrices in the SO basis.
  * @param Whether to damp the update or not
  * @return RMS density change
  */
  double
  DCFTSolver::update_scf_density(bool damp)
  {
      dcft_timer_on("DCFTSolver::update_scf_density");

      double dampingFactor = options_.get_double("DAMPING_PERCENTAGE");
      double newFraction = damp ? 1.0 : 1.0 - dampingFactor/100.0;
      size_t nElements = 0;
      double sumOfSquares = 0.0;
      Matrix old(kappa_so_a_);
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
      old.copy(kappa_so_b_);
      for (int h = 0; h < nirrep_; ++h) {
          for (int mu = 0; mu < nsopi_[h]; ++mu) {
              for (int nu = 0; nu < nsopi_[h]; ++nu) {
                  double val = 0.0;
                  for (int i = 0; i < nboccpi_[h]; ++i)
                      val += Cb_->get(h, mu, i) * Cb_->get(h, nu, i);
                  kappa_so_b_->set(h, mu, nu, newFraction*val + (1.0-newFraction) * kappa_so_b_->get(h, mu, nu));
                  ++nElements;
                  sumOfSquares += pow(val - old.get(h, mu, nu), 2.0);
              }
          }
      }
      // We're not converged until the RMS error vector *and* the RMS density
      // changes are below the threshold
      dcft_timer_off("DCFTSolver::update_scf_density");

      return sqrt(sumOfSquares / nElements);
  }


  /**
  * Builds the G matrix for Fock matrix, the external potential (tau contracted with the integrals)
  * and other tensors, if requested, out-of-core using the SO integrals.
  * Also builds the AO basis tensors for the <VV||VV>, <vv||vv>, and <Vv|Vv> terms in the G and X intermediates
  * All quantities are built simultaneously to reduce I/O.
  */
  void
  DCFTSolver::process_so_ints()
  {
      dcft_timer_on("DCFTSolver::process_so_ints");

      IWL *iwl = new IWL(psio_.get(), PSIF_SO_TEI, int_tolerance_, 1, 1);

      Label *lblptr = iwl->labels();
      Value *valptr = iwl->values();

      double *Da = init_array(ntriso_);
      double *Db = init_array(ntriso_);
      double *Ta = init_array(ntriso_);
      double *Tb = init_array(ntriso_);
      double *Ga = init_array(ntriso_);
      double *Gb = init_array(ntriso_);
      double *Va = init_array(ntriso_);
      double *Vb = init_array(ntriso_);

      int soOffset = 0;
      for(int h = 0; h < nirrep_; ++h){
          for(int mu = 0; mu < nsopi_[h]; ++ mu){
              for(int nu = 0; nu <= mu; ++ nu){
                  int muNu = INDEX((nu+soOffset), (mu+soOffset));
                  Da[muNu] = kappa_so_a_->get(h, mu, nu);
                  Db[muNu] = kappa_so_b_->get(h, mu, nu);
                  Ta[muNu] = tau_so_a_->get(h,mu,nu);
                  Tb[muNu] = tau_so_b_->get(h,mu,nu);
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
    dpdbuf4 tau1_AO_aa, tau2_AO_aa;
    dpdbuf4 tau1_AO_ab, tau2_AO_ab;
    dpdbuf4 tau1_AO_bb, tau2_AO_bb;

    bool buildTensors = (options_.get_str("AO_BASIS") == "DISK");

    if(buildTensors){

        counter = 0;

        //Build the offset arrays needed for the DGEMM in half_transform
        pq_row_start = init_int_matrix(nirrep_, nirrep_);
        CD_row_start = init_int_matrix(nirrep_, nirrep_);
        cd_row_start = init_int_matrix(nirrep_, nirrep_);
        Cd_row_start = init_int_matrix(nirrep_, nirrep_);
        for(h = 0; h < nirrep_; ++h){
            for(Gc = 0, offset = 0; Gc < nirrep_; ++Gc){
                Gd = Gc ^ h;
                pq_row_start[h][Gc] = offset;
                offset += nsopi_[Gc] * nsopi_[Gd];
            }
            for(Gc = 0, offset = 0; Gc < nirrep_; ++Gc){
                Gd = Gc ^ h;
                CD_row_start[h][Gc] = offset;
                offset += navirpi_[Gc] * navirpi_[Gd];
            }
            for(Gc = 0, offset = 0; Gc < nirrep_; ++Gc){
                Gd = Gc ^ h;
                Cd_row_start[h][Gc] = offset;
                offset += navirpi_[Gc] * nbvirpi_[Gd];
            }
            for(Gc = 0, offset = 0; Gc < nirrep_; ++Gc){
                Gd = Gc ^ h;
                cd_row_start[h][Gc] = offset;
                offset += nbvirpi_[Gc] * nbvirpi_[Gd];
            }
        }

        dpd_set_default(_ints->get_dpd_id());

        /********** AA ***********/
        global_dpd_->buf4_init(&lambda, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
        global_dpd_->buf4_init(&tau1_AO_aa, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[n,n]"),
                ID("[O,O]"), ID("[n,n]"), 0, "tau1AO <OO|nn>");
        global_dpd_->buf4_scm(&tau1_AO_aa, 0.0);
        half_transform(&tau1_AO_aa, &lambda, avir_c_, avir_c_, navirpi_, navirpi_,
                pq_row_start, CD_row_start, true, 1.0, 0.0);
        global_dpd_->buf4_close(&lambda);
        global_dpd_->buf4_close(&tau1_AO_aa);

        // Now sort for better memory access patterns
        global_dpd_->buf4_init(&tau1_AO_aa, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[n,n]"),
                ID("[O,O]"), ID("[n,n]"), 1, "tau1AO <OO|nn>");
        global_dpd_->buf4_sort(&tau1_AO_aa, PSIF_DCFT_DPD, rspq, ID("[n,n]"), ID("[O,O]"), "tau1AO <nn|OO>");

        global_dpd_->buf4_close(&tau1_AO_aa);

        // Now reopen the two AO dpd_buf4's
        global_dpd_->buf4_init(&tau1_AO_aa, PSIF_DCFT_DPD, 0, ID("[n,n]"), ID("[O,O]"),
                ID("[n,n]"), ID("[O,O]"), 0, "tau1AO <nn|OO>");
        global_dpd_->buf4_init(&tau2_AO_aa, PSIF_DCFT_DPD, 0, ID("[n,n]"), ID("[O,O]"),
                ID("[n,n]"), ID("[O,O]"), 0, "tau2AO <nn|OO>");
        global_dpd_->buf4_scm(&tau2_AO_aa, 0.0);


        /********** BB ***********/
        global_dpd_->buf4_init(&lambda, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");
        global_dpd_->buf4_init(&tau1_AO_bb, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[n,n]"),
                ID("[o,o]"), ID("[n,n]"), 0, "tau1AO <oo|nn>");
        global_dpd_->buf4_scm(&tau1_AO_bb, 0.0);
        half_transform(&tau1_AO_bb, &lambda, bvir_c_, bvir_c_, nbvirpi_, nbvirpi_,
                pq_row_start, cd_row_start, true, 1.0, 0.0);
        global_dpd_->buf4_close(&lambda);
        global_dpd_->buf4_close(&tau1_AO_bb);

        // Now sort for better memory access patterns
        global_dpd_->buf4_init(&tau1_AO_bb, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[n,n]"),
                ID("[o,o]"), ID("[n,n]"), 1, "tau1AO <oo|nn>");
        global_dpd_->buf4_sort(&tau1_AO_bb, PSIF_DCFT_DPD, rspq, ID("[n,n]"), ID("[o,o]"), "tau1AO <nn|oo>");
        global_dpd_->buf4_close(&tau1_AO_bb);

        // Now reopen the two AO dpd_buf4's
        global_dpd_->buf4_init(&tau1_AO_bb, PSIF_DCFT_DPD, 0, ID("[n,n]"), ID("[o,o]"),
                ID("[n,n]"), ID("[o,o]"), 0, "tau1AO <nn|oo>");
        global_dpd_->buf4_init(&tau2_AO_bb, PSIF_DCFT_DPD, 0, ID("[n,n]"), ID("[o,o]"),
                ID("[n,n]"), ID("[o,o]"), 0, "tau2AO <nn|oo>");
        global_dpd_->buf4_scm(&tau2_AO_bb, 0.0);


        /********** AB ***********/
        global_dpd_->buf4_init(&lambda, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
            ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
        global_dpd_->buf4_init(&tau1_AO_ab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[n,n]"),
            ID("[O,o]"), ID("[n,n]"), 0, "tau1AO <Oo|nn>");
        global_dpd_->buf4_scm(&tau1_AO_ab, 0.0);
        half_transform(&tau1_AO_ab, &lambda, avir_c_, bvir_c_, navirpi_, nbvirpi_,
                pq_row_start, Cd_row_start, true, 1.0, 0.0);
        global_dpd_->buf4_close(&lambda);
        global_dpd_->buf4_close(&tau1_AO_ab);

        // Now sort for better memory access patterns
        global_dpd_->buf4_init(&tau1_AO_ab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[n,n]"),
            ID("[O,o]"), ID("[n,n]"), 0, "tau1AO <Oo|nn>");
        global_dpd_->buf4_sort(&tau1_AO_ab, PSIF_DCFT_DPD, rspq, ID("[n,n]"), ID("[O,o]"), "tau1AO <nn|Oo>");
        global_dpd_->buf4_close(&tau1_AO_ab);

        // Reopen the two AO dpd_buf4's
        global_dpd_->buf4_init(&tau1_AO_ab, PSIF_DCFT_DPD, 0, ID("[n,n]"), ID("[O,o]"),
            ID("[n,n]"), ID("[O,o]"), 0, "tau1AO <nn|Oo>");
        global_dpd_->buf4_init(&tau2_AO_ab, PSIF_DCFT_DPD, 0, ID("[n,n]"), ID("[O,o]"),
            ID("[n,n]"), ID("[O,o]"), 0, "tau2AO <nn|Oo>");
        global_dpd_->buf4_scm(&tau2_AO_ab, 0.0);



        // Now put stuff in memory
        for(int h = 0; h < nirrep_; ++h){
            global_dpd_->buf4_mat_irrep_init(&tau1_AO_aa, h);
            global_dpd_->buf4_mat_irrep_rd(&tau1_AO_aa, h);
            global_dpd_->buf4_mat_irrep_init(&tau2_AO_aa, h);

            global_dpd_->buf4_mat_irrep_init(&tau1_AO_bb, h);
            global_dpd_->buf4_mat_irrep_rd(&tau1_AO_bb, h);
            global_dpd_->buf4_mat_irrep_init(&tau2_AO_bb, h);

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
                AO_contribute(&tau1_AO_aa, &tau2_AO_aa, p, q, r, s, value);
                AO_contribute(&tau1_AO_bb, &tau2_AO_bb, p, q, r, s, value);
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
            Ga[rsArr] += (Da[pqArr] + Db[pqArr]) * value;
            Gb[rsArr] += (Da[pqArr] + Db[pqArr]) * value;
            Va[rsArr] += (Ta[pqArr] + Tb[pqArr]) * value;
            Vb[rsArr] += (Ta[pqArr] + Tb[pqArr]) * value;
            if(q >= r){
                Ga[qrArr] -= Da[psArr] * value;
                Gb[qrArr] -= Db[psArr] * value;
                Va[qrArr] -= Ta[psArr] * value;
                Vb[qrArr] -= Tb[psArr] * value;
            }

            if(p!=q && r!=s && pqArr!=rsArr){
                /* (pq|sr) */
                if(s >= r){
                    Ga[srArr] += (Da[pqArr] + Db[pqArr]) * value;
                    Gb[srArr] += (Da[pqArr] + Db[pqArr]) * value;
                    Va[srArr] += (Ta[pqArr] + Tb[pqArr]) * value;
                    Vb[srArr] += (Ta[pqArr] + Tb[pqArr]) * value;
                }
                if(q >= s){
                    Ga[qsArr] -= Da[prArr] * value;
                    Gb[qsArr] -= Db[prArr] * value;
                    Va[qsArr] -= Ta[prArr] * value;
                    Vb[qsArr] -= Tb[prArr] * value;
                }

                /* (qp|rs) */
                if(r >= s){
                    Ga[rsArr] += (Da[qpArr] + Db[qpArr]) * value;
                    Gb[rsArr] += (Da[qpArr] + Db[qpArr]) * value;
                    Va[rsArr] += (Ta[qpArr] + Tb[qpArr]) * value;
                    Vb[rsArr] += (Ta[qpArr] + Tb[qpArr]) * value;
                }
                if(p >= r){
                    Ga[prArr] -= Da[qsArr] * value;
                    Gb[prArr] -= Db[qsArr] * value;
                    Va[prArr] -= Ta[qsArr] * value;
                    Vb[prArr] -= Tb[qsArr] * value;
                }

                /* (qp|sr) */
                if(s >= r){
                    Ga[srArr] += (Da[qpArr] + Db[qpArr]) * value;
                    Gb[srArr] += (Da[qpArr] + Db[qpArr]) * value;
                    Va[srArr] += (Ta[qpArr] + Tb[qpArr]) * value;
                    Vb[srArr] += (Ta[qpArr] + Tb[qpArr]) * value;
                }
                if(p >= s){
                    Ga[psArr] -= Da[qrArr] * value;
                    Gb[psArr] -= Db[qrArr] * value;
                    Va[psArr] -= Ta[qrArr] * value;
                    Vb[psArr] -= Tb[qrArr] * value;
                }

                /* (rs|pq) */
                if(p >= q){
                    Ga[pqArr] += (Da[rsArr] + Db[rsArr]) * value;
                    Gb[pqArr] += (Da[rsArr] + Db[rsArr]) * value;
                    Va[pqArr] += (Ta[rsArr] + Tb[rsArr]) * value;
                    Vb[pqArr] += (Ta[rsArr] + Tb[rsArr]) * value;
                }
                if(s >= p){
                    Ga[spArr] -= Da[rqArr] * value;
                    Gb[spArr] -= Db[rqArr] * value;
                    Va[spArr] -= Ta[rqArr] * value;
                    Vb[spArr] -= Tb[rqArr] * value;
                }

                /* (sr|pq) */
                if(p >= q){
                    Ga[pqArr] += (Da[srArr] + Db[srArr]) * value;
                    Gb[pqArr] += (Da[srArr] + Db[srArr]) * value;
                    Va[pqArr] += (Ta[srArr] + Tb[srArr]) * value;
                    Vb[pqArr] += (Ta[srArr] + Tb[srArr]) * value;
                }
                if(r >= p){
                    Ga[rpArr] -= Da[sqArr] * value;
                    Gb[rpArr] -= Db[sqArr] * value;
                    Va[rpArr] -= Ta[sqArr] * value;
                    Vb[rpArr] -= Tb[sqArr] * value;
                }

                /* (rs|qp) */
                if(q >= p){
                    Ga[qpArr] += (Da[rsArr] + Db[rsArr]) * value;
                    Gb[qpArr] += (Da[rsArr] + Db[rsArr]) * value;
                    Va[qpArr] += (Ta[rsArr] + Tb[rsArr]) * value;
                    Vb[qpArr] += (Ta[rsArr] + Tb[rsArr]) * value;
                }
                if(s >= q){
                    Ga[sqArr] -= Da[rpArr] * value;
                    Gb[sqArr] -= Db[rpArr] * value;
                    Va[sqArr] -= Ta[rpArr] * value;
                    Vb[sqArr] -= Tb[rpArr] * value;
                }

                /* (sr|qp) */
                if(q >= p){
                    Ga[qpArr] += (Da[srArr] + Db[srArr]) * value;
                    Gb[qpArr] += (Da[srArr] + Db[srArr]) * value;
                    Va[qpArr] += (Ta[srArr] + Tb[srArr]) * value;
                    Vb[qpArr] += (Ta[srArr] + Tb[srArr]) * value;
                }
                if(r >= q){
                    Ga[rqArr] -= Da[spArr] * value;
                    Gb[rqArr] -= Db[spArr] * value;
                    Va[rqArr] -= Ta[spArr] * value;
                    Vb[rqArr] -= Tb[spArr] * value;
                }
            }else if(p!=q && r!=s && pqArr==rsArr){
                /* (pq|sr) */
                if(s >= r){
                    Ga[srArr] += (Da[pqArr] + Db[pqArr]) * value;
                    Gb[srArr] += (Da[pqArr] + Db[pqArr]) * value;
                    Va[srArr] += (Ta[pqArr] + Tb[pqArr]) * value;
                    Vb[srArr] += (Ta[pqArr] + Tb[pqArr]) * value;
                }
                if(q >= s){
                    Ga[qsArr] -= Da[prArr] * value;
                    Gb[qsArr] -= Db[prArr] * value;
                    Va[qsArr] -= Ta[prArr] * value;
                    Vb[qsArr] -= Tb[prArr] * value;
                }
                /* (qp|rs) */
                if(r >= s){
                    Ga[rsArr] += (Da[qpArr] + Db[qpArr]) * value;
                    Gb[rsArr] += (Da[qpArr] + Db[qpArr]) * value;
                    Va[rsArr] += (Ta[qpArr] + Tb[qpArr]) * value;
                    Vb[rsArr] += (Ta[qpArr] + Tb[qpArr]) * value;
                }
                if(p >= r){
                    Ga[prArr] -= Da[qsArr] * value;
                    Gb[prArr] -= Db[qsArr] * value;
                    Va[prArr] -= Ta[qsArr] * value;
                    Vb[prArr] -= Tb[qsArr] * value;
                }

                /* (qp|sr) */
                if(s >= r){
                    Ga[srArr] += (Da[qpArr] + Db[qpArr]) * value;
                    Gb[srArr] += (Da[qpArr] + Db[qpArr]) * value;
                    Va[srArr] += (Ta[qpArr] + Tb[qpArr]) * value;
                    Vb[srArr] += (Ta[qpArr] + Tb[qpArr]) * value;
                }
                if(p >= s){
                    Ga[psArr] -= Da[qrArr] * value;
                    Gb[psArr] -= Db[qrArr] * value;
                    Va[psArr] -= Ta[qrArr] * value;
                    Vb[psArr] -= Tb[qrArr] * value;
                }
            }else if(p!=q && r==s){
                /* (qp|rs) */
                if(r >= s){
                    Ga[rsArr] += (Da[qpArr] + Db[qpArr]) * value;
                    Gb[rsArr] += (Da[qpArr] + Db[qpArr]) * value;
                    Va[rsArr] += (Ta[qpArr] + Tb[qpArr]) * value;
                    Vb[rsArr] += (Ta[qpArr] + Tb[qpArr]) * value;
                }
                if(p >= r){
                    Ga[prArr] -= Da[qsArr] * value;
                    Gb[prArr] -= Db[qsArr] * value;
                    Va[prArr] -= Ta[qsArr] * value;
                    Vb[prArr] -= Tb[qsArr] * value;
                }

                /* (rs|pq) */
                if(p >= q){
                    Ga[pqArr] += (Da[rsArr] + Db[rsArr]) * value;
                    Gb[pqArr] += (Da[rsArr] + Db[rsArr]) * value;
                    Va[pqArr] += (Ta[rsArr] + Tb[rsArr]) * value;
                    Vb[pqArr] += (Ta[rsArr] + Tb[rsArr]) * value;
                }
                if(s >= p){
                    Ga[spArr] -= Da[rqArr] * value;
                    Gb[spArr] -= Db[rqArr] * value;
                    Va[spArr] -= Ta[rqArr] * value;
                    Vb[spArr] -= Tb[rqArr] * value;
                }

                /* (rs|qp) */
                if(q >= p){
                    Ga[qpArr] += (Da[rsArr] + Db[rsArr]) * value;
                    Gb[qpArr] += (Da[rsArr] + Db[rsArr]) * value;
                    Va[qpArr] += (Ta[rsArr] + Tb[rsArr]) * value;
                    Vb[qpArr] += (Ta[rsArr] + Tb[rsArr]) * value;
                }
                if(s >= q){
                    Ga[sqArr] -= Da[rpArr] * value;
                    Gb[sqArr] -= Db[rpArr] * value;
                    Va[sqArr] -= Ta[rpArr] * value;
                    Vb[sqArr] -= Tb[rpArr] * value;
                }
            }else if(p==q && r!=s){
                /* (pq|sr) */
                if(s >= r){
                    Ga[srArr] += (Da[pqArr] + Db[pqArr]) * value;
                    Gb[srArr] += (Da[pqArr] + Db[pqArr]) * value;
                    Va[srArr] += (Ta[pqArr] + Tb[pqArr]) * value;
                    Vb[srArr] += (Ta[pqArr] + Tb[pqArr]) * value;
                }
                if(q >= s){
                    Ga[qsArr] -= Da[prArr] * value;
                    Gb[qsArr] -= Db[prArr] * value;
                    Va[qsArr] -= Ta[prArr] * value;
                    Vb[qsArr] -= Tb[prArr] * value;
                }

                /* (rs|pq) */
                if(p >= q){
                    Ga[pqArr] += (Da[rsArr] + Db[rsArr]) * value;
                    Gb[pqArr] += (Da[rsArr] + Db[rsArr]) * value;
                    Va[pqArr] += (Ta[rsArr] + Tb[rsArr]) * value;
                    Vb[pqArr] += (Ta[rsArr] + Tb[rsArr]) * value;
                }
                if(s >= p){
                    Ga[spArr] -= Da[rqArr] * value;
                    Gb[spArr] -= Db[rqArr] * value;
                    Va[spArr] -= Ta[rqArr] * value;
                    Vb[spArr] -= Tb[rqArr] * value;
                }

                /* (sr|pq) */
                if(p >= q){
                    Ga[pqArr] += (Da[srArr] + Db[srArr]) * value;
                    Gb[pqArr] += (Da[srArr] + Db[srArr]) * value;
                    Va[pqArr] += (Ta[srArr] + Tb[srArr]) * value;
                    Vb[pqArr] += (Ta[srArr] + Tb[srArr]) * value;
                }
                if(r >= p){
                    Ga[rpArr] -= Da[sqArr] * value;
                    Gb[rpArr] -= Db[sqArr] * value;
                    Va[rpArr] -= Ta[sqArr] * value;
                    Vb[rpArr] -= Tb[sqArr] * value;
                }
            }else if(p==q && r==s && pqArr!=rsArr){
                /* (rs|pq) */
                if(p >= q){
                    Ga[pqArr] += (Da[rsArr] + Db[rsArr]) * value;
                    Gb[pqArr] += (Da[rsArr] + Db[rsArr]) * value;
                    Va[pqArr] += (Ta[rsArr] + Tb[rsArr]) * value;
                    Vb[pqArr] += (Ta[rsArr] + Tb[rsArr]) * value;
                }
                if(s >= p){
                    Ga[spArr] -= Da[rqArr] * value;
                    Gb[spArr] -= Db[rqArr] * value;
                    Va[spArr] -= Ta[rqArr] * value;
                    Vb[spArr] -= Tb[rqArr] * value;
                }
            }
        } /* end loop through current buffer */
        if(!lastBuffer) iwl->fetch();
    }while(!lastBuffer);
    iwl->set_keep_flag(1);
    delete iwl;
    if(buildTensors){
        if(print_ > 1){
            outfile->Printf( "Processed %d SO integrals each for AA, BB, and AB\n", counter);
        }
        for(int h = 0; h < nirrep_; ++h){
            global_dpd_->buf4_mat_irrep_wrt(&tau2_AO_aa, h);
            global_dpd_->buf4_mat_irrep_close(&tau1_AO_aa, h);
            global_dpd_->buf4_mat_irrep_close(&tau2_AO_aa, h);

            global_dpd_->buf4_mat_irrep_wrt(&tau2_AO_bb, h);
            global_dpd_->buf4_mat_irrep_close(&tau1_AO_bb, h);
            global_dpd_->buf4_mat_irrep_close(&tau2_AO_bb, h);

            global_dpd_->buf4_mat_irrep_wrt(&tau2_AO_ab, h);
            global_dpd_->buf4_mat_irrep_close(&tau1_AO_ab, h);
            global_dpd_->buf4_mat_irrep_close(&tau2_AO_ab, h);
        }

        global_dpd_->buf4_close(&tau1_AO_aa);
        global_dpd_->buf4_close(&tau1_AO_bb);
        global_dpd_->buf4_close(&tau1_AO_ab);
        global_dpd_->buf4_close(&tau2_AO_aa);
        global_dpd_->buf4_close(&tau2_AO_bb);
        global_dpd_->buf4_close(&tau2_AO_ab);

        /********** AA ***********/
        global_dpd_->buf4_init(&tau2_AO_aa, PSIF_DCFT_DPD, 0, ID("[n,n]"), ID("[O,O]"),
               ID("[n,n]"), ID("[O,O]"), 0, "tau2AO <nn|OO>");
        global_dpd_->buf4_sort(&tau2_AO_aa, PSIF_DCFT_DPD, rspq, ID("[O,O]"), ID("[n,n]"), "tau2AO <OO|nn>");
        global_dpd_->buf4_close(&tau2_AO_aa);
        global_dpd_->buf4_init(&tau2_AO_aa, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[n,n]"),
               ID("[O,O]"), ID("[n,n]"), 0, "tau2AO <OO|nn>");
        global_dpd_->buf4_init(&tau_temp, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
               ID("[O,O]"), ID("[V,V]"), 0, "tau(temp) <OO|VV>");
        global_dpd_->buf4_scm(&tau_temp, 0.0);
        half_transform(&tau2_AO_aa, &tau_temp, avir_c_, avir_c_, navirpi_, navirpi_,
                pq_row_start, CD_row_start, false, 0.5, 0.0);
        global_dpd_->buf4_close(&tau2_AO_aa);
        global_dpd_->buf4_close(&tau_temp);

        /********** BB ***********/
        global_dpd_->buf4_init(&tau2_AO_bb, PSIF_DCFT_DPD, 0, ID("[n,n]"), ID("[o,o]"),
               ID("[n,n]"), ID("[o,o]"), 0, "tau2AO <nn|oo>");
        global_dpd_->buf4_sort(&tau2_AO_bb, PSIF_DCFT_DPD, rspq, ID("[o,o]"), ID("[n,n]"), "tau2AO <oo|nn>");
        global_dpd_->buf4_close(&tau2_AO_bb);
        global_dpd_->buf4_init(&tau2_AO_bb, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[n,n]"),
               ID("[o,o]"), ID("[n,n]"), 0, "tau2AO <oo|nn>");
        global_dpd_->buf4_init(&tau_temp, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                ID("[o,o]"), ID("[v,v]"), 0, "tau(temp) <oo|vv>");
        global_dpd_->buf4_scm(&tau_temp, 0.0);
        half_transform(&tau2_AO_bb, &tau_temp, bvir_c_, bvir_c_, nbvirpi_, nbvirpi_,
                pq_row_start, cd_row_start, false, 0.5, 0.0);
        global_dpd_->buf4_close(&tau2_AO_bb);
        global_dpd_->buf4_close(&tau_temp);

        /********** AB ***********/
        global_dpd_->buf4_init(&tau2_AO_ab, PSIF_DCFT_DPD, 0, ID("[n,n]"), ID("[O,o]"),
                ID("[n,n]"), ID("[O,o]"), 0, "tau2AO <nn|Oo>");
        global_dpd_->buf4_sort(&tau2_AO_ab, PSIF_DCFT_DPD, rspq, ID("[O,o]"), ID("[n,n]"), "tau2AO <Oo|nn>");
        global_dpd_->buf4_close(&tau2_AO_ab);
        global_dpd_->buf4_init(&tau2_AO_ab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[n,n]"),
               ID("[O,o]"), ID("[n,n]"), 0, "tau2AO <Oo|nn>");
        global_dpd_->buf4_init(&tau_temp, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                ID("[O,o]"), ID("[V,v]"), 0, "tau(temp) <Oo|Vv>");
        global_dpd_->buf4_scm(&tau_temp, 0.0);
        half_transform(&tau2_AO_ab, &tau_temp, avir_c_, bvir_c_, navirpi_, nbvirpi_,
                pq_row_start, Cd_row_start, false, 1.0, 0.0);
        global_dpd_->buf4_close(&tau2_AO_ab);
        global_dpd_->buf4_close(&tau_temp);


        free_int_matrix(pq_row_start);
        free_int_matrix(CD_row_start);
        free_int_matrix(cd_row_start);
        free_int_matrix(Cd_row_start);


    }



    // Build the Fock matrices from the H and G matrices
    soOffset = 0;
    for(int h = 0; h < nirrep_; ++h){
        for(int mu = 0; mu < nsopi_[h]; ++mu){
            for(int nu = 0; nu <= mu; ++nu){
                int muNu = INDEX((nu+soOffset), (mu+soOffset));
                double aVal   = Ga[muNu];
                double bVal   = Gb[muNu];
                double aGTVal = Va[muNu];
                double bGTVal = Vb[muNu];
                Fa_->add(h, mu, nu, aVal);
                Fb_->add(h, mu, nu, bVal);
                g_tau_a_->set(h, mu, nu, aGTVal);
                g_tau_b_->set(h, mu, nu, bGTVal);
                if(mu != nu){
                    Fa_->add(h, nu, mu, aVal);
                    Fb_->add(h, nu, mu, bVal);
                    g_tau_a_->set(h, nu, mu, aGTVal);
                    g_tau_b_->set(h, nu, mu, bGTVal);
                }
            }
        }
        soOffset += nsopi_[h];
    }

    delete [] Ta;
    delete [] Tb;
    delete [] Va;
    delete [] Vb;
    delete [] Da;
    delete [] Db;
    delete [] Ga;
    delete [] Gb;

    dcft_timer_off("DCFTSolver::process_so_ints");
}

  /*
  * Updates the Fock operator every lambda iteration using new Tau
  */
  void
  DCFTSolver::update_fock()
  {
      dcft_timer_on("DCFTSolver::update_fock");

      dpdfile2 Gtau;

      moFa_->copy(moF0a_);
      moFb_->copy(moF0b_);

      // Copy MO basis GTau to the memory

      psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

      // Alpha occupied
      global_dpd_->file2_init(&Gtau, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "GTau <O|O>");
      global_dpd_->file2_mat_init(&Gtau);
      global_dpd_->file2_mat_rd(&Gtau);
      for(int h = 0; h < nirrep_; ++h){
          for(int i = 0 ; i < naoccpi_[h]; ++i){
              for(int j = 0 ; j < naoccpi_[h]; ++j){
                  moG_tau_a_->set(h, frzcpi_[h] + i, frzcpi_[h] + j, Gtau.matrix[h][i][j]);
              }
          }
      }
      global_dpd_->file2_mat_close(&Gtau);
      global_dpd_->file2_close(&Gtau);

      // Alpha virtual
      global_dpd_->file2_init(&Gtau, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "GTau <V|V>");
      global_dpd_->file2_mat_init(&Gtau);
      global_dpd_->file2_mat_rd(&Gtau);
      for(int h = 0; h < nirrep_; ++h){
          for(int i = 0 ; i < navirpi_[h]; ++i){
              for(int j = 0; j < navirpi_[h]; ++j){
                  moG_tau_a_->set(h, naoccpi_[h] + i, naoccpi_[h] + j, Gtau.matrix[h][i][j]);
              }
          }
      }
      global_dpd_->file2_mat_close(&Gtau);
      global_dpd_->file2_close(&Gtau);

      // Beta occupied
      global_dpd_->file2_init(&Gtau, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "GTau <o|o>");
      global_dpd_->file2_mat_init(&Gtau);
      global_dpd_->file2_mat_rd(&Gtau);
      for(int h = 0; h < nirrep_; ++h){
          for(int i = 0 ; i < nboccpi_[h]; ++i){
              for(int j = 0 ; j < nboccpi_[h]; ++j){
                  moG_tau_b_->set(h, frzcpi_[h] + i, frzcpi_[h] + j, Gtau.matrix[h][i][j]);
              }
          }
      }
      global_dpd_->file2_mat_close(&Gtau);
      global_dpd_->file2_close(&Gtau);

      // Beta virtual
      global_dpd_->file2_init(&Gtau, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "GTau <v|v>");
      global_dpd_->file2_mat_init(&Gtau);
      global_dpd_->file2_mat_rd(&Gtau);
      for(int h = 0; h < nirrep_; ++h){
          for(int i = 0; i < nbvirpi_[h]; ++i){
              for(int j = 0; j < nbvirpi_[h]; ++j){
                  moG_tau_b_->set(h, nboccpi_[h] + i, nboccpi_[h] + j, Gtau.matrix[h][i][j]);
              }
          }
      }
      global_dpd_->file2_mat_close(&Gtau);
      global_dpd_->file2_close(&Gtau);

      // Add the GTau contribution to the Fock operator
      moFa_->add(moG_tau_a_);
      moFb_->add(moG_tau_b_);

      // Write the MO basis Fock operator to the DPD file and update the denominators
      build_denominators();

      psio_->close(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

      dcft_timer_off("DCFTSolver::update_fock");
  }


  /**
  * Builds the AO basis tensors for the <VV||VV>, <vv||vv>, and <Vv|Vv> terms in the G and X intermediates
  */
  void
  DCFTSolver::build_AO_tensors()
  {
      dcft_timer_on("DCFTSolver::build_tensors");

      IWL *iwl = new IWL(psio_.get(), PSIF_SO_TEI, int_tolerance_, 1, 1);

      Label *lblptr = iwl->labels();
      Value *valptr = iwl->values();

      double value;
      int Gc, Gd;
      int offset, labelIndex, p, q, r, s, h, counter;
      int **pq_row_start, **CD_row_start, **Cd_row_start, **cd_row_start;
      dpdbuf4 tau_temp, lambda;
      dpdbuf4 tau1_AO_aa, tau2_AO_aa;
      dpdbuf4 tau1_AO_ab, tau2_AO_ab;
      dpdbuf4 tau1_AO_bb, tau2_AO_bb;
      dpdfile2 s_aa_1, s_aa_2, s_aa_3, s_aa_4, tau;
      dpdfile2 s_bb_1, s_bb_2, s_bb_3, s_bb_4;

      counter = 0;

      //Build the offset arrays needed for the DGEMM in half_transform
      pq_row_start = init_int_matrix(nirrep_, nirrep_);
      CD_row_start = init_int_matrix(nirrep_, nirrep_);
      cd_row_start = init_int_matrix(nirrep_, nirrep_);
      Cd_row_start = init_int_matrix(nirrep_, nirrep_);
      for(h = 0; h < nirrep_; ++h){
          for(Gc = 0, offset = 0; Gc < nirrep_; ++Gc){
              Gd = Gc ^ h;
              pq_row_start[h][Gc] = offset;
              offset += nsopi_[Gc] * nsopi_[Gd];
          }
          for(Gc = 0, offset = 0; Gc < nirrep_; ++Gc){
              Gd = Gc ^ h;
              CD_row_start[h][Gc] = offset;
              offset += navirpi_[Gc] * navirpi_[Gd];
          }
          for(Gc = 0, offset = 0; Gc < nirrep_; ++Gc){
              Gd = Gc ^ h;
              Cd_row_start[h][Gc] = offset;
              offset += navirpi_[Gc] * nbvirpi_[Gd];
          }
          for(Gc = 0, offset = 0; Gc < nirrep_; ++Gc){
              Gd = Gc ^ h;
              cd_row_start[h][Gc] = offset;
              offset += nbvirpi_[Gc] * nbvirpi_[Gd];
          }
      }

      dpd_set_default(_ints->get_dpd_id());

      /********** AA ***********/
      global_dpd_->buf4_init(&lambda, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                    ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
      global_dpd_->buf4_init(&tau1_AO_aa, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[n,n]"),
                    ID("[O,O]"), ID("[n,n]"), 0, "tau1AO <OO|nn>");
      global_dpd_->buf4_scm(&tau1_AO_aa, 0.0);
      half_transform(&tau1_AO_aa, &lambda, avir_c_, avir_c_, navirpi_, navirpi_,
                     pq_row_start, CD_row_start, true, 1.0, 0.0);
      global_dpd_->buf4_close(&lambda);
      global_dpd_->buf4_close(&tau1_AO_aa);

      // Now sort for better memory access patterns
      global_dpd_->buf4_init(&tau1_AO_aa, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[n,n]"),
                    ID("[O,O]"), ID("[n,n]"), 1, "tau1AO <OO|nn>");
      global_dpd_->buf4_sort(&tau1_AO_aa, PSIF_DCFT_DPD, rspq, ID("[n,n]"), ID("[O,O]"), "tau1AO <nn|OO>");

      global_dpd_->buf4_close(&tau1_AO_aa);

      // Now reopen the two AO dpd_buf4's
      global_dpd_->buf4_init(&tau1_AO_aa, PSIF_DCFT_DPD, 0, ID("[n,n]"), ID("[O,O]"),
                    ID("[n,n]"), ID("[O,O]"), 0, "tau1AO <nn|OO>");
      global_dpd_->buf4_init(&tau2_AO_aa, PSIF_DCFT_DPD, 0, ID("[n,n]"), ID("[O,O]"),
                    ID("[n,n]"), ID("[O,O]"), 0, "tau2AO <nn|OO>");
      global_dpd_->buf4_scm(&tau2_AO_aa, 0.0);

      /********** BB ***********/
      global_dpd_->buf4_init(&lambda, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                    ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");
      global_dpd_->buf4_init(&tau1_AO_bb, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[n,n]"),
                    ID("[o,o]"), ID("[n,n]"), 0, "tau1AO <oo|nn>");
      global_dpd_->buf4_scm(&tau1_AO_bb, 0.0);
      half_transform(&tau1_AO_bb, &lambda, bvir_c_, bvir_c_, nbvirpi_, nbvirpi_,
                     pq_row_start, cd_row_start, true, 1.0, 0.0);
      global_dpd_->buf4_close(&lambda);
      global_dpd_->buf4_close(&tau1_AO_bb);

      // Now sort for better memory access patterns
      global_dpd_->buf4_init(&tau1_AO_bb, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[n,n]"),
                    ID("[o,o]"), ID("[n,n]"), 1, "tau1AO <oo|nn>");
      global_dpd_->buf4_sort(&tau1_AO_bb, PSIF_DCFT_DPD, rspq, ID("[n,n]"), ID("[o,o]"), "tau1AO <nn|oo>");
      global_dpd_->buf4_close(&tau1_AO_bb);

      // Now reopen the two AO dpd_buf4's
      global_dpd_->buf4_init(&tau1_AO_bb, PSIF_DCFT_DPD, 0, ID("[n,n]"), ID("[o,o]"),
                    ID("[n,n]"), ID("[o,o]"), 0, "tau1AO <nn|oo>");
      global_dpd_->buf4_init(&tau2_AO_bb, PSIF_DCFT_DPD, 0, ID("[n,n]"), ID("[o,o]"),
                    ID("[n,n]"), ID("[o,o]"), 0, "tau2AO <nn|oo>");
      global_dpd_->buf4_scm(&tau2_AO_bb, 0.0);

      /********** AB ***********/
      global_dpd_->buf4_init(&lambda, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                    ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
      global_dpd_->buf4_init(&tau1_AO_ab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[n,n]"),
                    ID("[O,o]"), ID("[n,n]"), 0, "tau1AO <Oo|nn>");
      global_dpd_->buf4_scm(&tau1_AO_ab, 0.0);
      half_transform(&tau1_AO_ab, &lambda, avir_c_, bvir_c_, navirpi_, nbvirpi_,
                     pq_row_start, Cd_row_start, true, 1.0, 0.0);
      global_dpd_->buf4_close(&lambda);
      global_dpd_->buf4_close(&tau1_AO_ab);

      // Now sort for better memory access patterns
      global_dpd_->buf4_init(&tau1_AO_ab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[n,n]"),
                    ID("[O,o]"), ID("[n,n]"), 0, "tau1AO <Oo|nn>");
      global_dpd_->buf4_sort(&tau1_AO_ab, PSIF_DCFT_DPD, rspq, ID("[n,n]"), ID("[O,o]"), "tau1AO <nn|Oo>");
      global_dpd_->buf4_close(&tau1_AO_ab);

      // Reopen the two AO dpd_buf4's
      global_dpd_->buf4_init(&tau1_AO_ab, PSIF_DCFT_DPD, 0, ID("[n,n]"), ID("[O,o]"),
                    ID("[n,n]"), ID("[O,o]"), 0, "tau1AO <nn|Oo>");
      global_dpd_->buf4_init(&tau2_AO_ab, PSIF_DCFT_DPD, 0, ID("[n,n]"), ID("[O,o]"),
                    ID("[n,n]"), ID("[O,o]"), 0, "tau2AO <nn|Oo>");
      global_dpd_->buf4_scm(&tau2_AO_ab, 0.0);

      // Preparing Tau for the GTau terms
      // This is where the SO-basis Tau will be
      global_dpd_->file2_init(&s_aa_1, PSIF_DCFT_DPD, 0, ID('n'), ID('n'), "s1(temp)A <n|n>");
      global_dpd_->file2_init(&s_bb_1, PSIF_DCFT_DPD, 0, ID('n'), ID('n'), "s1(temp)B <n|n>");

      // Write SO-basis Tau to disk
      global_dpd_->file2_mat_init(&s_aa_1);
      global_dpd_->file2_mat_init(&s_bb_1);
      for(int h = 0; h < nirrep_; ++h){
          for(int i = 0 ; i < nsopi_[h]; ++i){
              for(int j = 0 ; j < nsopi_[h]; ++j){
                  s_aa_1.matrix[h][i][j] = tau_so_a_->get(h, i, j);
                  s_bb_1.matrix[h][i][j] = tau_so_b_->get(h, i, j);
              }
          }
      }
      global_dpd_->file2_mat_wrt(&s_aa_1);
      global_dpd_->file2_mat_wrt(&s_bb_1);
      global_dpd_->file2_close(&s_aa_1);
      global_dpd_->file2_close(&s_bb_1);

      // Reopen the arrays
      global_dpd_->file2_init(&s_aa_1, PSIF_DCFT_DPD, 0, ID('n'), ID('n'), "s1(temp)A <n|n>");
      global_dpd_->file2_init(&s_bb_1, PSIF_DCFT_DPD, 0, ID('n'), ID('n'), "s1(temp)B <n|n>");

      // This is where GTau contribution will be placed
      global_dpd_->file2_init(&s_aa_2, PSIF_DCFT_DPD, 0, ID('n'), ID('n'), "s2(temp)A <n|n>");
      global_dpd_->file2_init(&s_bb_2, PSIF_DCFT_DPD, 0, ID('n'), ID('n'), "s2(temp)B <n|n>");
      global_dpd_->file2_scm(&s_aa_2, 0.0);
      global_dpd_->file2_scm(&s_bb_2, 0.0);

      // Now put stuff in memory
      global_dpd_->file2_mat_init(&s_aa_1);
      global_dpd_->file2_mat_init(&s_aa_2);
      global_dpd_->file2_mat_rd(&s_aa_1);
      global_dpd_->file2_mat_rd(&s_aa_2);
      global_dpd_->file2_mat_init(&s_bb_1);
      global_dpd_->file2_mat_init(&s_bb_2);
      global_dpd_->file2_mat_rd(&s_bb_1);
      global_dpd_->file2_mat_rd(&s_bb_2);

      for(int h = 0; h < nirrep_; ++h){
          global_dpd_->buf4_mat_irrep_init(&tau1_AO_aa, h);
          global_dpd_->buf4_mat_irrep_rd(&tau1_AO_aa, h);
          global_dpd_->buf4_mat_irrep_init(&tau2_AO_aa, h);

          global_dpd_->buf4_mat_irrep_init(&tau1_AO_bb, h);
          global_dpd_->buf4_mat_irrep_rd(&tau1_AO_bb, h);
          global_dpd_->buf4_mat_irrep_init(&tau2_AO_bb, h);

          global_dpd_->buf4_mat_irrep_init(&tau1_AO_ab, h);
          global_dpd_->buf4_mat_irrep_rd(&tau1_AO_ab, h);
          global_dpd_->buf4_mat_irrep_init(&tau2_AO_ab, h);

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
              AO_contribute(&tau1_AO_aa, &tau2_AO_aa, p, q, r, s, value, &s_aa_1, &s_bb_1, &s_aa_2);
              AO_contribute(&tau1_AO_bb, &tau2_AO_bb, p, q, r, s, value, &s_bb_1, &s_aa_1, &s_bb_2);
              AO_contribute(&tau1_AO_ab, &tau2_AO_ab, p, q, r, s, value);
              ++counter;

          } /* end loop through current buffer */
          if(!lastBuffer) iwl->fetch();
      }while(!lastBuffer);
      iwl->set_keep_flag(1);
      delete iwl;
      if(print_ > 1){
          outfile->Printf( "Processed %d SO integrals each for AA, BB, and AB\n", counter);
      }
      for(int h = 0; h < nirrep_; ++h){
          global_dpd_->buf4_mat_irrep_wrt(&tau2_AO_aa, h);
          global_dpd_->buf4_mat_irrep_close(&tau1_AO_aa, h);
          global_dpd_->buf4_mat_irrep_close(&tau2_AO_aa, h);

          global_dpd_->buf4_mat_irrep_wrt(&tau2_AO_bb, h);
          global_dpd_->buf4_mat_irrep_close(&tau1_AO_bb, h);
          global_dpd_->buf4_mat_irrep_close(&tau2_AO_bb, h);

          global_dpd_->buf4_mat_irrep_wrt(&tau2_AO_ab, h);
          global_dpd_->buf4_mat_irrep_close(&tau1_AO_ab, h);
          global_dpd_->buf4_mat_irrep_close(&tau2_AO_ab, h);
      }

      global_dpd_->buf4_close(&tau1_AO_aa);
      global_dpd_->buf4_close(&tau1_AO_bb);
      global_dpd_->buf4_close(&tau1_AO_ab);
      global_dpd_->buf4_close(&tau2_AO_aa);
      global_dpd_->buf4_close(&tau2_AO_bb);
      global_dpd_->buf4_close(&tau2_AO_ab);
      global_dpd_->file2_mat_wrt(&s_aa_2);
      global_dpd_->file2_mat_wrt(&s_bb_2);
      global_dpd_->file2_mat_close(&s_aa_2);
      global_dpd_->file2_mat_close(&s_bb_2);


      /********** AA ***********/
      global_dpd_->buf4_init(&tau2_AO_aa, PSIF_DCFT_DPD, 0, ID("[n,n]"), ID("[O,O]"),
                    ID("[n,n]"), ID("[O,O]"), 0, "tau2AO <nn|OO>");
      global_dpd_->buf4_sort(&tau2_AO_aa, PSIF_DCFT_DPD, rspq, ID("[O,O]"), ID("[n,n]"), "tau2AO <OO|nn>");
      global_dpd_->buf4_close(&tau2_AO_aa);
      global_dpd_->buf4_init(&tau2_AO_aa, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[n,n]"),
                    ID("[O,O]"), ID("[n,n]"), 0, "tau2AO <OO|nn>");
      global_dpd_->buf4_init(&tau_temp, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                    ID("[O,O]"), ID("[V,V]"), 0, "tau(temp) <OO|VV>");
      global_dpd_->buf4_scm(&tau_temp, 0.0);
      half_transform(&tau2_AO_aa, &tau_temp, avir_c_, avir_c_, navirpi_, navirpi_,
                     pq_row_start, CD_row_start, false, 0.5, 0.0);
      global_dpd_->buf4_close(&tau2_AO_aa);
      global_dpd_->buf4_close(&tau_temp);


      /********** BB ***********/
      global_dpd_->buf4_init(&tau2_AO_bb, PSIF_DCFT_DPD, 0, ID("[n,n]"), ID("[o,o]"),
                    ID("[n,n]"), ID("[o,o]"), 0, "tau2AO <nn|oo>");
      global_dpd_->buf4_sort(&tau2_AO_bb, PSIF_DCFT_DPD, rspq, ID("[o,o]"), ID("[n,n]"), "tau2AO <oo|nn>");
      global_dpd_->buf4_close(&tau2_AO_bb);
      global_dpd_->buf4_init(&tau2_AO_bb, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[n,n]"),
                    ID("[o,o]"), ID("[n,n]"), 0, "tau2AO <oo|nn>");
      global_dpd_->buf4_init(&tau_temp, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                    ID("[o,o]"), ID("[v,v]"), 0, "tau(temp) <oo|vv>");
      global_dpd_->buf4_scm(&tau_temp, 0.0);
      half_transform(&tau2_AO_bb, &tau_temp, bvir_c_, bvir_c_, nbvirpi_, nbvirpi_,
                     pq_row_start, cd_row_start, false, 0.5, 0.0);
      global_dpd_->buf4_close(&tau2_AO_bb);
      global_dpd_->buf4_close(&tau_temp);


      /********** AB ***********/
      global_dpd_->buf4_init(&tau2_AO_ab, PSIF_DCFT_DPD, 0, ID("[n,n]"), ID("[O,o]"),
                    ID("[n,n]"), ID("[O,o]"), 0, "tau2AO <nn|Oo>");
      global_dpd_->buf4_sort(&tau2_AO_ab, PSIF_DCFT_DPD, rspq, ID("[O,o]"), ID("[n,n]"), "tau2AO <Oo|nn>");
      global_dpd_->buf4_close(&tau2_AO_ab);
      global_dpd_->buf4_init(&tau2_AO_ab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[n,n]"),
                    ID("[O,o]"), ID("[n,n]"), 0, "tau2AO <Oo|nn>");
      global_dpd_->buf4_init(&tau_temp, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                    ID("[O,o]"), ID("[V,v]"), 0, "tau(temp) <Oo|Vv>");
      global_dpd_->buf4_scm(&tau_temp, 0.0);
      half_transform(&tau2_AO_ab, &tau_temp, avir_c_, bvir_c_, navirpi_, nbvirpi_,
                     pq_row_start, Cd_row_start, false, 1.0, 0.0);
      global_dpd_->buf4_close(&tau2_AO_ab);
      global_dpd_->buf4_close(&tau_temp);

      // Transform the GTau contribution from MO to SO basis for the future use in the Fock operator
      // Alpha occupied
      global_dpd_->file2_init(&tau, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "GTau <O|O>");
      global_dpd_->file2_scm(&tau, 0.0);
      file2_transform(&s_aa_2, &tau, aocc_c_, false);
      global_dpd_->file2_close(&tau);

      // Alpha virtual
      global_dpd_->file2_init(&tau, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "GTau <V|V>");
      global_dpd_->file2_scm(&tau, 0.0);
      file2_transform(&s_aa_2, &tau, avir_c_, false);
      global_dpd_->file2_close(&tau);

      // Beta occupied
      global_dpd_->file2_init(&tau, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "GTau <o|o>");
      global_dpd_->file2_scm(&tau, 0.0);
      file2_transform(&s_bb_2, &tau, bocc_c_, false);
      global_dpd_->file2_close(&tau);

      // Beta virtual
      global_dpd_->file2_init(&tau, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "GTau <v|v>");
      global_dpd_->file2_scm(&tau, 0.0);
      file2_transform(&s_bb_2, &tau, bvir_c_, false);
      global_dpd_->file2_close(&tau);

      global_dpd_->file2_close(&s_aa_1);
      global_dpd_->file2_close(&s_aa_2);
      global_dpd_->file2_close(&s_bb_1);
      global_dpd_->file2_close(&s_bb_2);

      free_int_matrix(pq_row_start);
      free_int_matrix(CD_row_start);
      free_int_matrix(cd_row_start);
      free_int_matrix(Cd_row_start);

      dcft_timer_off("DCFTSolver::build_tensors");
  }



  /**
   * Builds the G matrix for Fock matrix, the external potential (tau contracted with the integrals)
   * and other tensors, if requested, out-of-core using the SO integrals.
   * All quantities are built simultaneously to reduce I/O.
   */
  void
  DCFTSolver::build_G()
  {
      dcft_timer_on("DCFTSolver::build_G");

      IWL *iwl = new IWL(psio_.get(), PSIF_SO_TEI, int_tolerance_, 1, 1);

      Label *lblptr = iwl->labels();
      Value *valptr = iwl->values();

      double *Da = init_array(ntriso_);
      double *Db = init_array(ntriso_);
      double *Ta = init_array(ntriso_);
      double *Tb = init_array(ntriso_);
      double *Ga = init_array(ntriso_);
      double *Gb = init_array(ntriso_);
      double *Va = init_array(ntriso_);
      double *Vb = init_array(ntriso_);
      int soOffset = 0;
      for(int h = 0; h < nirrep_; ++h){
          for(int mu = 0; mu < nsopi_[h]; ++ mu){
              for(int nu = 0; nu <= mu; ++ nu){
                  int muNu = INDEX((nu+soOffset), (mu+soOffset));
                  Da[muNu] = kappa_so_a_->get(h, mu, nu);
                  Db[muNu] = kappa_so_b_->get(h, mu, nu);
                  Ta[muNu] = tau_so_a_->get(h,mu,nu);
                  Tb[muNu] = tau_so_b_->get(h,mu,nu);
              }
          }
          soOffset += nsopi_[h];
      }

      double value;
      int Gpr, Grp, Gps, Gsp, Gsq, Gqs, Gqr, Grq, Gc, Gd;
      int pr, rp, ps, sp, qr, rq, qs, sq;
      int pqArr, qpArr, rsArr, srArr, qrArr, rqArr;
      int qsArr, sqArr, psArr, spArr, prArr, rpArr;
      int offset, labelIndex, p, q, r, s, h;


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

              qpArr = pqArr = INDEX(p, q);
              srArr = rsArr = INDEX(r, s);
              prArr = rpArr = INDEX(p, r);
              qsArr = sqArr = INDEX(q, s);
              spArr = psArr = INDEX(p, s);
              qrArr = rqArr = INDEX(q, r);

              /* (pq|rs) */
              Ga[rsArr] += (Da[pqArr] + Db[pqArr]) * value;
              Gb[rsArr] += (Da[pqArr] + Db[pqArr]) * value;
              Va[rsArr] += (Ta[pqArr] + Tb[pqArr]) * value;
              Vb[rsArr] += (Ta[pqArr] + Tb[pqArr]) * value;
              if(q >= r){
                  Ga[qrArr] -= Da[psArr] * value;
                  Gb[qrArr] -= Db[psArr] * value;
                  Va[qrArr] -= Ta[psArr] * value;
                  Vb[qrArr] -= Tb[psArr] * value;
              }

              if(p!=q && r!=s && pqArr!=rsArr){
                  /* (pq|sr) */
                  if(s >= r){
                      Ga[srArr] += (Da[pqArr] + Db[pqArr]) * value;
                      Gb[srArr] += (Da[pqArr] + Db[pqArr]) * value;
                      Va[srArr] += (Ta[pqArr] + Tb[pqArr]) * value;
                      Vb[srArr] += (Ta[pqArr] + Tb[pqArr]) * value;
                  }
                  if(q >= s){
                      Ga[qsArr] -= Da[prArr] * value;
                      Gb[qsArr] -= Db[prArr] * value;
                      Va[qsArr] -= Ta[prArr] * value;
                      Vb[qsArr] -= Tb[prArr] * value;
                  }

                  /* (qp|rs) */
                  if(r >= s){
                      Ga[rsArr] += (Da[qpArr] + Db[qpArr]) * value;
                      Gb[rsArr] += (Da[qpArr] + Db[qpArr]) * value;
                      Va[rsArr] += (Ta[qpArr] + Tb[qpArr]) * value;
                      Vb[rsArr] += (Ta[qpArr] + Tb[qpArr]) * value;
                  }
                  if(p >= r){
                      Ga[prArr] -= Da[qsArr] * value;
                      Gb[prArr] -= Db[qsArr] * value;
                      Va[prArr] -= Ta[qsArr] * value;
                      Vb[prArr] -= Tb[qsArr] * value;
                  }

                  /* (qp|sr) */
                  if(s >= r){
                      Ga[srArr] += (Da[qpArr] + Db[qpArr]) * value;
                      Gb[srArr] += (Da[qpArr] + Db[qpArr]) * value;
                      Va[srArr] += (Ta[qpArr] + Tb[qpArr]) * value;
                      Vb[srArr] += (Ta[qpArr] + Tb[qpArr]) * value;
                  }
                  if(p >= s){
                      Ga[psArr] -= Da[qrArr] * value;
                      Gb[psArr] -= Db[qrArr] * value;
                      Va[psArr] -= Ta[qrArr] * value;
                      Vb[psArr] -= Tb[qrArr] * value;
                  }

                  /* (rs|pq) */
                  if(p >= q){
                      Ga[pqArr] += (Da[rsArr] + Db[rsArr]) * value;
                      Gb[pqArr] += (Da[rsArr] + Db[rsArr]) * value;
                      Va[pqArr] += (Ta[rsArr] + Tb[rsArr]) * value;
                      Vb[pqArr] += (Ta[rsArr] + Tb[rsArr]) * value;
                  }
                  if(s >= p){
                      Ga[spArr] -= Da[rqArr] * value;
                      Gb[spArr] -= Db[rqArr] * value;
                      Va[spArr] -= Ta[rqArr] * value;
                      Vb[spArr] -= Tb[rqArr] * value;
                  }

                  /* (sr|pq) */
                  if(p >= q){
                      Ga[pqArr] += (Da[srArr] + Db[srArr]) * value;
                      Gb[pqArr] += (Da[srArr] + Db[srArr]) * value;
                      Va[pqArr] += (Ta[srArr] + Tb[srArr]) * value;
                      Vb[pqArr] += (Ta[srArr] + Tb[srArr]) * value;
                  }
                  if(r >= p){
                      Ga[rpArr] -= Da[sqArr] * value;
                      Gb[rpArr] -= Db[sqArr] * value;
                      Va[rpArr] -= Ta[sqArr] * value;
                      Vb[rpArr] -= Tb[sqArr] * value;
                  }

                  /* (rs|qp) */
                  if(q >= p){
                      Ga[qpArr] += (Da[rsArr] + Db[rsArr]) * value;
                      Gb[qpArr] += (Da[rsArr] + Db[rsArr]) * value;
                      Va[qpArr] += (Ta[rsArr] + Tb[rsArr]) * value;
                      Vb[qpArr] += (Ta[rsArr] + Tb[rsArr]) * value;
                  }
                  if(s >= q){
                      Ga[sqArr] -= Da[rpArr] * value;
                      Gb[sqArr] -= Db[rpArr] * value;
                      Va[sqArr] -= Ta[rpArr] * value;
                      Vb[sqArr] -= Tb[rpArr] * value;
                  }

                  /* (sr|qp) */
                  if(q >= p){
                      Ga[qpArr] += (Da[srArr] + Db[srArr]) * value;
                      Gb[qpArr] += (Da[srArr] + Db[srArr]) * value;
                      Va[qpArr] += (Ta[srArr] + Tb[srArr]) * value;
                      Vb[qpArr] += (Ta[srArr] + Tb[srArr]) * value;
                  }
                  if(r >= q){
                      Ga[rqArr] -= Da[spArr] * value;
                      Gb[rqArr] -= Db[spArr] * value;
                      Va[rqArr] -= Ta[spArr] * value;
                      Vb[rqArr] -= Tb[spArr] * value;
                  }
              }else if(p!=q && r!=s && pqArr==rsArr){
                  /* (pq|sr) */
                  if(s >= r){
                      Ga[srArr] += (Da[pqArr] + Db[pqArr]) * value;
                      Gb[srArr] += (Da[pqArr] + Db[pqArr]) * value;
                      Va[srArr] += (Ta[pqArr] + Tb[pqArr]) * value;
                      Vb[srArr] += (Ta[pqArr] + Tb[pqArr]) * value;
                  }
                  if(q >= s){
                      Ga[qsArr] -= Da[prArr] * value;
                      Gb[qsArr] -= Db[prArr] * value;
                      Va[qsArr] -= Ta[prArr] * value;
                      Vb[qsArr] -= Tb[prArr] * value;
                  }
                  /* (qp|rs) */
                  if(r >= s){
                      Ga[rsArr] += (Da[qpArr] + Db[qpArr]) * value;
                      Gb[rsArr] += (Da[qpArr] + Db[qpArr]) * value;
                      Va[rsArr] += (Ta[qpArr] + Tb[qpArr]) * value;
                      Vb[rsArr] += (Ta[qpArr] + Tb[qpArr]) * value;
                  }
                  if(p >= r){
                      Ga[prArr] -= Da[qsArr] * value;
                      Gb[prArr] -= Db[qsArr] * value;
                      Va[prArr] -= Ta[qsArr] * value;
                      Vb[prArr] -= Tb[qsArr] * value;
                  }

                  /* (qp|sr) */
                  if(s >= r){
                      Ga[srArr] += (Da[qpArr] + Db[qpArr]) * value;
                      Gb[srArr] += (Da[qpArr] + Db[qpArr]) * value;
                      Va[srArr] += (Ta[qpArr] + Tb[qpArr]) * value;
                      Vb[srArr] += (Ta[qpArr] + Tb[qpArr]) * value;
                  }
                  if(p >= s){
                      Ga[psArr] -= Da[qrArr] * value;
                      Gb[psArr] -= Db[qrArr] * value;
                      Va[psArr] -= Ta[qrArr] * value;
                      Vb[psArr] -= Tb[qrArr] * value;
                  }
              }else if(p!=q && r==s){
                  /* (qp|rs) */
                  if(r >= s){
                      Ga[rsArr] += (Da[qpArr] + Db[qpArr]) * value;
                      Gb[rsArr] += (Da[qpArr] + Db[qpArr]) * value;
                      Va[rsArr] += (Ta[qpArr] + Tb[qpArr]) * value;
                      Vb[rsArr] += (Ta[qpArr] + Tb[qpArr]) * value;
                  }
                  if(p >= r){
                      Ga[prArr] -= Da[qsArr] * value;
                      Gb[prArr] -= Db[qsArr] * value;
                      Va[prArr] -= Ta[qsArr] * value;
                      Vb[prArr] -= Tb[qsArr] * value;
                  }

                  /* (rs|pq) */
                  if(p >= q){
                      Ga[pqArr] += (Da[rsArr] + Db[rsArr]) * value;
                      Gb[pqArr] += (Da[rsArr] + Db[rsArr]) * value;
                      Va[pqArr] += (Ta[rsArr] + Tb[rsArr]) * value;
                      Vb[pqArr] += (Ta[rsArr] + Tb[rsArr]) * value;
                  }
                  if(s >= p){
                      Ga[spArr] -= Da[rqArr] * value;
                      Gb[spArr] -= Db[rqArr] * value;
                      Va[spArr] -= Ta[rqArr] * value;
                      Vb[spArr] -= Tb[rqArr] * value;
                  }

                  /* (rs|qp) */
                  if(q >= p){
                      Ga[qpArr] += (Da[rsArr] + Db[rsArr]) * value;
                      Gb[qpArr] += (Da[rsArr] + Db[rsArr]) * value;
                      Va[qpArr] += (Ta[rsArr] + Tb[rsArr]) * value;
                      Vb[qpArr] += (Ta[rsArr] + Tb[rsArr]) * value;
                  }
                  if(s >= q){
                      Ga[sqArr] -= Da[rpArr] * value;
                      Gb[sqArr] -= Db[rpArr] * value;
                      Va[sqArr] -= Ta[rpArr] * value;
                      Vb[sqArr] -= Tb[rpArr] * value;
                  }
              }else if(p==q && r!=s){
                  /* (pq|sr) */
                  if(s >= r){
                      Ga[srArr] += (Da[pqArr] + Db[pqArr]) * value;
                      Gb[srArr] += (Da[pqArr] + Db[pqArr]) * value;
                      Va[srArr] += (Ta[pqArr] + Tb[pqArr]) * value;
                      Vb[srArr] += (Ta[pqArr] + Tb[pqArr]) * value;
                  }
                  if(q >= s){
                      Ga[qsArr] -= Da[prArr] * value;
                      Gb[qsArr] -= Db[prArr] * value;
                      Va[qsArr] -= Ta[prArr] * value;
                      Vb[qsArr] -= Tb[prArr] * value;
                  }

                  /* (rs|pq) */
                  if(p >= q){
                      Ga[pqArr] += (Da[rsArr] + Db[rsArr]) * value;
                      Gb[pqArr] += (Da[rsArr] + Db[rsArr]) * value;
                      Va[pqArr] += (Ta[rsArr] + Tb[rsArr]) * value;
                      Vb[pqArr] += (Ta[rsArr] + Tb[rsArr]) * value;
                  }
                  if(s >= p){
                      Ga[spArr] -= Da[rqArr] * value;
                      Gb[spArr] -= Db[rqArr] * value;
                      Va[spArr] -= Ta[rqArr] * value;
                      Vb[spArr] -= Tb[rqArr] * value;
                  }

                  /* (sr|pq) */
                  if(p >= q){
                      Ga[pqArr] += (Da[srArr] + Db[srArr]) * value;
                      Gb[pqArr] += (Da[srArr] + Db[srArr]) * value;
                      Va[pqArr] += (Ta[srArr] + Tb[srArr]) * value;
                      Vb[pqArr] += (Ta[srArr] + Tb[srArr]) * value;
                  }
                  if(r >= p){
                      Ga[rpArr] -= Da[sqArr] * value;
                      Gb[rpArr] -= Db[sqArr] * value;
                      Va[rpArr] -= Ta[sqArr] * value;
                      Vb[rpArr] -= Tb[sqArr] * value;
                  }
              }else if(p==q && r==s && pqArr!=rsArr){
                  /* (rs|pq) */
                  if(p >= q){
                      Ga[pqArr] += (Da[rsArr] + Db[rsArr]) * value;
                      Gb[pqArr] += (Da[rsArr] + Db[rsArr]) * value;
                      Va[pqArr] += (Ta[rsArr] + Tb[rsArr]) * value;
                      Vb[pqArr] += (Ta[rsArr] + Tb[rsArr]) * value;
                  }
                  if(s >= p){
                      Ga[spArr] -= Da[rqArr] * value;
                      Gb[spArr] -= Db[rqArr] * value;
                      Va[spArr] -= Ta[rqArr] * value;
                      Vb[spArr] -= Tb[rqArr] * value;
                  }
              }
          } /* end loop through current buffer */
          if(!lastBuffer) iwl->fetch();
      }while(!lastBuffer);
      iwl->set_keep_flag(1);
      delete iwl;

      // Build the Fock matrices from the H and G matrices
      soOffset = 0;
      for(int h = 0; h < nirrep_; ++h){
          for(int mu = 0; mu < nsopi_[h]; ++mu){
              for(int nu = 0; nu <= mu; ++nu){
                  int muNu = INDEX((nu+soOffset), (mu+soOffset));
                  double aVal   = Ga[muNu];
                  double bVal   = Gb[muNu];
                  double aGTVal = Va[muNu];
                  double bGTVal = Vb[muNu];
                  Fa_->add(h, mu, nu, aVal);
                  Fb_->add(h, mu, nu, bVal);
                  g_tau_a_->set(h, mu, nu, aGTVal);
                  g_tau_b_->set(h, mu, nu, bGTVal);
                  if(mu != nu){
                      Fa_->add(h, nu, mu, aVal);
                      Fb_->add(h, nu, mu, bVal);
                      g_tau_a_->set(h, nu, mu, aGTVal);
                      g_tau_b_->set(h, nu, mu, bGTVal);
                  }
              }
          }
          soOffset += nsopi_[h];
      }

      free(Ta);
      free(Tb);
      free(Va);
      free(Vb);
      free(Da);
      free(Db);
      free(Ga);
      free(Gb);

      dcft_timer_off("DCFTSolver::build_G");
  }

}} // Namespaces
