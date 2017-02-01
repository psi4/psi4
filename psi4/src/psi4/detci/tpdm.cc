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

/*! \file
**  \ingroup DETCI
**  \brief Compute the two-particle density matrix (TPDM)
**
**  C. David Sherrill
*/

#include <cstdio>
#include <cstdlib>
#include <cmath>
/* may no longer need #include <libc.h> */
#include "psi4/psifiles.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libqt/qt.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"
#include "psi4/detci/structs.h"
#include "psi4/detci/civect.h"
#include "psi4/detci/ciwave.h"

namespace psi { namespace detci {

#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))

// DGAS this is still awkward, I think the TPDM code can be less general than the OPDM one for now.
void CIWavefunction::form_tpdm(void) {
  SharedCIVector Ivec = new_civector(Parameters_->num_roots, Parameters_->d_filenum);
  Ivec->init_io_files(true);
  SharedCIVector Jvec = new_civector(Parameters_->num_roots, Parameters_->d_filenum);
  Jvec->init_io_files(true);

  std::vector<SharedMatrix> tpdm_list;
  std::vector<std::tuple<int, int, double> > states_vec;
  for (int root_idx = 0; root_idx < Parameters_->average_num; root_idx++){
    states_vec.push_back(std::make_tuple(Parameters_->average_states[root_idx],
                                         Parameters_->average_states[root_idx],
                                         Parameters_->average_weights[root_idx]));
  }

  tpdm_list = tpdm(Ivec, Jvec, states_vec);

  Ivec->close_io_files(true); // Closes Jvec too

  tpdm_aa_ = tpdm_list[0];
  tpdm_ab_ = tpdm_list[1];
  tpdm_bb_ = tpdm_list[2];
  tpdm_    = tpdm_list[3];
  tpdm_called_ = true;
}

/*
Compute a single root without weighting
 */
std::vector<SharedMatrix> CIWavefunction::tpdm(SharedCIVector Ivec,
                                               SharedCIVector Jvec, int Iroot,
                                               int Jroot) {
  std::vector<std::tuple<int, int, double> > states_vec;
  states_vec.push_back(std::make_tuple(Iroot, Jroot, 1.0));
  return tpdm(Ivec, Jvec, states_vec);
}

/*
** Computes the two-particle density matrix between CIVectors I and J.
** The states vec tuple represents (Iroot, Jroot, weight).
** \Gamma_{tuvw} = < Iroot | \hat{e}^{tuvw} | Jroot >
*/
std::vector<SharedMatrix> CIWavefunction::tpdm(SharedCIVector Ivec, SharedCIVector Jvec,
                                               std::vector<std::tuple<int, int, double> > states_vec) {
  double **transp_tmp = nullptr;
  double **transp_tmp2 = nullptr;
  int Iblock, Iblock2, Ibuf, Iac, Ibc, Inas, Inbs, Iairr;
  int Jblock, Jblock2, Jbuf, Jac, Jbc, Jnas, Jnbs, Jairr;
  int do_Jblock, do_Jblock2;

  timer_on("CIWave: tpdm");
  if (!CalcInfo_->sigma_initialized) sigma_init(*(Ivec).get(), *(Jvec).get());

  int nact = CalcInfo_->num_ci_orbs;
  int nact2 = nact * nact;
  int ntri2 = (nact2 * (nact2 + 1)) / 2;

  SharedVector tpdm_aa(new Vector("MO-basis TPDM AA (ci order)", ntri2));
  SharedVector tpdm_ab(
      new Vector("MO-basis TPDM AB (ci order)", nact2 * nact2));
  SharedVector tpdm_bb(new Vector("MO-basis TPDM BB (ci order)", ntri2));
  double *tpdm_aap = tpdm_aa->pointer();
  double *tpdm_abp = tpdm_ab->pointer();
  double *tpdm_bbp = tpdm_bb->pointer();

  if ((Ivec->icore_ == 2 && Ivec->Ms0_ && CalcInfo_->ref_sym != 0) ||
      (Ivec->icore_ == 0 && Ivec->Ms0_)) {
    int maxrows = 0, maxcols = 0;
    for (int i = 0; i < Ivec->num_blocks_; i++) {
      if (Ivec->Ia_size_[i] > maxrows) maxrows = Ivec->Ia_size_[i];
      if (Ivec->Ib_size_[i] > maxcols) maxcols = Ivec->Ib_size_[i];
    }
    if (maxcols > maxrows) maxrows = maxcols;
    transp_tmp = (double **)malloc(maxrows * sizeof(double *));
    transp_tmp2 = (double **)malloc(maxrows * sizeof(double *));
    if (transp_tmp == nullptr || transp_tmp2 == nullptr) {
     outfile->Printf("(tpdm): Trouble with malloc'ing transp_tmp\n");
    }
    unsigned long bufsz = Ivec->get_max_blk_size();
    transp_tmp[0] = init_array(bufsz);
    transp_tmp2[0] = init_array(bufsz);
    if (transp_tmp[0] == nullptr || transp_tmp2[0] == nullptr) {
     outfile->Printf("(tpdm): Trouble with malloc'ing transp_tmp[0]\n");
    }
  }

  if (Parameters_->icore == 0) {
    /* loop over all the roots requested */
    for (int root_idx = 0; root_idx < states_vec.size(); root_idx++) {
      int Iroot = std::get<0>(states_vec[root_idx]);
      int Jroot = std::get<1>(states_vec[root_idx]);
      double weight = std::get<2>(states_vec[root_idx]);

      for (Ibuf = 0; Ibuf < Ivec->buf_per_vect_; Ibuf++) {
        Ivec->read(Iroot, Ibuf);
        Iblock = Ivec->buf2blk_[Ibuf];
        Iac = Ivec->Ia_code_[Iblock];
        Ibc = Ivec->Ib_code_[Iblock];
        Inas = Ivec->Ia_size_[Iblock];
        Inbs = Ivec->Ib_size_[Iblock];

        for (Jbuf = 0; Jbuf < Jvec->buf_per_vect_; Jbuf++) {
          do_Jblock = 0;
          do_Jblock2 = 0;
          Jblock = Jvec->buf2blk_[Jbuf];
          Jblock2 = -1;
          Jac = Jvec->Ia_code_[Jblock];
          Jbc = Jvec->Ib_code_[Jblock];
          if (Jvec->Ms0_) Jblock2 = Jvec->decode_[Jbc][Jac];
          Jnas = Jvec->Ia_size_[Jblock];
          Jnbs = Jvec->Ib_size_[Jblock];
          if (s1_contrib_[Iblock][Jblock] || s2_contrib_[Iblock][Jblock] ||
              s3_contrib_[Iblock][Jblock])
            do_Jblock = 1;
          if (Jvec->buf_offdiag_[Jbuf] &&
              (s1_contrib_[Iblock][Jblock2] || s2_contrib_[Iblock][Jblock2] ||
               s3_contrib_[Iblock][Jblock2]))
            do_Jblock2 = 1;
          if (!do_Jblock && !do_Jblock2) continue;

          Jvec->read(Jroot, Jbuf);

          if (do_Jblock) {
            tpdm_block(alplist_, betlist_, CalcInfo_->num_ci_orbs,
                       Ivec->num_alpcodes_, Ivec->num_betcodes_, tpdm_aap,
                       tpdm_bbp, tpdm_abp, Jvec->blocks_[Jblock],
                       Ivec->blocks_[Iblock], Jac, Jbc, Jnas, Jnbs, Iac, Ibc,
                       Inas, Inbs, weight);
          }

          if (do_Jblock2) {
            Jvec->transp_block(Jblock, transp_tmp);
            tpdm_block(alplist_, betlist_, CalcInfo_->num_ci_orbs,
                       Ivec->num_alpcodes_, Ivec->num_betcodes_, tpdm_aap,
                       tpdm_bbp, tpdm_abp, transp_tmp, Ivec->blocks_[Iblock],
                       Jbc, Jac, Jnbs, Jnas, Iac, Ibc, Inas, Inbs, weight);
          }

        } /* end loop over Jbuf */

        if (Ivec->buf_offdiag_[Ibuf]) { /* need to get contrib of transpose */
          Iblock2 = Ivec->decode_[Ibc][Iac];
          Iac = Ivec->Ia_code_[Iblock2];
          Ibc = Ivec->Ib_code_[Iblock2];
          Inas = Ivec->Ia_size_[Iblock2];
          Inbs = Ivec->Ib_size_[Iblock2];

          Ivec->transp_block(Iblock, transp_tmp2);

          for (Jbuf = 0; Jbuf < Jvec->buf_per_vect_; Jbuf++) {
            do_Jblock = 0;
            do_Jblock2 = 0;
            Jblock = Jvec->buf2blk_[Jbuf];
            Jblock2 = -1;
            Jac = Jvec->Ia_code_[Jblock];
            Jbc = Jvec->Ib_code_[Jblock];
            if (Jvec->Ms0_) Jblock2 = Jvec->decode_[Jbc][Jac];
            Jnas = Jvec->Ia_size_[Jblock];
            Jnbs = Jvec->Ib_size_[Jblock];
            if (s1_contrib_[Iblock2][Jblock] || s2_contrib_[Iblock2][Jblock] ||
                s3_contrib_[Iblock2][Jblock])
              do_Jblock = 1;
            if (Jvec->buf_offdiag_[Jbuf] && (s1_contrib_[Iblock2][Jblock2] ||
                                             s2_contrib_[Iblock2][Jblock2] ||
                                             s3_contrib_[Iblock2][Jblock2]))
              do_Jblock2 = 1;
            if (!do_Jblock && !do_Jblock2) continue;

            Jvec->read(Jroot, Jbuf);

            if (do_Jblock) {
              tpdm_block(alplist_, betlist_, CalcInfo_->num_ci_orbs,
                         Ivec->num_alpcodes_, Ivec->num_betcodes_, tpdm_aap,
                         tpdm_bbp, tpdm_abp, Jvec->blocks_[Jblock], transp_tmp2,
                         Jac, Jbc, Jnas, Jnbs, Iac, Ibc, Inas, Inbs, weight);
            }

            if (do_Jblock2) {
              Jvec->transp_block(Jblock, transp_tmp);
              tpdm_block(alplist_, betlist_, CalcInfo_->num_ci_orbs,
                         Ivec->num_alpcodes_, Ivec->num_betcodes_, tpdm_aap,
                         tpdm_bbp, tpdm_abp, transp_tmp, transp_tmp2, Jbc, Jac,
                         Jnbs, Jnas, Iac, Ibc, Inas, Inbs, weight);
            }
          } /* end loop over Jbuf */

        } /* end loop over Ibuf transpose */
      }   /* end loop over Ibuf */
    }     /* end loop over roots */
  }       /* end icore==0 */

  else if (Parameters_->icore == 1) { /* whole vectors in-core */
    for (int root_idx = 0; root_idx < states_vec.size(); root_idx++) {
      int Iroot = std::get<0>(states_vec[root_idx]);
      int Jroot = std::get<1>(states_vec[root_idx]);
      double weight = std::get<2>(states_vec[root_idx]);

      Ivec->read(Iroot, 0);
      Jvec->read(Jroot, 0);
      for (Iblock = 0; Iblock < Ivec->num_blocks_; Iblock++) {
        Iac = Ivec->Ia_code_[Iblock];
        Ibc = Ivec->Ib_code_[Iblock];
        Inas = Ivec->Ia_size_[Iblock];
        Inbs = Ivec->Ib_size_[Iblock];
        if (Inas == 0 || Inbs == 0) continue;
        for (Jblock = 0; Jblock < Jvec->num_blocks_; Jblock++) {
          Jac = Jvec->Ia_code_[Jblock];
          Jbc = Jvec->Ib_code_[Jblock];
          Jnas = Jvec->Ia_size_[Jblock];
          Jnbs = Jvec->Ib_size_[Jblock];
          if (s1_contrib_[Iblock][Jblock] || s2_contrib_[Iblock][Jblock] || s3_contrib_[Iblock][Jblock])
            tpdm_block(alplist_, betlist_, CalcInfo_->num_ci_orbs,
                       Ivec->num_alpcodes_, Ivec->num_betcodes_, tpdm_aap,
                       tpdm_bbp, tpdm_abp, Jvec->blocks_[Jblock],
                       Ivec->blocks_[Iblock], Jac, Jbc, Jnas, Jnbs, Iac, Ibc,
                       Inas, Inbs, weight);
        }
      } /* end loop over Iblock */
    }   /* end loop over roots */
  }     /* end icore==1 */

  else if (Parameters_->icore == 2) { /* icore==2 */
    for (int root_idx = 0; root_idx < states_vec.size(); root_idx++) {
      int Iroot = std::get<0>(states_vec[root_idx]);
      int Jroot = std::get<1>(states_vec[root_idx]);
      double weight = std::get<2>(states_vec[root_idx]);

      for (Ibuf = 0; Ibuf < Ivec->buf_per_vect_; Ibuf++) {
        Ivec->read(Iroot, Ibuf);
        Iairr = Ivec->buf2blk_[Ibuf];

        for (Jbuf = 0; Jbuf < Jvec->buf_per_vect_; Jbuf++) {
          Jvec->read(Jroot, Jbuf);
          Jairr = Jvec->buf2blk_[Jbuf];

          for (Iblock = Ivec->first_ablk_[Iairr];
               Iblock <= Ivec->last_ablk_[Iairr]; Iblock++) {
            Iac = Ivec->Ia_code_[Iblock];
            Ibc = Ivec->Ib_code_[Iblock];
            Inas = Ivec->Ia_size_[Iblock];
            Inbs = Ivec->Ib_size_[Iblock];

            for (Jblock = Jvec->first_ablk_[Jairr];
                 Jblock <= Jvec->last_ablk_[Jairr]; Jblock++) {
              Jac = Jvec->Ia_code_[Jblock];
              Jbc = Jvec->Ib_code_[Jblock];
              Jnas = Jvec->Ia_size_[Jblock];
              Jnbs = Jvec->Ib_size_[Jblock];

              if (s1_contrib_[Iblock][Jblock] || s2_contrib_[Iblock][Jblock] ||
                  s3_contrib_[Iblock][Jblock])
                tpdm_block(alplist_, betlist_, CalcInfo_->num_ci_orbs,
                           Ivec->num_alpcodes_, Ivec->num_betcodes_, tpdm_aap,
                           tpdm_bbp, tpdm_abp, Jvec->blocks_[Jblock],
                           Ivec->blocks_[Iblock], Jac, Jbc, Jnas, Jnbs, Iac,
                           Ibc, Inas, Inbs, weight);

              if (Jvec->buf_offdiag_[Jbuf]) {
                Jblock2 = Jvec->decode_[Jbc][Jac];
                if (s1_contrib_[Iblock][Jblock2] ||
                    s2_contrib_[Iblock][Jblock2] ||
                    s3_contrib_[Iblock][Jblock2]) {
                  Jvec->transp_block(Jblock, transp_tmp);
                  tpdm_block(alplist_, betlist_, CalcInfo_->num_ci_orbs,
                             Ivec->num_alpcodes_, Ivec->num_betcodes_, tpdm_aap,
                             tpdm_bbp, tpdm_abp, transp_tmp,
                             Ivec->blocks_[Iblock], Jbc, Jac, Jnbs, Jnas, Iac,
                             Ibc, Inas, Inbs, weight);
                }
              }

            } /* end loop over Jblock */

            if (Ivec->buf_offdiag_[Ibuf]) {
              Iblock2 = Ivec->decode_[Ibc][Iac];
              Ivec->transp_block(Iblock, transp_tmp2);
              Iac = Ivec->Ia_code_[Iblock2];
              Ibc = Ivec->Ib_code_[Iblock2];
              Inas = Ivec->Ia_size_[Iblock2];
              Inbs = Ivec->Ib_size_[Iblock2];

              for (Jblock = Jvec->first_ablk_[Jairr];
                   Jblock <= Jvec->last_ablk_[Jairr]; Jblock++) {
                Jac = Jvec->Ia_code_[Jblock];
                Jbc = Jvec->Ib_code_[Jblock];
                Jnas = Jvec->Ia_size_[Jblock];
                Jnbs = Jvec->Ib_size_[Jblock];

                if (s1_contrib_[Iblock2][Jblock] ||
                    s2_contrib_[Iblock2][Jblock] ||
                    s3_contrib_[Iblock2][Jblock])
                  tpdm_block(alplist_, betlist_, CalcInfo_->num_ci_orbs,
                             Ivec->num_alpcodes_, Ivec->num_betcodes_, tpdm_aap,
                             tpdm_bbp, tpdm_abp, Jvec->blocks_[Jblock],
                             transp_tmp2, Jac, Jbc, Jnas, Jnbs, Iac, Ibc, Inas,
                             Inbs, weight);

                if (Jvec->buf_offdiag_[Jbuf]) {
                  Jblock2 = Jvec->decode_[Jbc][Jac];
                  if (s1_contrib_[Iblock][Jblock2] ||
                      s2_contrib_[Iblock][Jblock2] ||
                      s3_contrib_[Iblock][Jblock2]) {
                    Jvec->transp_block(Jblock, transp_tmp);
                    tpdm_block(alplist_, betlist_, CalcInfo_->num_ci_orbs,
                               Ivec->num_alpcodes_, Ivec->num_betcodes_,
                               tpdm_aap, tpdm_bbp, tpdm_abp, transp_tmp,
                               transp_tmp2, Jbc, Jac, Jnbs, Jnas, Iac, Ibc,
                               Inas, Inbs, weight);
                  }
                }

              } /* end loop over Jblock */
            }   /* end Ivec offdiag */

          } /* end loop over Iblock */
        }   /* end loop over Jbuf */
      }     /* end loop over Ibuf */
    }       /* end loop over roots */
  }         /* end icore==2 */

  else {
    throw PSIEXCEPTION("CIWavefunction::tpdm: unrecognized core option!\n");
  }

  // Symmetrize and reorder the TPDM
  SharedMatrix tpdm_aam(new Matrix("MO-basis TPDM AA", nact2, nact2));
  SharedMatrix tpdm_abm(new Matrix("MO-basis TPDM AB", nact2, nact2));
  SharedMatrix tpdm_bbm(new Matrix("MO-basis TPDM BB", nact2, nact2));
  double **tpdm_aamp = tpdm_aam->pointer();
  double **tpdm_abmp = tpdm_abm->pointer();
  double **tpdm_bbmp = tpdm_bbm->pointer();

  // Reorder our density matrices
  for (int p = 0; p < nact; p++) {
    for (int q = 0; q < nact; q++) {
      for (int r = 0; r < nact; r++) {
        for (int s = 0; s < nact; s++) {
          // Reorder index
          int r_p = CalcInfo_->act_order[p];
          int r_q = CalcInfo_->act_order[q];
          int r_r = CalcInfo_->act_order[r];
          int r_s = CalcInfo_->act_order[s];
          int r_pq = r_p * nact + r_q;
          int r_rs = r_r * nact + r_s;

          // aa/bb index
          int pq = p * nact + q;
          int rs = r * nact + s;
          size_t pqrs = INDEX(pq, rs);

          tpdm_aamp[r_pq][r_rs] = tpdm_aap[pqrs];
          tpdm_abmp[r_pq][r_rs] = tpdm_abp[pq * nact2 + rs];
          tpdm_bbmp[r_pq][r_rs] = tpdm_bbp[pqrs];
        }
      }
    }
  }
  tpdm_aa.reset();
  tpdm_ab.reset();
  tpdm_bb.reset();

  // Build our spin summed density matrix
  SharedMatrix tpdm(new Matrix("MO-basis TPDM", nact2, nact2));
  double **tpdmp = tpdm->pointer();

  for (int p = 0; p < nact; p++) {
    for (int q = 0; q < nact; q++) {
      for (int r = 0; r < nact; r++) {
        for (int s = 0; s < nact; s++) {
          int pq = p * nact + q;
          int qp = q * nact + p;
          int rs = r * nact + s;
          int sr = s * nact + r;

          tpdmp[pq][rs] = 0.5 * (tpdm_aamp[pq][rs] + tpdm_bbmp[pq][rs] +
                                 tpdm_abmp[pq][rs] + tpdm_abmp[rs][pq]);
        }
      }
    }
  }

  // Ivec->buf_unlock();
  // Jvec->buf_unlock();
  if (transp_tmp != nullptr) free_block(transp_tmp);
  if (transp_tmp2 != nullptr) free_block(transp_tmp2);

  std::vector<int> nshape{nact, nact, nact, nact};
  tpdm_aam->set_numpy_shape(nshape);
  tpdm_abm->set_numpy_shape(nshape);
  tpdm_bbm->set_numpy_shape(nshape);
  tpdm->set_numpy_shape(nshape);

  std::vector<SharedMatrix> ret_list;
  ret_list.push_back(tpdm_aam);
  ret_list.push_back(tpdm_abm);
  ret_list.push_back(tpdm_bbm);
  ret_list.push_back(tpdm);

  timer_off("CIWave: tpdm");

  return ret_list;
}

void CIWavefunction::tpdm_block(struct stringwr **alplist,
                                struct stringwr **betlist, int nbf,
                                int nalplists, int nbetlists, double *twopdm_aa,
                                double *twopdm_bb, double *twopdm_ab,
                                double **CJ, double **CI, int Ja_list,
                                int Jb_list, int Jnas, int Jnbs, int Ia_list,
                                int Ib_list, int Inas, int Inbs,
                                double weight) {
  const int nbf2 = nbf * nbf;
  int Ia_idx, Ib_idx, Ja_idx, Jb_idx, Ja_ex, Jb_ex, Jbcnt, Jacnt;
  int Kbcnt, Kacnt, Kb_ex, Ka_ex, Kb_list, Ka_list, Kb_idx, Ka_idx;
  struct stringwr *Jb, *Ja, *Kb, *Ka;
  signed char *Jbsgn, *Jasgn, *Kbsgn, *Kasgn;
  unsigned int *Jbridx, *Jaridx, *Kbridx, *Karidx;
  double C1, C2, Ib_sgn, Ia_sgn, Kb_sgn, Ka_sgn, tval;
  int i, j, k, l, ij, kl, ijkl, oij, okl, *Jboij, *Jaoij, *Kboij, *Kaoij;

  /* loop over Ia in Ia_list */
  if (Ia_list == Ja_list) {
    for (Ia_idx = 0; Ia_idx < Inas; Ia_idx++) {
      for (Jb = betlist[Jb_list], Jb_idx = 0; Jb_idx < Jnbs; Jb_idx++, Jb++) {
        C1 = CJ[Ia_idx][Jb_idx] * weight;

        /* loop over excitations E^b_{kl} from |B(J_b)> */
        for (Kb_list = 0; Kb_list < nbetlists; Kb_list++) {
          Jbcnt = Jb->cnt[Kb_list];
          Jbridx = Jb->ridx[Kb_list];
          Jbsgn = Jb->sgn[Kb_list];
          Jboij = Jb->oij[Kb_list];
          for (Jb_ex = 0; Jb_ex < Jbcnt; Jb_ex++) {
            okl = *Jboij++;
            Kb_idx = *Jbridx++;
            Kb_sgn = (double)*Jbsgn++;

            Kb = betlist[Kb_list] + Kb_idx;
            if (Kb_list == Ib_list) {
              C2 = CI[Ia_idx][Kb_idx];
              i = okl / nbf;
              l = okl % nbf;
              for (j = 0; j < nbf && j <= i; j++) {
                ij = i * nbf + j;
                kl = j * nbf + l;
                if (ij >= kl) {
                  ijkl = INDEX(ij, kl);
                  twopdm_bb[ijkl] -= Kb_sgn * C1 * C2;
                }
              }
            }

            /* loop over excitations E^b_{ij} from |B(K_b)> */
            /* Ib_list pre-determined because of C blocking */
            Kbcnt = Kb->cnt[Ib_list];
            Kbridx = Kb->ridx[Ib_list];
            Kbsgn = Kb->sgn[Ib_list];
            Kboij = Kb->oij[Ib_list];
            for (Kb_ex = 0; Kb_ex < Kbcnt; Kb_ex++) {
              Ib_idx = *Kbridx++;
              Ib_sgn = (double)*Kbsgn++;
              oij = *Kboij++;
              if (oij >= okl) {
                C2 = CI[Ia_idx][Ib_idx];
                ijkl = INDEX(oij, okl);
                twopdm_bb[ijkl] += Ib_sgn * Kb_sgn * C1 * C2;
              }
            }

          } /* end loop over Jb_ex */
        }   /* end loop over Kb_list */
      }     /* end loop over Jb_idx */
    }       /* end loop over Ia_idx */
  }         /* end case Ia_list == Ja_list */

  /* loop over Ib in Ib_list */
  if (Ib_list == Jb_list) {
    for (Ib_idx = 0; Ib_idx < Inbs; Ib_idx++) {
      for (Ja = alplist[Ja_list], Ja_idx = 0; Ja_idx < Jnas; Ja_idx++, Ja++) {
        C1 = CJ[Ja_idx][Ib_idx] * weight;

        /* loop over excitations E^a_{kl} from |A(J_a)> */
        for (Ka_list = 0; Ka_list < nalplists; Ka_list++) {
          Jacnt = Ja->cnt[Ka_list];
          Jaridx = Ja->ridx[Ka_list];
          Jasgn = Ja->sgn[Ka_list];
          Jaoij = Ja->oij[Ka_list];
          for (Ja_ex = 0; Ja_ex < Jacnt; Ja_ex++) {
            okl = *Jaoij++;
            Ka_idx = *Jaridx++;
            Ka_sgn = (double)*Jasgn++;

            Ka = alplist[Ka_list] + Ka_idx;
            if (Ka_list == Ia_list) {
              C2 = CI[Ka_idx][Ib_idx];
              i = okl / nbf;
              l = okl % nbf;
              for (j = 0; j < nbf && j <= i; j++) {
                ij = i * nbf + j;
                kl = j * nbf + l;
                if (ij >= kl) {
                  ijkl = INDEX(ij, kl);
                  twopdm_aa[ijkl] -= Ka_sgn * C1 * C2;
                }
              }
            }

            /* loop over excitations E^a_{ij} from |A(K_a)> */
            /* Ia_list pre-determined because of C blocking */
            Kacnt = Ka->cnt[Ia_list];
            Karidx = Ka->ridx[Ia_list];
            Kasgn = Ka->sgn[Ia_list];
            Kaoij = Ka->oij[Ia_list];
            for (Ka_ex = 0; Ka_ex < Kacnt; Ka_ex++) {
              Ia_idx = *Karidx++;
              Ia_sgn = (double)*Kasgn++;
              oij = *Kaoij++;
              if (oij >= okl) {
                C2 = CI[Ia_idx][Ib_idx];
                ijkl = INDEX(oij, okl);
                twopdm_aa[ijkl] += Ia_sgn * Ka_sgn * C1 * C2;
              }
            }

          } /* end loop over Ja_ex */
        }   /* end loop over Ka_list */
      }     /* end loop over Ja_idx */
    }       /* end loop over Ib_idx */
  }         /* end case Ib_list == Jb_list */

  /* now do the sigma3 looking (alpha-beta) part */
  /* loop over Ja                                */
  for (Ja = alplist[Ja_list], Ja_idx = 0; Ja_idx < Jnas; Ja_idx++, Ja++) {
    /* loop over excitations E^a_{kl} from |A(I_a)> */
    Jacnt = Ja->cnt[Ia_list];
    Jaridx = Ja->ridx[Ia_list];
    Jasgn = Ja->sgn[Ia_list];
    Jaoij = Ja->oij[Ia_list];
    for (Ja_ex = 0; Ja_ex < Jacnt; Ja_ex++) {
      okl = *Jaoij++;
      Ia_idx = *Jaridx++;
      Ia_sgn = (double)*Jasgn++;

      /* loop over Jb */
      for (Jb = betlist[Jb_list], Jb_idx = 0; Jb_idx < Jnbs; Jb_idx++, Jb++) {
        C1 = CJ[Ja_idx][Jb_idx] * weight;

        /* loop over excitations E^b_{ij} from |B(J_b)> */
        Jbcnt = Jb->cnt[Ib_list];
        Jbridx = Jb->ridx[Ib_list];
        Jbsgn = Jb->sgn[Ib_list];
        Jboij = Jb->oij[Ib_list];

        for (Jb_ex = 0; Jb_ex < Jbcnt; Jb_ex++) {
          oij = *Jboij++;
          Ib_idx = *Jbridx++;
          Ib_sgn = (double)*Jbsgn++;
          C2 = CI[Ia_idx][Ib_idx];
          // alpha-beta matrix is stored without packing bra and ket together
          ijkl = oij * nbf2 + okl;
          tval = Ib_sgn * Ia_sgn * C1 * C2;
          // in orbital (i.e. non-spi-orbital) code had to scale the diagonal by 2
          // because d(ij,kl) += d_ab(ij,kl) + d_ab(kl,ij), hence
          // d(ij,ij) += 2 d_ab(ij,ij)
          // if (oij == okl) tval *= 2.0;
          twopdm_ab[ijkl] += tval;
        }
      } /* end loop over Jb */
    }   /* end loop over Ja_ex */
  }     /* end loop over Ja */
}

}} // namespace psi::detci
