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
    \ingroup DETCI
    \brief Enter brief description of file here
*/

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libqt/qt.h"
#include "psi4/libmints/oeprop.h"
#include "psi4/psifiles.h"
#include "psi4/physconst.h"
#include "psi4/detci/structs.h"
#include "psi4/detci/civect.h"
#include "psi4/detci/ciwave.h"

namespace psi { namespace detci {

#define MIN0(a,b) (((a)<(b)) ? (a) : (b))
#define MAX0(a,b) (((a)>(b)) ? (a) : (b))
#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))
#define TOL 1E-14

void CIWavefunction::form_opdm(void) {
    // if we're trying to follow a root, figure out which one here
    // CDS help: Why is this here and where can we move it?
    if (Parameters_->follow_vec_num > 0) {
        CIvect Ivec(Parameters_->icore, Parameters_->num_roots, 1,
                    Parameters_->d_filenum, CIblks_, CalcInfo_, Parameters_,
                    H0block_, false);
        Ivec.init_io_files(true);
        double max_overlap = 0.0, overlap = 0.0;
        int j = 0;
        for (int i = 0; i < Parameters_->num_roots; i++) {
            overlap = Ivec.compute_follow_overlap(
                i, Parameters_->follow_vec_num, Parameters_->follow_vec_coef.data(),
                Parameters_->follow_vec_Iac.data(), Parameters_->follow_vec_Iaridx.data(),
                Parameters_->follow_vec_Ibc.data(), Parameters_->follow_vec_Ibridx.data());

            if (overlap > max_overlap) {
                max_overlap = overlap;
                j = i;
            }
        }

        Parameters_->root = j;
        outfile->Printf("DETCI following root %d (overlap %6.4lf)\n",
                        Parameters_->root + 1, max_overlap);

        for (int i = 0; i < Parameters_->average_num; i++) {
            Parameters_->average_states[i] = 0;
            Parameters_->average_weights[i] = 0.0;
        }
        Parameters_->average_num = 1;
        Parameters_->average_states[0] = Parameters_->root;
        Parameters_->average_weights[0] = 1.0;

        /* correct "the" energy in the checkpoint file */
        Process::environment.globals["CURRENT ENERGY"] = overlap;
        Process::environment.globals["CI TOTAL ENERGY"] = overlap;

        // eref is wrong for open-shells so replace it with escf until
        // I fix it, CDS 11/5/11
        Process::environment.globals["CI CORRELATION ENERGY"] = overlap - CalcInfo_->escf;
        Process::environment.globals["CURRENT CORRELATION ENERGY"] = overlap - CalcInfo_->escf;
        Process::environment.globals["CURRENT REFERENCE ENERGY"] = CalcInfo_->escf;
    }

    SharedCIVector Ivec = new_civector(Parameters_->num_roots, Parameters_->d_filenum);
    Ivec->init_io_files(true);
    SharedCIVector Jvec = new_civector(Parameters_->num_roots, Parameters_->d_filenum);
    Jvec->init_io_files(true);

    std::vector<std::vector<SharedMatrix> > opdm_list;
    std::vector<std::tuple<int, int> > states_vec;

    // Transition-OPDM's
    if (Parameters_->transdens) {
        states_vec.clear();
        for (int i = 0; i < Parameters_->num_roots - 1; ++i) {
            states_vec.push_back(std::make_tuple(0, i + 1));
        }
        opdm_list = opdm(Ivec, Jvec, states_vec);
        for (int i = 0; i < Parameters_->num_roots - 1; i++) {
            opdm_map_[opdm_list[i][0]->name()] = opdm_list[i][0];
            opdm_map_[opdm_list[i][1]->name()] = opdm_list[i][1];
            opdm_map_[opdm_list[i][2]->name()] = opdm_list[i][2];
        }
    }

    // OPDM's
    states_vec.clear();
    for (int i = 0; i < Parameters_->num_roots; ++i) {
        states_vec.push_back(std::make_tuple(i, i));
    }
    opdm_list = opdm(Ivec, Jvec, states_vec);
    for (int i = 0; i < Parameters_->num_roots; i++) {
        opdm_map_[opdm_list[i][0]->name()] = opdm_list[i][0];
        opdm_map_[opdm_list[i][1]->name()] = opdm_list[i][1];
        opdm_map_[opdm_list[i][2]->name()] = opdm_list[i][2];
    }
    Ivec->close_io_files(true); // Closes Jvec too

    // Figure out which OPDM should be current
    if (Parameters_->opdm_ave) {
        Dimension act_dim = get_dimension("ACT");
        opdm_a_ = SharedMatrix(new Matrix("MO-basis Alpha OPDM", nirrep_, act_dim, act_dim));
        opdm_b_ = SharedMatrix(new Matrix("MO-basis Beta OPDM", nirrep_, act_dim, act_dim));
        opdm_ = SharedMatrix(new Matrix("MO-basis OPDM", nirrep_, act_dim, act_dim));

        for (int i = 0; i < Parameters_->average_num; i++) {
            int croot = Parameters_->average_states[i];
            double weight = Parameters_->average_weights[i];
            opdm_a_->axpy(weight, opdm_list[croot][0]);
            opdm_b_->axpy(weight, opdm_list[croot][1]);
            opdm_->axpy(weight, opdm_list[croot][2]);
        }
    } else {
        int croot = Parameters_->root;
        opdm_a_ = opdm_list[croot][0]->clone();
        opdm_b_ = opdm_list[croot][1]->clone();
        opdm_ = opdm_list[croot][2]->clone();
    }

    SharedMatrix MO_Da = opdm_add_inactive(opdm_a_, 1.0, true);
    SharedMatrix MO_Db = opdm_add_inactive(opdm_b_, 1.0, true);
    Da_ = Matrix::triplet(Ca_, MO_Da, Ca_, false, false, true);
    Db_ = Matrix::triplet(Cb_, MO_Db, Cb_, false, false, true);

    opdm_called_ = true;
}
/*
** Resizes an act by act OPDM matrix to include the inactive porition.
** The diagonal of the inactive portion is set to value. If virt is true
** the return OPDM with of size nmo x nmo.
*/
SharedMatrix CIWavefunction::opdm_add_inactive(SharedMatrix opdm, double value,
                                               bool virt) {
    Dimension dim_drc = get_dimension("DRC");
    Dimension dim_act = get_dimension("ACT");
    Dimension dim_ia = dim_drc + dim_act;
    Dimension dim_ret;

    if (virt) {
        dim_ret = dim_ia + get_dimension("DRV");
    } else {
        dim_ret = dim_ia;
    }

    SharedMatrix ret(new Matrix(opdm->name(), dim_ret, dim_ret));

    for (int h = 0; h < nirrep_; h++) {
        if (!dim_ia[h]) continue;

        double **retp = ret->pointer(h);
        double **opdmp = opdm->pointer(h);

        // Diagonal inact part
        for (int i = 0; i < dim_drc[h]; i++) {
            retp[i][i] = value;
        }

        // Non-diagonal act part
        for (int t = 0; t < dim_act[h]; t++) {
            for (int u = 0; u < dim_act[h]; u++) {
                retp[dim_drc[h] + t][dim_drc[h] + u] = opdmp[t][u];
            }
        }
    }
    return ret;
}

std::vector<SharedMatrix> CIWavefunction::opdm(SharedCIVector Ivec,
                                               SharedCIVector Jvec, int Iroot,
                                               int Jroot) {
    std::vector<std::tuple<int, int> > states_vec;
    states_vec.push_back(std::make_tuple(Iroot, Jroot));
    return opdm(Ivec, Jvec, states_vec)[0];
}

/*
** Computes the one-particle density matrix for all nroots starting at root_start.
** If transden is true, opdm will then compute transition matrices for Iroot = root_start
** and Jroot from Iroot to Iroot + Nroots.
** \gamma_{tu} = < Iroot | \hat{E}^{tu} | Jroot >
*/
std::vector<std::vector<SharedMatrix> > CIWavefunction::opdm(SharedCIVector Ivec, SharedCIVector Jvec,
                                                            std::vector<std::tuple<int, int> > states_vec)
{
  double **transp_tmp = nullptr;
  double **transp_tmp2 = nullptr;
  int Iblock, Iblock2, Ibuf, Iac, Ibc, Inas, Inbs, Iairr;
  int Jblock, Jblock2, Jbuf, Jac, Jbc, Jnas, Jnbs, Jairr;
  int do_Jblock, do_Jblock2;

  timer_on("CIWave: opdm");
  if (!CalcInfo_->sigma_initialized) sigma_init(*(Ivec).get(), *(Jvec).get());

  std::vector<std::vector<SharedMatrix> > opdm_list;

  // Alloc trans_tmp arrays if needed
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
     outfile->Printf("(opdm): Trouble with malloc'ing transp_tmp\n");
    }
    unsigned long bufsz = Ivec->get_max_blk_size();
    transp_tmp[0] = init_array(bufsz);
    transp_tmp2[0] = init_array(bufsz);
    if (transp_tmp[0] == nullptr || transp_tmp2[0] == nullptr) {
     outfile->Printf("(opdm): Trouble with malloc'ing transp_tmp[0]\n");
    }
  }

  int nci = CalcInfo_->num_ci_orbs;
  Dimension act_dim = get_dimension("ACT");
  SharedMatrix scratch_a(new Matrix("OPDM A Scratch", nci, nci));
  SharedMatrix scratch_b(new Matrix("OPDM B Scratch", nci, nci));
  double** scratch_ap = scratch_a->pointer();
  double** scratch_bp = scratch_b->pointer();

  for (int root_idx=0; root_idx<states_vec.size(); root_idx++)  {
    int Iroot = std::get<0>(states_vec[root_idx]);
    int Jroot = std::get<1>(states_vec[root_idx]);

    scratch_a->zero();
    scratch_b->zero();

    if (Parameters_->icore == 0) {

      for (Ibuf=0; Ibuf<Ivec->buf_per_vect_; Ibuf++) {
        Ivec->read(Iroot, Ibuf);
        Iblock = Ivec->buf2blk_[Ibuf];
        Iac = Ivec->Ia_code_[Iblock];
        Ibc = Ivec->Ib_code_[Iblock];
        Inas = Ivec->Ia_size_[Iblock];
        Inbs = Ivec->Ib_size_[Iblock];

        for (Jbuf=0; Jbuf<Jvec->buf_per_vect_; Jbuf++) {
          do_Jblock=0; do_Jblock2=0;
          Jblock = Jvec->buf2blk_[Jbuf];
          Jblock2 = -1;
          Jac = Jvec->Ia_code_[Jblock];
          Jbc = Jvec->Ib_code_[Jblock];
          if (Jvec->Ms0_) Jblock2 = Jvec->decode_[Jbc][Jac];
          Jnas = Jvec->Ia_size_[Jblock];
          Jnbs = Jvec->Ib_size_[Jblock];
          if (s1_contrib_[Iblock][Jblock] || s2_contrib_[Iblock][Jblock])
            do_Jblock = 1;
          if (Jvec->buf_offdiag_[Jbuf] && (s1_contrib_[Iblock][Jblock2] ||
                                         s2_contrib_[Iblock][Jblock2]))
            do_Jblock2 = 1;
          if (!do_Jblock && !do_Jblock2) continue;

          Jvec->read(Jroot, Jbuf);

          if (do_Jblock) {
            opdm_block(alplist_, betlist_, scratch_ap, scratch_bp, Jvec->blocks_[Jblock],
                       Ivec->blocks_[Iblock], Jac, Jbc, Jnas,
                       Jnbs, Iac, Ibc, Inas, Inbs);
            }

          if (do_Jblock2) {
            Jvec->transp_block(Jblock, transp_tmp);
            opdm_block(alplist_, betlist_, scratch_ap, scratch_bp, transp_tmp,
                       Ivec->blocks_[Iblock], Jbc, Jac, Jnbs,
                       Jnas, Iac, Ibc, Inas, Inbs);
          }

        } /* end loop over Jbuf */

        if (Ivec->buf_offdiag_[Ibuf]) { /* need to get contrib of transpose */
          Iblock2 = Ivec->decode_[Ibc][Iac];
          Iac = Ivec->Ia_code_[Iblock2];
          Ibc = Ivec->Ib_code_[Iblock2];
          Inas = Ivec->Ia_size_[Iblock2];
          Inbs = Ivec->Ib_size_[Iblock2];

          Ivec->transp_block(Iblock, transp_tmp2);

          for (Jbuf=0; Jbuf<Jvec->buf_per_vect_; Jbuf++) {
            do_Jblock=0; do_Jblock2=0;
            Jblock = Jvec->buf2blk_[Jbuf];
            Jblock2 = -1;
            Jac = Jvec->Ia_code_[Jblock];
            Jbc = Jvec->Ib_code_[Jblock];
            if (Jvec->Ms0_) Jblock2 = Jvec->decode_[Jbc][Jac];
            Jnas = Jvec->Ia_size_[Jblock];
            Jnbs = Jvec->Ib_size_[Jblock];
            if (s1_contrib_[Iblock2][Jblock] || s2_contrib_[Iblock2][Jblock])
              do_Jblock = 1;
            if (Jvec->buf_offdiag_[Jbuf] && (s1_contrib_[Iblock2][Jblock2] ||
                                           s2_contrib_[Iblock2][Jblock2]))
              do_Jblock2 = 1;
            if (!do_Jblock && !do_Jblock2) continue;

            Jvec->read(Jroot, Jbuf);

            if (do_Jblock) {
              opdm_block(alplist_, betlist_, scratch_ap, scratch_bp, Jvec->blocks_[Jblock],
                         transp_tmp2, Jac, Jbc, Jnas,
                         Jnbs, Iac, Ibc, Inas, Inbs);
            }

            if (do_Jblock2) {
              Jvec->transp_block(Jblock, transp_tmp);
              opdm_block(alplist_, betlist_, scratch_ap, scratch_bp, transp_tmp,
                         transp_tmp2, Jbc, Jac, Jnbs,
                         Jnas, Iac, Ibc, Inas, Inbs);
            }
          } /* end loop over Jbuf */
        } /* end loop over Ibuf transpose */
      } /* end loop over Ibuf */
    } /* end icore==0 */

    else if (Parameters_->icore==1) { /* whole vectors in-core */
      Ivec->read(Iroot, 0);
      Jvec->read(Jroot, 0);
      for (Iblock=0; Iblock<Ivec->num_blocks_; Iblock++) {
        Iac = Ivec->Ia_code_[Iblock];
        Ibc = Ivec->Ib_code_[Iblock];
        Inas = Ivec->Ia_size_[Iblock];
        Inbs = Ivec->Ib_size_[Iblock];
        if (Inas==0 || Inbs==0) continue;
        for (Jblock=0; Jblock<Jvec->num_blocks_; Jblock++) {
          Jac = Jvec->Ia_code_[Jblock];
          Jbc = Jvec->Ib_code_[Jblock];
          Jnas = Jvec->Ia_size_[Jblock];
          Jnbs = Jvec->Ib_size_[Jblock];
          if (s1_contrib_[Iblock][Jblock] || s2_contrib_[Iblock][Jblock])
            opdm_block(alplist_, betlist_, scratch_ap, scratch_bp, Jvec->blocks_[Jblock],
                       Ivec->blocks_[Iblock], Jac, Jbc, Jnas,
                       Jnbs, Iac, Ibc, Inas, Inbs);
        }
      } /* end loop over Iblock */
    } /* end icore==1 */

    else if (Parameters_->icore==2) { /* icore==2 */
      for (Ibuf=0; Ibuf<Ivec->buf_per_vect_; Ibuf++) {
        Ivec->read(Iroot, Ibuf);
        Iairr = Ivec->buf2blk_[Ibuf];

        for (Jbuf=0; Jbuf<Jvec->buf_per_vect_; Jbuf++) {
          Jvec->read(Jroot, Jbuf);
          Jairr = Jvec->buf2blk_[Jbuf];

        for (Iblock=Ivec->first_ablk_[Iairr]; Iblock<=Ivec->last_ablk_[Iairr];
             Iblock++) {
          Iac = Ivec->Ia_code_[Iblock];
          Ibc = Ivec->Ib_code_[Iblock];
          Inas = Ivec->Ia_size_[Iblock];
          Inbs = Ivec->Ib_size_[Iblock];

          for (Jblock=Jvec->first_ablk_[Jairr]; Jblock<=Jvec->last_ablk_[Jairr];
               Jblock++) {
            Jac = Jvec->Ia_code_[Jblock];
            Jbc = Jvec->Ib_code_[Jblock];
            Jnas = Jvec->Ia_size_[Jblock];
            Jnbs = Jvec->Ib_size_[Jblock];

            if (s1_contrib_[Iblock][Jblock] || s2_contrib_[Iblock][Jblock])
              opdm_block(alplist_, betlist_, scratch_ap, scratch_bp, Jvec->blocks_[Jblock],
                         Ivec->blocks_[Iblock], Jac, Jbc, Jnas,
                         Jnbs, Iac, Ibc, Inas, Inbs);

            if (Jvec->buf_offdiag_[Jbuf]) {
              Jblock2 = Jvec->decode_[Jbc][Jac];
              if (s1_contrib_[Iblock][Jblock2] ||
                  s2_contrib_[Iblock][Jblock2]) {
              Jvec->transp_block(Jblock, transp_tmp);
                opdm_block(alplist_, betlist_, scratch_ap, scratch_bp, transp_tmp,
                  Ivec->blocks_[Iblock], Jbc, Jac,
                  Jnbs, Jnas, Iac, Ibc, Inas, Inbs);
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

            for (Jblock=Jvec->first_ablk_[Jairr]; Jblock<=Jvec->last_ablk_[Jairr];
              Jblock++) {
              Jac = Jvec->Ia_code_[Jblock];
              Jbc = Jvec->Ib_code_[Jblock];
              Jnas = Jvec->Ia_size_[Jblock];
              Jnbs = Jvec->Ib_size_[Jblock];

              if (s1_contrib_[Iblock2][Jblock] || s2_contrib_[Iblock2][Jblock])
                opdm_block(alplist_, betlist_, scratch_ap, scratch_bp, Jvec->blocks_[Jblock],
                           transp_tmp2, Jac, Jbc, Jnas, Jnbs, Iac, Ibc,
                           Inas, Inbs);

                if (Jvec->buf_offdiag_[Jbuf]) {
                  Jblock2 = Jvec->decode_[Jbc][Jac];
                  if (s1_contrib_[Iblock][Jblock2] ||
                    s2_contrib_[Iblock][Jblock2]) {
                    Jvec->transp_block(Jblock, transp_tmp);
                    opdm_block(alplist_, betlist_, scratch_ap, scratch_bp, transp_tmp,
                      transp_tmp2, Jbc, Jac, Jnbs, Jnas, Iac, Ibc, Inas, Inbs);
                  }
          }

        } /* end loop over Jblock */
            } /* end Ivec offdiag */

          } /* end loop over Iblock */
        } /* end loop over Jbuf */
      } /* end loop over Ibuf */
    } /* end icore==2 */

    else {
      throw PSIEXCEPTION("CIWavefunction::opdm: unrecognized core option!\n");
    }

    std::stringstream opdm_name;
    opdm_name << "MO-basis Alpha OPDM <" << Iroot << "| Etu |" << Jroot << ">";
    SharedMatrix new_OPDM_a(new Matrix(opdm_name.str(), nirrep_, act_dim, act_dim));

    opdm_name.str(std::string());
    opdm_name << "MO-basis Beta OPDM <" << Iroot << "| Etu |" << Jroot << ">";
    SharedMatrix new_OPDM_b(new Matrix(opdm_name.str(), nirrep_, act_dim, act_dim));

    opdm_name.str(std::string());
    opdm_name << "MO-basis OPDM <" << Iroot << "| Etu |" << Jroot << ">";
    SharedMatrix new_OPDM(new Matrix(opdm_name.str(), nirrep_, act_dim, act_dim));

    int offset = 0;
    for (int h=0; h<nirrep_; h++){
      if (!CalcInfo_->ci_orbs[h]) continue;

      double* opdm_a = new_OPDM_a->pointer(h)[0];
      double* opdm_b = new_OPDM_b->pointer(h)[0];
      double* opdm  = new_OPDM->pointer(h)[0];

      for (int i=0, target=0; i<CalcInfo_->ci_orbs[h]; i++){
        int ni = CalcInfo_->act_reorder[i + offset];
        for (int j=0; j<CalcInfo_->ci_orbs[h]; j++){
          int nj = CalcInfo_->act_reorder[j + offset];

          opdm_a[target] = scratch_ap[ni][nj];
          opdm_b[target] = scratch_bp[ni][nj];
          opdm[target++] = scratch_ap[ni][nj] + scratch_bp[ni][nj];
        }
      }
      offset += CalcInfo_->ci_orbs[h];
    }

    std::vector<SharedMatrix> opdm_root_vec;
    opdm_root_vec.push_back(new_OPDM_a);
    opdm_root_vec.push_back(new_OPDM_b);
    opdm_root_vec.push_back(new_OPDM);

    opdm_list.push_back(opdm_root_vec);


  } /* end loop over states_vec */

  if (transp_tmp != nullptr) free_block(transp_tmp);
  if (transp_tmp2 != nullptr) free_block(transp_tmp2);

  scratch_a.reset();
  scratch_b.reset();
  timer_off("CIWave: opdm");

  return opdm_list;
}

void CIWavefunction::opdm_block(struct stringwr **alplist, struct stringwr **betlist,
    double **onepdm_a, double **onepdm_b, double **CJ, double **CI, int Ja_list,
    int Jb_list, int Jnas, int Jnbs, int Ia_list, int Ib_list,
    int Inas, int Inbs)
{
  int Ia_idx, Ib_idx, Ja_idx, Jb_idx, Ja_ex, Jb_ex, Jbcnt, Jacnt;
  struct stringwr *Jb, *Ja;
  signed char *Jbsgn, *Jasgn;
  unsigned int *Jbridx, *Jaridx;
  double C1, C2, Ib_sgn, Ia_sgn;
  int i, j, oij, *Jboij, *Jaoij;

  /* loop over Ia in Ia_list */
  if (Ia_list == Ja_list) {
    for (Ia_idx = 0; Ia_idx < Inas; Ia_idx++) {
      for (Jb = betlist[Jb_list], Jb_idx = 0; Jb_idx < Jnbs; Jb_idx++, Jb++) {
        C1 = CJ[Ia_idx][Jb_idx];

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
          i = oij / CalcInfo_->num_ci_orbs;
          j = oij % CalcInfo_->num_ci_orbs;
          onepdm_b[i][j] += C1 * C2 * Ib_sgn;
        }
      }
    }
  }

  /* loop over Ib in Ib_list */
  if (Ib_list == Jb_list) {
    for (Ib_idx = 0; Ib_idx < Inbs; Ib_idx++) {
      for (Ja = alplist[Ja_list], Ja_idx = 0; Ja_idx < Jnas; Ja_idx++, Ja++) {
        C1 = CJ[Ja_idx][Ib_idx];

        /* loop over excitations */
        Jacnt = Ja->cnt[Ia_list];
        Jaridx = Ja->ridx[Ia_list];
        Jasgn = Ja->sgn[Ia_list];
        Jaoij = Ja->oij[Ia_list];
        for (Ja_ex = 0; Ja_ex < Jacnt; Ja_ex++) {
          oij = *Jaoij++;
          Ia_idx = *Jaridx++;
          Ia_sgn = (double)*Jasgn++;
          C2 = CI[Ia_idx][Ib_idx];
          i = oij / CalcInfo_->num_ci_orbs;
          j = oij % CalcInfo_->num_ci_orbs;
          onepdm_a[i][j] += C1 * C2 * Ia_sgn;
        }
      }
    }
  }

}// End OPDM Block

void CIWavefunction::ci_nat_orbs() {
  // We can only do restricted orbitals
  outfile->Printf("\n   Computing CI Natural Orbitals\n");
  outfile->Printf("   !Warning: New orbitals will be sorted by occuption number,\n");
  outfile->Printf("   orbital spaces (occ/act/vir) may change.\n");

  SharedMatrix NO_vecs(new Matrix("OPDM Eigvecs", Da_->rowspi(), Da_->colspi()));
  SharedVector NO_occ(new Vector("OPDM Occuption", Da_->rowspi()));

  SharedMatrix D = opdm_add_inactive(opdm_, 2.0, true);
  D->diagonalize(NO_vecs, NO_occ, descending);

  Ca_ = Matrix::doublet(Ca_, NO_vecs);
  Cb_ = Ca_;
}

  // for the edification of the autodoc-er. these set py-side or through oeprop calls
  /*- Process::environment.globals["CIn TOTAL ENERGY"] -*/
  /*- Process::environment.globals["CIn CORRELATION ENERGY"] -*/
  /*- Process::environment.globals["CI ROOT n TOTAL ENERGY"] -*/
  /*- Process::environment.globals["CI ROOT n CORRELATION ENERGY"] -*/
  /*- Process::environment.globals["CI ROOT n DIPOLE X"] -*/
  /*- Process::environment.globals["CI ROOT n DIPOLE Y"] -*/
  /*- Process::environment.globals["CI ROOT n DIPOLE Z"] -*/
  /*- Process::environment.globals["CI ROOT n QUADRUPOLE XX"] -*/
  /*- Process::environment.globals["CI ROOT n QUADRUPOLE XY"] -*/
  /*- Process::environment.globals["CI ROOT n QUADRUPOLE XZ"] -*/
  /*- Process::environment.globals["CI ROOT n QUADRUPOLE YY"] -*/
  /*- Process::environment.globals["CI ROOT n QUADRUPOLE YZ"] -*/
  /*- Process::environment.globals["CI ROOT n QUADRUPOLE ZZ"] -*/
  /*- Process::environment.globals["CI ROOT n -> ROOT m DIPOLE X"] -*/
  /*- Process::environment.globals["CI ROOT n -> ROOT m DIPOLE Y"] -*/
  /*- Process::environment.globals["CI ROOT n -> ROOT m DIPOLE Z"] -*/
  /*- Process::environment.globals["CI ROOT n -> ROOT m QUADRUPOLE XX"] -*/
  /*- Process::environment.globals["CI ROOT n -> ROOT m QUADRUPOLE XY"] -*/
  /*- Process::environment.globals["CI ROOT n -> ROOT m QUADRUPOLE XZ"] -*/
  /*- Process::environment.globals["CI ROOT n -> ROOT m QUADRUPOLE YY"] -*/
  /*- Process::environment.globals["CI ROOT n -> ROOT m QUADRUPOLE YZ"] -*/
  /*- Process::environment.globals["CI ROOT n -> ROOT m QUADRUPOLE ZZ"] -*/

}} // namespace psi::detci
