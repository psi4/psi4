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

#include <cstdio>
#include <cmath>
#include <cstring>
#include "psi4/psifiles.h"
#include "psi4/libqt/qt.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libqt/slaterdset.h"

#include "psi4/detci/structs.h"
#include "psi4/detci/slaterd.h"
#include "psi4/detci/civect.h"
#include "psi4/detci/ciwave.h"

namespace psi {
namespace detci {

/*
** diag_h(): Function diagonalizes the hamiltonian
**
** Parameters:
**    alplist = list of alpha strings
**    betlist = list of beta strings
**
** Returns: none
*/
int CIWavefunction::diag_h(double conv_e, double conv_rms) {
    BIGINT size;
    int nroots, i, j;
    double *evals, **evecs, nucrep, edrc, tval;
    double *cbuf;
    char e_label[PSIO_KEYLEN]; /* 80... */

    nroots = Parameters_->num_roots;
    Process::environment.globals["DETCI AVG DVEC NORM"] = 0.0;

    if (print_){
        outfile->Printf("\n   ==> Starting CI iterations <==\n\n");
    }

    if (conv_rms < 0) conv_rms = Parameters_->convergence;
    if (conv_e < 0) conv_e = Parameters_->energy_convergence;

    Parameters_->diag_h_converged = false;
    Parameters_->diag_iters_taken = 0;

    size = CIblks_->vectlen;
    if ((BIGINT)Parameters_->nprint > size) Parameters_->nprint = (int)size;
    nucrep = CalcInfo_->enuc;
    edrc = CalcInfo_->edrc;

    if (Parameters_->bendazzoli)
        outfile->Printf("    Bendazzoli algorithm selected for sigma3\n");

    /* Direct Method --- use RSP diagonalization routine */
    if (Parameters_->diag_method == METHOD_RSP) {

        double h_size = (double)(8 * size * size);
        if (h_size > (Process::environment.get_memory() * 0.4)) {
            outfile->Printf(
                "CIWave::Requsted size of the hamiltonian is %4.2lf GB!\n",
                h_size / 1E9);
            throw PSIEXCEPTION(
                "CIWave::hamiltonian: Size is too large for explicit "
                "hamiltonian build. Increase memory to avoid this error.");
        }

        SharedCIVector Dvec = D_vector();
        Dvec->init_io_files(true);

        if (print_) {
            outfile->Printf("    Exact Diagonalization Algorithm\n");
        }

        SharedMatrix H = hamiltonian();
        SharedMatrix evecs(new Matrix("CI Eigenvectors", (size_t)size, (size_t)size));
        SharedVector evals_v(new Vector("CI Eigenvalues", (size_t)size));

        if (print_ > 4 && size < 200) {
            outfile->Printf("    Hamiltonian matrix:\n");
            H->print();
        }

        H->diagonalize(evecs, evals_v, ascending);

        // Write evecs to Dvec
        evals = init_array(nroots);
        double* tmp_buff = new double[size];
        double** evecsp = evecs->pointer();
        double* evals_vp = evals_v->pointer();
        for (size_t root = 0; root < nroots; root++) {
            for (size_t i = 0; i < size; i++) {
                tmp_buff[i] = evecsp[i][root];
            }
            Dvec->setarray(tmp_buff, size);
            Dvec->write(root, 0);

            // evals doesnt want edrc, but H *has* drc in it
            evals[root] = evals_vp[root] - edrc;
        }
        delete[] tmp_buff;

        // Init sigma for the future
        sigma_init(*(Dvec.get()), *(Dvec.get()));
        Dvec->close_io_files(1);

        H.reset();
        evecs.reset();
        evals_v.reset();
        Parameters_->diag_h_converged = true;
    } /* end RSP section */

    /* RSP test of Davidson/Liu (SEM) diagonalization routine */
    else if (Parameters_->diag_method == METHOD_RSPTEST_OF_SEM) {
        // in-core CIvectors, shouldn't need to open files
        CIvect Cvec(1, 1, 0, 0, CIblks_, CalcInfo_, Parameters_, H0block_);
        CIvect Hd(1, 1, 0, 0, CIblks_, CalcInfo_, Parameters_, H0block_);

        double **H, **b;
        int Iarel, Ialist, Ibrel, Iblist, ij, k, l, tmpi, L;
        unsigned long int ii, jj;
        SlaterDeterminant I, J;
        int *mi_iac, *mi_ibc, *mi_iaidx, *mi_ibidx;
        double *mi_coeff;
        int sm_tridim;
        double *sm_evals, *sm_mat, **sm_evecs, tval;

        if (print_) {
            outfile->Printf("    Using the SEM Test Method\n");
            outfile->Printf("    (n.b. this is for debugging purposes only!)\n");
        }

        /* get the diagonal elements of H into an array Hd */

        Hd.diag_mat_els(alplist_, betlist_, CalcInfo_->onel_ints->pointer(),
                        CalcInfo_->twoel_ints->pointer(), edrc, CalcInfo_->num_alp_expl,
                        CalcInfo_->num_bet_expl, CalcInfo_->num_ci_orbs,
                        Parameters_->hd_ave);

        /* get the biggest elements and put in H0block */
        if (H0block_->size) {
            Hd.max_abs_vals(H0block_->size, H0block_->alplist,
                            H0block_->betlist, H0block_->alpidx,
                            H0block_->betidx, H0block_->H00,
                            Parameters_->neg_only);
        }

        H0block_setup(CIblks_->num_blocks, CIblks_->Ia_code, CIblks_->Ib_code);
        if (Parameters_->hd_ave) {
            H0block_spin_cpl_chk();
            if (H0block_->osize - H0block_->size) {
                outfile->Printf(
                    "H0block size reduced by %d to %d to ensure"
                    "completion of spin-coupling sets\n",
                    (H0block_->osize - H0block_->size), H0block_->size);
            }
        }
        if (Parameters_->Ms0) {
            H0block_pairup(0);
            if (H0block_->osize - H0block_->size) {
                outfile->Printf(
                    "H0block size reduced by %d to ensure pairing.\n",
                    (H0block_->osize - H0block_->size));
            }
        }

        if (print_ > 4 && Parameters_->hd_otf == FALSE) {
            outfile->Printf("\n    Diagonal elements of the Hamiltonian\n");
            Hd.print();
        }

        if (H0block_->size) {
            H0block_fill();
        }

        if (print_ > 2 && H0block_->size) {
            H0block_print();
        }

        if (print_ > 3 && H0block_->size) {
            outfile->Printf("\n\nH0 Block:\n");
            print_mat(H0block_->H0b, H0block_->size, H0block_->size, "outfile");
        }

        H = init_matrix(size, size);
        evals = init_array(size);
        for (ii = 0; ii < size; ii++) {
            Cvec.det2strings(ii, &Ialist, &Iarel, &Iblist, &Ibrel);
            I.set(CalcInfo_->num_alp_expl, alplist_[Ialist][Iarel].occs,
                  CalcInfo_->num_bet_expl, betlist_[Iblist][Ibrel].occs);
            /* introduce symmetry or other restrictions here */
            for (jj = 0; jj < ii; jj++) {
                Cvec.det2strings(jj, &Ialist, &Iarel, &Iblist, &Ibrel);
                J.set(CalcInfo_->num_alp_expl, alplist_[Ialist][Iarel].occs,
                      CalcInfo_->num_bet_expl, betlist_[Iblist][Ibrel].occs);
                H[ii][jj] = H[jj][ii] = matrix_element(&I, &J);
            }
            Cvec.det2strings(jj, &Ialist, &Iarel, &Iblist, &Ibrel);
            J.set(CalcInfo_->num_alp_expl, alplist_[Ialist][Iarel].occs,
                  CalcInfo_->num_bet_expl, betlist_[Iblist][Ibrel].occs);
            H[jj][jj] = matrix_element(&J, &J) + CalcInfo_->edrc;
        }

        /* obtain a set of L orthonormal trial vectors, L > nroots */
        b = (double **)malloc(Parameters_->maxnvect * sizeof(double *));
        for (i = 0; i < Parameters_->maxnvect; i++) {
            if (i < Parameters_->num_init_vecs)
                b[i] = init_array(size);
            else
                b[i] = NULL;
        }

        evecs = init_matrix(Parameters_->num_roots, size);

        L = H0block_->size;
        sm_tridim = L * (L + 1) / 2;
        sm_mat = init_array(sm_tridim);
        sm_evals = init_array(L);
        sm_evecs = init_matrix(L, L);
        for (i = 0, ij = 0; i < L; i++)
            for (j = 0; j <= i; j++, ij++) sm_mat[ij] = H0block_->H0b[i][j];
        rsp(L, L, sm_tridim, sm_mat, sm_evals, 1, sm_evecs, 1E-14);

        /*
        if (Parameters_->precon == PRECON_GEN_DAVIDSON) {
          for (i=0; i<H0block_->size; i++) {
             H0block_->H0b_eigvals[i] = sm_evals[i];
             for (j=0; j<H0block_->size; i++)
                H0block_->H0b_diag[i][j] = sm_evecs[i][j];
             }
          }
        */

        /* need to fill out sm_evecs into b (pad w/ 0's) */
        cbuf = *(Cvec.blockptr(0));
        Cvec.buf_unlock();
        for (i = 0, k = 0; i < L && k < Parameters_->num_init_vecs; i++) {
            /* check sm_evecs[i] to see if it has the correct spin symmetry */
            for (j = 0, tmpi = 0; Parameters_->Ms0 && j < L && !tmpi; j++) {
                l = H0block_->pair[j];
                if (l == -1) {
                    throw PsiException(
                        "(diag_h sem_test): unpaired H0block member!", __FILE__,
                        __LINE__);
                }
                tval = sm_evecs[l][i];
                if ((int)Parameters_->S % 2) tval = -tval;
                if (sm_evecs[j][i] - tval > 1.0E-12) tmpi = 1;
            }
            if (tmpi) continue;

            for (j = 0; j < L; j++) sm_evals[j] = sm_evecs[j][i];

            Cvec.buf_lock(b[k]);
            Cvec.init_vals(0, L, H0block_->alplist, H0block_->alpidx,
                           H0block_->betlist, H0block_->betidx,
                           H0block_->blknum, sm_evals);
            Cvec.buf_unlock();
            k++;
        }
        Cvec.buf_lock(cbuf);
        free(sm_mat);
        free(sm_evals);
        free_matrix(sm_evecs, L);
        for (i = k; i < Parameters_->num_init_vecs; i++) free(b[i]);
        L = k;
        if (L < Parameters_->num_roots) {
            throw PsiException("(diag_h sem_test): Ooops! L < num_roots!",
                               __FILE__, __LINE__);
        }
        sem_test(H, size, Parameters_->num_roots, L, evecs, evals, b, conv_e,
                 conv_rms, Parameters_->maxiter, (nucrep + CalcInfo_->edrc), &i,
                 Parameters_->maxnvect);

        outfile->Printf("SEM used %d expansion vectors\n", i);

        if (print_ > 4) {
            eivout(evecs, evals, size, nroots, "outfile");
        }
        free_matrix(H, size);

        if (print_) {
            mi_iac = init_int_array(Parameters_->nprint);
            mi_ibc = init_int_array(Parameters_->nprint);
            mi_iaidx = init_int_array(Parameters_->nprint);
            mi_ibidx = init_int_array(Parameters_->nprint);
            mi_coeff = init_array(Parameters_->nprint);

            for (i = 0; i < nroots; i++) {
                outfile->Printf("\n\n* ROOT %2d CI total energy = %17.13lf\n",
                                i, evals[i] + nucrep);
                Cvec.setarray(evecs[i], size);
                zero_arr(mi_coeff, Parameters_->nprint);
                Cvec.max_abs_vals(Parameters_->nprint, mi_iac, mi_ibc, mi_iaidx,
                                  mi_ibidx, mi_coeff, Parameters_->neg_only);
                print_vec(Parameters_->nprint, mi_iac, mi_ibc, mi_iaidx,
                          mi_ibidx, mi_coeff);
            }
            free(mi_iac);
            free(mi_ibc);
            free(mi_iaidx);
            free(mi_ibidx);
            free(mi_coeff);
            free_matrix(evecs, Parameters_->num_roots);
        }
        Parameters_->diag_h_converged = true;

    } /* end Test of Davidson/Liu section */

    /*
     * Davidson/Liu Simultaneous Expansion Method OR
     * Mitrushenkov's Olsen-modified Davidson Algorithm
     */

    else {
        /* prepare the H0 block */

        CIvect Hd(Parameters_->icore, 1, 1,
                  Parameters_->hd_filenum, CIblks_, CalcInfo_, Parameters_,
                  H0block_);

        bool open_old = false;
        if (Parameters_->restart) open_old = true;
        Hd.init_io_files(open_old);

        /* get the diagonal elements of H into an array Hd */
        if (!Parameters_->restart ||
            (Parameters_->restart && Parameters_->hd_otf)) {
            if (print_ > 1) {
                outfile->Printf("\n    Forming diagonal elements of H\n");
            }
            Hd.diag_mat_els(alplist_, betlist_, CalcInfo_->onel_ints->pointer(),
                            CalcInfo_->twoel_ints->pointer(), edrc,
                            CalcInfo_->num_alp_expl, CalcInfo_->num_bet_expl,
                            CalcInfo_->num_ci_orbs, Parameters_->hd_ave);
        } else {
            Hd.read(0, 0);
        }

        /* get the biggest elements and put in H0block */
        if (H0block_->size) {
            if (print_ > 1) {
                outfile->Printf("\n    Forming H0 block\n");
            }

            if (!Parameters_->hd_otf)
                Hd.max_abs_vals(H0block_->size + H0block_->coupling_size,
                                H0block_->alplist, H0block_->betlist,
                                H0block_->alpidx, H0block_->betidx,
                                H0block_->H00, Parameters_->neg_only);
        }

        if (Parameters_->hd_otf) psio_close(Parameters_->hd_filenum, 1);

        H0block_setup(CIblks_->num_blocks, CIblks_->Ia_code, CIblks_->Ib_code);
        if (Parameters_->filter_guess) H0block_filter_setup();
        if (Parameters_->hd_ave) {
            H0block_spin_cpl_chk();
            if ((H0block_->osize - H0block_->size) &&
                print_ > 1) {
                outfile->Printf(
                    "H0block size reduced by %d to ensure "
                    "completion of spin-coupling sets\n",
                    (H0block_->osize - H0block_->size));
                H0block_->osize = H0block_->size;
            }
            if ((H0block_->oguess_size - H0block_->guess_size) &&
                print_ > 1) {
                outfile->Printf(
                    "    H0block guess size reduced by %d to ensure "
                    "completion of spin-coupling sets\n",
                    (H0block_->oguess_size - H0block_->guess_size));
                H0block_->oguess_size = H0block_->guess_size;
            }
            if ((H0block_->ocoupling_size - H0block_->coupling_size) &&
                print_ > 1) {
                outfile->Printf(
                    "    H0block coupling size reduced by %d to ensure "
                    "completion of spin-coupling sets\n",
                    (H0block_->ocoupling_size - H0block_->coupling_size));
                H0block_->ocoupling_size = H0block_->coupling_size;
            }
        }
        if (Parameters_->Ms0) {
            /* if (H0block_->guess_size < H0block_->size) */
            H0block_pairup(0); /* pairup h0block size */
            H0block_pairup(1); /* pairup guess_size */
            H0block_pairup(2); /* pairup coupling size */
            if ((H0block_->osize - H0block_->size) &&
                print_ > 1) {
                outfile->Printf(
                    "    H0block size reduced by %d to ensure pairing"
                    "and spin-coupling.\n",
                    (H0block_->osize - H0block_->size));
            }
            if ((H0block_->oguess_size - H0block_->guess_size) &&
                print_ > 1) {
                outfile->Printf(
                    "    H0block guess size reduced by %d to "
                    "ensure pairing and spin-coupling.\n",
                    (H0block_->oguess_size - H0block_->guess_size));
            }
            if ((H0block_->ocoupling_size - H0block_->coupling_size) &&
                print_ > 1) {
                outfile->Printf(
                    "    H0block coupling size reduced by %d to "
                    "ensure pairing and spin-coupling.\n",
                    (H0block_->ocoupling_size - H0block_->coupling_size));
            }
        }

        Parameters_->neg_only = 0; /* MLL 7-2-97 */
        if (print_ > 4) {
            outfile->Printf("\n    Diagonal elements of the Hamiltonian\n");
            Hd.print();
        }

        if (H0block_->size) {
            H0block_fill();
        }

        if (print_ > 2 && H0block_->size) {
            H0block_print();
        }

        if (print_ > 3 && H0block_->size) {
            outfile->Printf("\n\nH0 Block:\n");
            print_mat(H0block_->H0b, H0block_->size, H0block_->size, "outfile");
        }

        /* Davidson/Liu Simultaneous Expansion Method */
        if (Parameters_->diag_method == METHOD_DAVIDSON_LIU_SEM) {
            if (print_) {
                outfile->Printf("\n    Simultaneous Expansion Method (Block Davidson Method)\n");
            }

            evals = init_array(nroots);

            sem_iter(Hd, alplist_, betlist_, evals, conv_e, conv_rms, nucrep,
                     edrc, nroots, Parameters_->maxiter, Parameters_->maxnvect);
        }

        /* Mitrushenkov's Olsen Method */
        else {
            if (print_) {
                if (Parameters_->diag_method == METHOD_MITRUSHENKOV)
                    outfile->Printf("\n    Mitrushenkov's two vector algorithm\n");
                else if (Parameters_->diag_method == METHOD_OLSEN)
                    outfile->Printf("\n    Olsen's single vector algorithm\n");
            }

            evals = init_array(nroots);

            mitrush_iter(Hd, alplist_, betlist_, nroots, evals, conv_rms,
                         conv_e, nucrep, edrc, Parameters_->maxiter,
                         Parameters_->maxnvect);

        }

    } /* end the Davidson-Liu/Mitrushenkov-Olsen-Davidson section */

    // Check convergence
    if (!Parameters_->diag_h_converged){
       convergence_death();
       if (print_){
           outfile->Printf("\n    Warning! CI diagonalization did not fully converge!\n\n");
       }
    }



    // Write out energy to psivars
    tval = evals[Parameters_->root] + edrc + nucrep;

    Process::environment.globals["CURRENT ENERGY"] = tval;
    Process::environment.globals["CURRENT CORRELATION ENERGY"] =
        tval - CalcInfo_->escf;
    Process::environment.globals["CURRENT REFERENCE ENERGY"] = CalcInfo_->escf;
    Process::environment.globals["CI TOTAL ENERGY"] = tval;
    Process::environment.globals["CI CORRELATION ENERGY"] =
        tval - CalcInfo_->escf;

    if (Parameters_->fci) {
        Process::environment.globals["FCI TOTAL ENERGY"] = tval;
        Process::environment.globals["FCI CORRELATION ENERGY"] =
            tval - CalcInfo_->escf;
    } else {
        if (Parameters_->ex_lvl == 2) {
            Process::environment.globals["CISD TOTAL ENERGY"] = tval;
            Process::environment.globals["CISD CORRELATION ENERGY"] =
                tval - CalcInfo_->escf;
        } else if (Parameters_->ex_lvl == 3) {
            Process::environment.globals["CISDT TOTAL ENERGY"] = tval;
            Process::environment.globals["CISDT CORRELATION ENERGY"] =
                tval - CalcInfo_->escf;
        } else if (Parameters_->ex_lvl == 4) {
            Process::environment.globals["CISDTQ TOTAL ENERGY"] = tval;
            Process::environment.globals["CISDTQ CORRELATION ENERGY"] =
                tval - CalcInfo_->escf;
        } else {
            std::stringstream s;
            s << "CI" << Parameters_->ex_lvl << " TOTAL ENERGY";
            Process::environment.globals[s.str()] = tval;
            s.str(std::string());
            s << "CI" << Parameters_->ex_lvl << " CORRELATION ENERGY";
            Process::environment.globals[s.str()] = tval - CalcInfo_->escf;
        }
    }

    for (i = 0; i < nroots; i++) {
        sprintf(e_label, "Root %2d energy", i);
        tval = evals[i] + edrc + nucrep;

        std::stringstream s;
        s << "CI ROOT " << i << " TOTAL ENERGY";
        Process::environment.globals[s.str()] = tval;
        s.str(std::string());
        s << "CI ROOT " << i << " CORRELATION ENERGY";
        Process::environment.globals[s.str()] = tval - CalcInfo_->escf;
    }

    if (Parameters_->average_num > 1) {
        tval = 0.0;
        for (i = 0; i < Parameters_->average_num; i++) {
            tval += Parameters_->average_weights[i] *
                    (edrc + nucrep + evals[Parameters_->average_states[i]]);
        }
        Process::environment.globals["CI STATE-AVERAGED TOTAL ENERGY"] = tval;
        Process::environment.globals["CI STATE-AVERAGED CORRELATION ENERGY"] =
            tval - CalcInfo_->escf;
        Process::environment.globals["CURRENT CORRELATION ENERGY"] =
            Process::environment.globals["CI STATE-AVERAGED CORRELATION ENERGY"];
    }

    // Set the energy as MCSCF would find it
    if (Parameters_->average_num > 1) {  // state average
        Process::environment.globals["MCSCF TOTAL ENERGY"] =
            Process::environment.globals["CI STATE-AVERAGED TOTAL ENERGY"];
    } else if (Parameters_->root != 0) {  // follow some specific root != lowest
        std::stringstream s;
        s << "CI ROOT " << Parameters_->root << " TOTAL ENERGY";
        Process::environment.globals["MCSCF TOTAL ENERGY"] =
            Process::environment.globals[s.str()];
    } else {
        Process::environment.globals["MCSCF TOTAL ENERGY"] =
            Process::environment.globals["CI TOTAL ENERGY"];
    }

    return Parameters_->diag_iters_taken;


}  // end CIWave::diag_h
}}  // namespace psi::detci
