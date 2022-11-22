/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2022 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include "psi4/psifiles.h"
#include "psi4/libqt/qt.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libqt/slaterdset.h"
#include "psi4/libpsi4util/process.h"

#include "psi4/detci/structs.h"
#include "psi4/detci/slaterd.h"
#include "psi4/detci/civect.h"
#include "psi4/detci/ciwave.h"

#include <cstdio>
#include <cmath>
#include <cstring>

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
    size_t size;
    int nroots, i, j;
    double *evals, **evecs, nucrep, edrc, tval;
    double *cbuf;
    char e_label[PSIO_KEYLEN]; /* 80... */

    nroots = Parameters_->num_roots;
    set_scalar_variable("DETCI AVG DVEC NORM", 0.0);

    if (print_) {
        outfile->Printf("\n   ==> Starting CI iterations <==\n\n");
    }

    if (conv_rms < 0) conv_rms = Parameters_->convergence;
    if (conv_e < 0) conv_e = Parameters_->energy_convergence;

    Parameters_->diag_h_converged = false;
    Parameters_->diag_iters_taken = 0;

    size = CIblks_->vectlen;
    if ((size_t)Parameters_->nprint > size) Parameters_->nprint = (int)size;
    nucrep = CalcInfo_->enuc;
    edrc = CalcInfo_->edrc;

    if (Parameters_->bendazzoli) outfile->Printf("    Bendazzoli algorithm selected for sigma3\n");

    /* Direct Method --- use RSP diagonalization routine */
    if (Parameters_->diag_method == METHOD_RSP) {
        double h_size = (double)(8 * size * size);
        if (h_size > (Process::environment.get_memory() * 0.4)) {
            outfile->Printf("CIWave::Requsted size of the hamiltonian is %4.2lf GB!\n", h_size / 1E9);
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
        auto evecs = std::make_shared<Matrix>("CI Eigenvectors", (size_t)size, (size_t)size);
        auto evals_v = std::make_shared<Vector>("CI Eigenvalues", (size_t)size);

        if (print_ > 4 && size < 200) {
            outfile->Printf("    Hamiltonian matrix:\n");
            H->print();
        }

        H->diagonalize(evecs, evals_v, ascending);

        // Write evecs to Dvec
        evals = init_array(nroots);
        auto *tmp_buff = new double[size];
        double **evecsp = evecs->pointer();
        double *evals_vp = evals_v->pointer();
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

    /*
     * Davidson/Liu Simultaneous Expansion Method OR
     * Mitrushenkov's Olsen-modified Davidson Algorithm
     */

    else {
        /* prepare the H0 block */

        CIvect Hd(Parameters_->icore, 1, 1, Parameters_->hd_filenum, CIblks_, CalcInfo_, Parameters_, H0block_);

        bool open_old = false;
        if (Parameters_->restart) open_old = true;
        Hd.init_io_files(open_old);

        /* get the diagonal elements of H into an array Hd */
        if (!Parameters_->restart || (Parameters_->restart && Parameters_->hd_otf)) {
            if (print_ > 1) {
                outfile->Printf("\n    Forming diagonal elements of H\n");
            }
            Hd.diag_mat_els(alplist_, betlist_, CalcInfo_->onel_ints->pointer(), CalcInfo_->twoel_ints->pointer(), edrc,
                            CalcInfo_->num_alp_expl, CalcInfo_->num_bet_expl, CalcInfo_->num_ci_orbs,
                            Parameters_->hd_ave);
        } else {
            Hd.read(0, 0);
        }

        /* get the biggest elements and put in H0block */
        if (H0block_->size) {
            if (print_ > 1) {
                outfile->Printf("\n    Forming H0 block\n");
            }

            if (!Parameters_->hd_otf)
                Hd.max_abs_vals(H0block_->size + H0block_->coupling_size, H0block_->alplist, H0block_->betlist,
                                H0block_->alpidx, H0block_->betidx, H0block_->H00, Parameters_->neg_only);
        }

        if (Parameters_->hd_otf) psio_close(Parameters_->hd_filenum, 1);

        H0block_setup(CIblks_->num_blocks, CIblks_->Ia_code, CIblks_->Ib_code);
        if (Parameters_->filter_guess) H0block_filter_setup();
        if (Parameters_->hd_ave) {
            H0block_spin_cpl_chk();
            if ((H0block_->osize - H0block_->size) && print_ > 1) {
                outfile->Printf(
                    "H0block size reduced by %d to ensure "
                    "completion of spin-coupling sets\n",
                    (H0block_->osize - H0block_->size));
                H0block_->osize = H0block_->size;
            }
            if ((H0block_->oguess_size - H0block_->guess_size) && print_ > 1) {
                outfile->Printf(
                    "    H0block guess size reduced by %d to ensure "
                    "completion of spin-coupling sets\n",
                    (H0block_->oguess_size - H0block_->guess_size));
                H0block_->oguess_size = H0block_->guess_size;
            }
            if ((H0block_->ocoupling_size - H0block_->coupling_size) && print_ > 1) {
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
            if ((H0block_->osize - H0block_->size) && print_ > 1) {
                outfile->Printf(
                    "    H0block size reduced by %d to ensure pairing"
                    "and spin-coupling.\n",
                    (H0block_->osize - H0block_->size));
            }
            if ((H0block_->oguess_size - H0block_->guess_size) && print_ > 1) {
                outfile->Printf(
                    "    H0block guess size reduced by %d to "
                    "ensure pairing and spin-coupling.\n",
                    (H0block_->oguess_size - H0block_->guess_size));
            }
            if ((H0block_->ocoupling_size - H0block_->coupling_size) && print_ > 1) {
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

            sem_iter(Hd, alplist_, betlist_, evals, conv_e, conv_rms, nucrep, edrc, nroots, Parameters_->maxiter,
                     Parameters_->maxnvect);
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

            mitrush_iter(Hd, alplist_, betlist_, nroots, evals, conv_rms, conv_e, nucrep, edrc, Parameters_->maxiter,
                         Parameters_->maxnvect);
        }

    } /* end the Davidson-Liu/Mitrushenkov-Olsen-Davidson section */

    // Check convergence
    if (!Parameters_->diag_h_converged) {
        convergence_death();
        if (print_) {
            outfile->Printf("\n    Warning! CI diagonalization did not fully converge!\n\n");
        }
    }

    // Write out energy to psivars
    tval = evals[Parameters_->root] + edrc + nucrep;

    set_energy(tval);
    set_scalar_variable("CURRENT ENERGY", tval);
    set_scalar_variable("CURRENT CORRELATION ENERGY", tval - CalcInfo_->escf);
    set_scalar_variable("CURRENT REFERENCE ENERGY", CalcInfo_->escf);
    set_scalar_variable("CI TOTAL ENERGY", tval);
    set_scalar_variable("CI CORRELATION ENERGY", tval - CalcInfo_->escf);

    if (Parameters_->fci) {
        set_scalar_variable("FCI TOTAL ENERGY", tval);
        set_scalar_variable("FCI CORRELATION ENERGY", tval - CalcInfo_->escf);
    } else {
        if (Parameters_->ex_lvl == 2) {
            set_scalar_variable("CISD TOTAL ENERGY", tval);
            set_scalar_variable("CISD CORRELATION ENERGY", tval - CalcInfo_->escf);
        } else if (Parameters_->ex_lvl == 3) {
            set_scalar_variable("CISDT TOTAL ENERGY", tval);
            set_scalar_variable("CISDT CORRELATION ENERGY", tval - CalcInfo_->escf);
        } else if (Parameters_->ex_lvl == 4) {
            set_scalar_variable("CISDTQ TOTAL ENERGY", tval);
            set_scalar_variable("CISDTQ CORRELATION ENERGY", tval - CalcInfo_->escf);
        } else {
            std::stringstream s;
            s << "CI" << Parameters_->ex_lvl << " TOTAL ENERGY";
            set_scalar_variable(s.str(), tval);
            s.str(std::string());
            s << "CI" << Parameters_->ex_lvl << " CORRELATION ENERGY";
            set_scalar_variable(s.str(), tval - CalcInfo_->escf);
        }
    }

    for (i = 0; i < nroots; i++) {
        sprintf(e_label, "Root %2d energy", i);
        tval = evals[i] + edrc + nucrep;

        std::stringstream s;
        s << "CI ROOT " << i << " TOTAL ENERGY";
        set_scalar_variable(s.str(), tval);
        s.str(std::string());
        s << "CI ROOT " << i << " CORRELATION ENERGY";
        set_scalar_variable(s.str(), tval - CalcInfo_->escf);
    }

    if (Parameters_->average_num > 1) {
        tval = 0.0;
        for (i = 0; i < Parameters_->average_num; i++) {
            tval += Parameters_->average_weights[i] * (edrc + nucrep + evals[Parameters_->average_states[i]]);
        }
        set_scalar_variable("CI STATE-AVERAGED TOTAL ENERGY", tval);
        set_scalar_variable("CI STATE-AVERAGED CORRELATION ENERGY", tval - CalcInfo_->escf);
        set_scalar_variable("CURRENT CORRELATION ENERGY", scalar_variable("CI STATE-AVERAGED CORRELATION ENERGY"));
    }

    // Set the energy as MCSCF would find it
    if (Parameters_->average_num > 1) {  // state average
        set_scalar_variable("MCSCF TOTAL ENERGY", scalar_variable("CI STATE-AVERAGED TOTAL ENERGY"));
    } else if (Parameters_->root != 0) {  // follow some specific root != lowest
        std::stringstream s;
        s << "CI ROOT " << Parameters_->root << " TOTAL ENERGY";
        set_scalar_variable("MCSCF TOTAL ENERGY", scalar_variable(s.str()));
    } else {
        set_scalar_variable("MCSCF TOTAL ENERGY", scalar_variable("CI TOTAL ENERGY"));
    }

    return Parameters_->diag_iters_taken;

}  // end CIWave::diag_h
}  // namespace detci
}  // namespace psi
