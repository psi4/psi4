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

/** Standard library includes */
#include <fstream>
#include "psi4/psifiles.h"
#include "psi4/libiwl/iwl.hpp"
#include "psi4/libqt/qt.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/sieve.h"
#include "psi4/libmints/integral.h"
#include "dfocc.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace psi;
using namespace std;

namespace psi{ namespace dfoccwave{

void DFOCC::oei_grad()
{

/********************************************************************************************/
/************************** Gradient ********************************************************/
/********************************************************************************************/
    //outfile->Printf("\tComputing analytic gradients...\n");
    //

    // Pointers
    double** Dp = G1ao->to_block_matrix();
    double** Wp = GFao->to_block_matrix();

/********************************************************************************************/
/************************** Nuclear Gradient ************************************************/
/********************************************************************************************/
    // => Nuclear Gradient <= //
    gradients["Nuclear"] = SharedMatrix(molecule_->nuclear_repulsion_energy_deriv1().clone());
    gradients["Nuclear"]->set_name("Nuclear Gradient");
    gradients["Nuclear"]->print_atom_vector();

/********************************************************************************************/
/************************** Kinetic Gradient ************************************************/
/********************************************************************************************/
    // => Kinetic Gradient <= //
    timer_on("Grad: T");
    {

        gradients["Kinetic"] = SharedMatrix(gradients["Nuclear"]->clone());
        gradients["Kinetic"]->set_name("Kinetic Gradient");
        gradients["Kinetic"]->zero();
        double** Tp = gradients["Kinetic"]->pointer();

        // Kinetic derivatives
        std::shared_ptr<OneBodyAOInt> Tint(integral_->ao_kinetic(1));
        const double* buffer = Tint->buffer();

        for (int P = 0; P < basisset_->nshell(); P++) {
            for (int Q = 0; Q <= P; Q++) {

                Tint->compute_shell_deriv1(P,Q);

                int nP = basisset_->shell(P).nfunction();
                int oP = basisset_->shell(P).function_index();
                int aP = basisset_->shell(P).ncenter();

                int nQ = basisset_->shell(Q).nfunction();
                int oQ = basisset_->shell(Q).function_index();
                int aQ = basisset_->shell(Q).ncenter();

                int offset = nP * nQ;
                const double* ref = buffer;
                double perm = (P == Q ? 1.0 : 2.0);

                // Px
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        Tp[aP][0] += perm * Dp[p + oP][q + oQ] * (*ref++);
                    }
                }

                // Py
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        Tp[aP][1] += perm * Dp[p + oP][q + oQ] * (*ref++);
                    }
                }

                // Pz
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        Tp[aP][2] += perm * Dp[p + oP][q + oQ] * (*ref++);
                    }
                }

                // Qx
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        Tp[aQ][0] += perm * Dp[p + oP][q + oQ] * (*ref++);
                    }
                }

                // Qy
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        Tp[aQ][1] += perm * Dp[p + oP][q + oQ] * (*ref++);
                    }
                }

                // Qz
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        Tp[aQ][2] += perm * Dp[p + oP][q + oQ] * (*ref++);
                    }
                }
            }
        }
    }
    gradients["Kinetic"]->print_atom_vector();
    timer_off("Grad: T");

/********************************************************************************************/
/************************** Potential Gradient **********************************************/
/********************************************************************************************/
    // => Potential Gradient <= //
    timer_on("Grad: V");
    {
        gradients["Potential"] = SharedMatrix(gradients["Nuclear"]->clone());
        gradients["Potential"]->set_name("Potential Gradient");
        gradients["Potential"]->zero();

        // Thread count
        int threads = 1;
        #ifdef _OPENMP
            threads = Process::environment.get_n_threads();
        #endif

        // Potential derivatives
        std::vector<std::shared_ptr<OneBodyAOInt> > Vint;
        std::vector<SharedMatrix> Vtemps;
        for (int t = 0; t < threads; t++) {
            Vint.push_back(std::shared_ptr<OneBodyAOInt>(integral_->ao_potential(1)));
            Vtemps.push_back(SharedMatrix(gradients["Potential"]->clone()));
        }

        // Lower Triangle
        std::vector<std::pair<int,int> > PQ_pairs;
        for (int P = 0; P < basisset_->nshell(); P++) {
            for (int Q = 0; Q <= P; Q++) {
                PQ_pairs.push_back(std::pair<int,int>(P,Q));
            }
        }

        #pragma omp parallel for schedule(dynamic) num_threads(threads)
        for (long int PQ = 0L; PQ < PQ_pairs.size(); PQ++) {

            int P = PQ_pairs[PQ].first;
            int Q = PQ_pairs[PQ].second;

            int thread = 0;
            #ifdef _OPENMP
                thread = omp_get_thread_num();
            #endif

            Vint[thread]->compute_shell_deriv1(P,Q);
            const double* buffer = Vint[thread]->buffer();

            int nP = basisset_->shell(P).nfunction();
            int oP = basisset_->shell(P).function_index();
            int aP = basisset_->shell(P).ncenter();

            int nQ = basisset_->shell(Q).nfunction();
            int oQ = basisset_->shell(Q).function_index();
            int aQ = basisset_->shell(Q).ncenter();

            double perm = (P == Q ? 1.0 : 2.0);

            double** Vp = Vtemps[thread]->pointer();

            for (int A = 0; A < natom; A++) {
                const double* ref0 = &buffer[3 * A * nP * nQ + 0 * nP * nQ];
                const double* ref1 = &buffer[3 * A * nP * nQ + 1 * nP * nQ];
                const double* ref2 = &buffer[3 * A * nP * nQ + 2 * nP * nQ];
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        double Vval = perm * Dp[p + oP][q + oQ];
                        Vp[A][0] += Vval * (*ref0++);
                        Vp[A][1] += Vval * (*ref1++);
                        Vp[A][2] += Vval * (*ref2++);
                    }
                }
            }
        }

        for (int t = 0; t < threads; t++) {
            gradients["Potential"]->add(Vtemps[t]);
        }
    }
    gradients["Potential"]->print_atom_vector();
    timer_off("Grad: V");

/********************************************************************************************/
/************************** Overlap Gradient ************************************************/
/********************************************************************************************/
    // => Overlap Gradient <= //
    timer_on("Grad: S");
    {
        gradients["Overlap"] = SharedMatrix(gradients["Nuclear"]->clone());
        gradients["Overlap"]->set_name("Overlap Gradient");
        gradients["Overlap"]->zero();
        double** Sp = gradients["Overlap"]->pointer();

        // Overlap derivatives
        std::shared_ptr<OneBodyAOInt> Sint(integral_->ao_overlap(1));
        const double* buffer = Sint->buffer();

        for (int P = 0; P < basisset_->nshell(); P++) {
            for (int Q = 0; Q <= P; Q++) {

                Sint->compute_shell_deriv1(P,Q);

                int nP = basisset_->shell(P).nfunction();
                int oP = basisset_->shell(P).function_index();
                int aP = basisset_->shell(P).ncenter();

                int nQ = basisset_->shell(Q).nfunction();
                int oQ = basisset_->shell(Q).function_index();
                int aQ = basisset_->shell(Q).ncenter();

                int offset = nP * nQ;
                const double* ref = buffer;
                double perm = (P == Q ? 1.0 : 2.0);

                // Px
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        Sp[aP][0] -= perm * Wp[p + oP][q + oQ] * (*ref++);
                    }
                }

                // Py
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        Sp[aP][1] -= perm * Wp[p + oP][q + oQ] * (*ref++);
                    }
                }

                // Pz
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        Sp[aP][2] -= perm * Wp[p + oP][q + oQ] * (*ref++);
                    }
                }

                // Qx
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        Sp[aQ][0] -= perm * Wp[p + oP][q + oQ] * (*ref++);
                    }
                }

                // Qy
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        Sp[aQ][1] -= perm * Wp[p + oP][q + oQ] * (*ref++);
                    }
                }

                // Qz
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        Sp[aQ][2] -= perm * Wp[p + oP][q + oQ] * (*ref++);
                    }
                }
            }
        }
    }
    gradients["Overlap"]->print_atom_vector();
    timer_off("Grad: S");

    // mem free
    free_block(Dp);
    free_block(Wp);
/********************************************************************************************/
/********************************************************************************************/
/********************************************************************************************/

//outfile->Printf("\toei_grad is done. \n");
}// end


}} // End Namespaces


