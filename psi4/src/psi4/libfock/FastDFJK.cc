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


#include "psi4/lib3index/3index.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/aiohandler.h"
#include "psi4/libqt/qt.h"
#include "psi4/psi4-dec.h"
#include "psi4/psifiles.h"
#include "psi4/libmints/sieve.h"
#include "psi4/libiwl/iwl.hpp"
#include "jk.h"
#include "jk_independent.h"
#include "link.h"
#include "direct_screening.h"
#include "cubature.h"
#include "points.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/twobody.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/3coverlap.h"
#include "psi4/lib3index/cholesky.h"

#include <sstream>
#include "psi4/libparallel/ParallelPrinter.h"
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace psi;

namespace psi {


FastDFJK::FastDFJK(std::shared_ptr<BasisSet> primary,
   std::shared_ptr<BasisSet> auxiliary) :
   JK(primary), auxiliary_(auxiliary)
{
    common_init();
}
FastDFJK::~FastDFJK()
{
}
void FastDFJK::common_init()
{
    df_ints_num_threads_ = 1;
    #ifdef _OPENMP
        df_ints_num_threads_ = Process::environment.get_n_threads();
    #endif
    df_ints_io_ = "NONE";
    condition_ = 1.0E-12;
    unit_ = PSIF_DFSCF_BJ;
    is_core_ = true;
    psio_ = PSIO::shared_object();
    metric_ = "COULOMB";
    theta_ = 1.0;
    domains_ = "DIATOMIC";
    bump_R0_ = 0.0;
    bump_R1_ = 0.0;

}
void FastDFJK::print_header() const
{
    if (print_) {
        outfile->Printf( "  ==> FastDFJK: Density-Fitted J/K Matrices <==\n\n");

        outfile->Printf( "    J tasked:          %11s\n", (do_J_ ? "Yes" : "No"));
        outfile->Printf( "    K tasked:          %11s\n", (do_K_ ? "Yes" : "No"));
        outfile->Printf( "    wK tasked:         %11s\n", (do_wK_ ? "Yes" : "No"));
        if (do_wK_)
            outfile->Printf( "    Omega:             %11.3E\n", omega_);
        outfile->Printf( "    OpenMP threads:    %11d\n", omp_nthread_);
        outfile->Printf( "    Integrals threads: %11d\n", df_ints_num_threads_);
        outfile->Printf( "    Memory (MB):       %11ld\n", (memory_ *8L) / (1024L * 1024L));
        outfile->Printf( "    Algorithm:         %11s\n",  (is_core_ ? "Core" : "Disk"));
        outfile->Printf( "    Integral Cache:    %11s\n",  df_ints_io_.c_str());
        outfile->Printf( "    Schwarz Cutoff:    %11.0E\n", cutoff_);
        outfile->Printf( "    Fitting Condition: %11.0E\n", condition_);
        outfile->Printf( "    Fitting Metric:    %11s\n", metric_.c_str());
        if (metric_ == "EWALD")
            outfile->Printf( "    Theta:             %11.3E\n", theta_);
        outfile->Printf( "    Fitting Domains:   %11s\n", domains_.c_str());
        if (domains_ != "DIATOMIC") {
            outfile->Printf( "    Bump R0:           %11.3E\n", bump_R0_);
            outfile->Printf( "    Bump R1:           %11.3E\n", bump_R1_);
        }

        outfile->Printf( "\n");

        outfile->Printf( "   => Auxiliary Basis Set <=\n\n");
        auxiliary_->print_by_level("outfile", print_);
    }
}
void FastDFJK::preiterations()
{
    // DF requires constant sieve, must be static throughout object life
    if (!sieve_) {
        sieve_ = std::shared_ptr<ERISieve>(new ERISieve(primary_, cutoff_));
    }

    build_atom_pairs();
    build_shell_pairs();
    build_auxiliary_partition();
    build_Bpq();

    Z_ = build_Z(0.0);
    if (do_wK_)
        Z_LR_ = build_Z(omega_);

    if (print_ > 1) {
        outfile->Printf("  ==> Atom Pair Tasks <==\n\n");
        for (size_t pair = 0L; pair < atom_pairs_.size(); pair++) {
            outfile->Printf("  Task %8zu: Atom Pair (%4d,%4d)\n", pair, atom_pairs_[pair].first, atom_pairs_[pair].second);
            const std::vector<int>& auxiliary_atoms = auxiliary_atoms_[pair];
            const std::vector<double>& bump_atoms = bump_atoms_[pair];
            outfile->Printf("   Auxiliary Atoms/Bump Functions:\n");
            for (size_t A = 0; A < auxiliary_atoms.size(); A++) {
                outfile->Printf("    %4d: %11.3E\n", auxiliary_atoms[A],bump_atoms[A]);
            }
            if (debug_ > 1) {
                outfile->Printf("   Primary Shell Pairs:\n");
                const std::vector<std::pair<int,int> >& shell_pairs = shell_pairs_[pair];
                for (size_t A = 0; A < shell_pairs.size(); A++) {
                    outfile->Printf("    (%4d,%4d)\n", shell_pairs[A].first, shell_pairs[A].second);
                }
            }
        }
        outfile->Printf( "\n");
    }

    if (debug_ > 3) {
        outfile->Printf("  ==> Atom Pair Tensors <==\n\n");
        for (size_t pair = 0L; pair < atom_pairs_.size(); pair++) {
            outfile->Printf("  Task %8zu: Atom Pair (%4d,%4d)\n", pair, atom_pairs_[pair].first, atom_pairs_[pair].second);
            Bpq_[pair]->print();
        }
        outfile->Printf( "\n");
        Z_->print();
        if (do_wK_) Z_LR_->print();
    }
}
void FastDFJK::compute_JK()
{
    if (do_J_) {
        build_J(Z_,D_ao_,J_ao_);
    }
    if (do_K_) {
        build_K(Z_,D_ao_,K_ao_);
    }
    if (do_wK_) {
        build_K(Z_LR_,D_ao_,wK_ao_);
    }
}
void FastDFJK::postiterations()
{
    Z_.reset();
    Z_LR_.reset();
    atom_pairs_.clear();
    shell_pairs_.clear();
    auxiliary_atoms_.clear();
    bump_atoms_.clear();
    Bpq_.clear();
}
std::shared_ptr<Matrix> FastDFJK::build_Z(double omega)
{
    int naux = auxiliary_->nbf();
    int nshell = auxiliary_->nshell();
    int nthread = 1;
    #ifdef _OPENMP
        nthread = df_ints_num_threads_;
    #endif

    std::shared_ptr<IntegralFactory> fact(new IntegralFactory(auxiliary_,BasisSet::zero_ao_basis_set(),auxiliary_,BasisSet::zero_ao_basis_set()));
    std::vector<std::shared_ptr<TwoBodyAOInt> > ints;
    for (int thread = 0; thread < nthread; thread++) {
        if (omega != 0.0) {
            ints.push_back(std::shared_ptr<TwoBodyAOInt>(fact->erf_eri(omega)));
        } else {
            ints.push_back(std::shared_ptr<TwoBodyAOInt>(fact->eri()));
        }
    }

    std::shared_ptr<Matrix> Z(new Matrix((omega != 0.0 ? "Z_LR" : "Z"), naux, naux));
    double** Zp = Z->pointer();

    #pragma omp parallel for schedule(dynamic) num_threads(nthread)
    for (int P = 0; P < nshell; P++) {
        int thread = 0;
        #ifdef _OPENMP
            thread = omp_get_thread_num();
        #endif
        const double* buffer = ints[thread]->buffer();
        for (int Q = 0; Q <= P; Q++) {
            ints[thread]->compute_shell(P,0,Q,0);
            int nP = auxiliary_->shell(P).nfunction();
            int nQ = auxiliary_->shell(Q).nfunction();
            int oP = auxiliary_->shell(P).function_index();
            int oQ = auxiliary_->shell(Q).function_index();
            for (int p = 0, index = 0; p < nP; p++) {
                for (int q = 0; q < nQ; q++, index++) {
                    Zp[p + oP][q + oQ] = Zp[q + oQ][p + oP] = buffer[index];
                }
            }
        }
    }
    return Z;
}
void FastDFJK::build_atom_pairs()
{
    atom_pairs_.clear();

    int natom = auxiliary_->molecule()->natom();
    for (int A = 0; A < natom; A++) {
        int nA = primary_->nshell_on_center(A);
        int oA = primary_->shell_on_center(A,0);
        for (int B = 0; B <= A; B++) {
            int nB = primary_->nshell_on_center(B);
            int oB = primary_->shell_on_center(B,0);
            for (int P = 0; P < nA; P++) {
                bool found = false;
                for (int Q = 0; Q < nB; Q++) {
                    if (sieve_->shell_pair_significant(P + oA, Q + oB)) {
                        found = true;
                        atom_pairs_.push_back(std::pair<int,int>(A,B));
                        break;
                    }
                }
                if (found) break;
            }
        }
    }
}
void FastDFJK::build_shell_pairs()
{
    shell_pairs_.clear();

    for (size_t AB = 0L; AB < atom_pairs_.size(); AB++) {
        int A = atom_pairs_[AB].first;
        int B = atom_pairs_[AB].second;
        std::vector<std::pair<int,int> > pairs;
        int nA = primary_->nshell_on_center(A);
        int oA = primary_->shell_on_center(A,0);
        int nB = primary_->nshell_on_center(B);
        int oB = primary_->shell_on_center(B,0);
        for (int P = 0; P < nA; P++) {
            for (int Q = 0; Q < nB; Q++) {
                if (Q + oB > P + oA) continue;
                if (sieve_->shell_pair_significant(P + oA, Q + oB)) {
                    pairs.push_back(std::pair<int,int>(P + oA, Q + oB));
                }
            }
        }
        shell_pairs_.push_back(pairs);
    }
}
void FastDFJK::build_auxiliary_partition()
{
    auxiliary_atoms_.clear();
    bump_atoms_.clear();

    // Atomic distance matrix
    std::shared_ptr<Molecule> mol = auxiliary_->molecule();
    int natom = mol->natom();
    std::shared_ptr<Matrix> RAB(new Matrix("RAB",natom,natom));
    double** RABp = RAB->pointer();
    for (int A = 0; A < natom; A++) {
        for (int B = A; B < natom; B++) {
            RABp[A][B] = RABp[B][A] =
                mol->xyz(A).distance(mol->xyz(B));
        }
    }
    double R0 = bump_R0_;
    double R1 = bump_R1_;

    if (domains_ == "DIATOMIC") {

        for (size_t pair = 0L; pair < atom_pairs_.size(); pair++) {
            int A = atom_pairs_[pair].first;
            int B = atom_pairs_[pair].second;
            std::vector<int> task;
            std::vector<double> bump;
            if (A == B) {
                task.push_back(A);
                bump.push_back(1.0);
            } else {
                task.push_back(A);
                task.push_back(B);
                bump.push_back(1.0);
                bump.push_back(1.0);
            }
            auxiliary_atoms_.push_back(task);
            bump_atoms_.push_back(bump);
        }

    } else if (domains_ == "SPHERES") {

        // Significant atom pairs list
        std::vector<std::vector<int> > sig;
        for (int A = 0; A < natom; A++) {
            std::vector<int> task;
            for (int B = 0; B < natom; B++) {
                if (RABp[A][B] < R1) {
                    task.push_back(B);
                }
            }
            sig.push_back(task);
        }

        for (size_t pair = 0L; pair < atom_pairs_.size(); pair++) {

            int A = atom_pairs_[pair].first;
            int B = atom_pairs_[pair].second;
            std::vector<int> task;
            std::vector<double> bump;

            // Diatomic frame
            if (A == B) {
                task.push_back(A);
                bump.push_back(1.0);
            } else {
                task.push_back(A);
                task.push_back(B);
                bump.push_back(1.0);
                bump.push_back(1.0);
            }

            // Unique augmentation atoms
            std::set<int> aug;
            const std::vector<int>& sigA = sig[A];
            const std::vector<int>& sigB = sig[B];
            for (size_t C2 = 0; C2 < sigA.size(); C2++) {
                int C = sigA[C2];
                if (! (C == A || C == B)) {
                    aug.insert(C);
                }
            }
            for (size_t C2 = 0; C2 < sigB.size(); C2++) {
                int C = sigB[C2];
                if (! (C == A || C == B)) {
                    aug.insert(C);
                }
            }

            // Augmentation and bump computation
            for (std::set<int>::const_iterator it = aug.begin();
                it != aug.end(); ++it) {
                int C = (*it);
                double bAC = 0.0;
                if (RABp[A][C] <= R0) {
                    bAC = 1.0;
                } else {
                    double R = RABp[A][C];
                    bAC = 1.0 / (1.0 + exp((R1 - R0) / (R1 - R) - (R1 - R0) / (R - R0)));
                }
                double bBC = 0.0;
                if (RABp[B][C] <= R0) {
                    bBC = 1.0;
                } else {
                    double R = RABp[B][C];
                    bBC = 1.0 / (1.0 + exp((R1 - R0) / (R1 - R) - (R1 - R0) / (R - R0)));
                }
                double b = sqrt(bAC * bBC); // Looks more like MHG in an atom-pair frame
                task.push_back(C);
                bump.push_back(b);
            }

            auxiliary_atoms_.push_back(task);
            bump_atoms_.push_back(bump);
        }

    } else {
        throw PSIEXCEPTION("Unrecognized fitting domain algorithm");
    }
}
void FastDFJK::build_Bpq()
{
    Bpq_.clear();
    Bpq_.resize(atom_pairs_.size());

    int nthread = 1;
    #ifdef _OPENMP
        nthread = df_ints_num_threads_;
    #endif

    std::shared_ptr<IntegralFactory>  fact1(new IntegralFactory(auxiliary_,BasisSet::zero_ao_basis_set(),primary_,primary_));
    std::shared_ptr<IntegralFactory> Jfact1(new IntegralFactory(auxiliary_,BasisSet::zero_ao_basis_set(),auxiliary_,BasisSet::zero_ao_basis_set()));
    std::vector<std::shared_ptr<TwoBodyAOInt> >  ints1;
    std::vector<std::shared_ptr<TwoBodyAOInt> > Jints1;

    std::shared_ptr<IntegralFactory>  fact2(new IntegralFactory(auxiliary_,primary_,primary_,primary_));
    std::shared_ptr<IntegralFactory> Jfact2(new IntegralFactory(auxiliary_));
    std::vector<std::shared_ptr<ThreeCenterOverlapInt> >  ints2;
    std::vector<std::shared_ptr<OneBodyAOInt> >          Jints2;

    bool algorithm;
    for (int thread = 0; thread < nthread; thread++) {
        if (metric_ == "COULOMB") {
            ints1.push_back(std::shared_ptr<TwoBodyAOInt>(fact1->eri()));
            Jints1.push_back(std::shared_ptr<TwoBodyAOInt>(Jfact1->eri()));
            algorithm = false;
        } else if (metric_ == "EWALD") {
            ints1.push_back(std::shared_ptr<TwoBodyAOInt>(fact1->erf_complement_eri(theta_)));
            Jints1.push_back(std::shared_ptr<TwoBodyAOInt>(Jfact1->erf_complement_eri(theta_)));
            algorithm = false;
        } else if (metric_ == "OVERLAP") {
            ints2.push_back(std::shared_ptr<ThreeCenterOverlapInt>(fact2->overlap_3c()));
            Jints2.push_back(std::shared_ptr<OneBodyAOInt>(Jfact2->ao_overlap()));
            algorithm = true;
        } else {
            throw PSIEXCEPTION("Unknown fitting metric");
        }
    }

    #pragma omp parallel for schedule(dynamic) num_threads(nthread)
    for (size_t pair = 0L; pair < atom_pairs_.size(); pair++) {

        int thread = 0;
        #ifdef _OPENMP
            thread = omp_get_thread_num();
        #endif
        const double* buffer  = (algorithm ?  ints2[thread]->buffer() :  ints1[thread]->buffer());
        const double* Jbuffer = (algorithm ? Jints2[thread]->buffer() : Jints1[thread]->buffer());

        const std::vector<std::pair<int,int> >& shell_pairs = shell_pairs_[pair];
        const std::vector<int>& auxiliary_atoms = auxiliary_atoms_[pair];
        const std::vector<double>& bump_atoms = bump_atoms_[pair];

        // => Sizing <= //

        int npq = 0;
        for (size_t PQ = 0; PQ < shell_pairs.size(); PQ++) {
            int P = shell_pairs[PQ].first;
            int Q = shell_pairs[PQ].second;
            npq += primary_->shell(P).nfunction() * primary_->shell(Q).nfunction();;
        }

        int naux = 0;
        for (size_t C = 0; C < auxiliary_atoms.size(); C++) {
            int C2 = auxiliary_atoms[C];
            int nC = auxiliary_->nshell_on_center(C2);
            int oC = auxiliary_->shell_on_center(C2,0);
            for (int A = 0; A < nC; A++) {
                naux += auxiliary_->shell(A + oC).nfunction();
            }
        }

        // => Tensor Allocation <= //

        std::shared_ptr<Matrix> Apq(new Matrix("Apq", naux, npq));
        double** Ap = Apq->pointer();
        std::shared_ptr<Matrix> Bpq(new Matrix("Bpq", naux, npq));
        double** Bp = Bpq->pointer();
        std::shared_ptr<Matrix> J(new Matrix("J", naux, naux));
        double** Jp = J->pointer();

        // => Generate Integrals <= //

        for (size_t C = 0,dA=0; C < auxiliary_atoms.size(); C++) {
            int C2 = auxiliary_atoms[C];
            int nC = auxiliary_->nshell_on_center(C2);
            int oC = auxiliary_->shell_on_center(C2,0);
            for (int A = oC; A < oC + nC; A++) {
                int nA = auxiliary_->shell(A).nfunction();
                for (size_t PQ = 0, dPQ = 0; PQ < shell_pairs.size(); PQ++) {
                    int P = shell_pairs[PQ].first;
                    int Q = shell_pairs[PQ].second;
                    if (algorithm) {
                        ints2[thread]->compute_shell(A,P,Q);
                    } else {
                        ints1[thread]->compute_shell(A,0,P,Q);
                    }
                    int nP = primary_->shell(P).nfunction();
                    int nQ = primary_->shell(Q).nfunction();
                    for (int a = 0, index = 0; a < nA; a++) {
                        for (int p = 0; p < nP; p++) {
                            for (int q = 0; q < nQ; q++,index++) {
                                Ap[a + dA][p * nQ + q + dPQ] = buffer[index];
                            }
                        }
                    }
                    dPQ += nP * nQ;
                }
                dA += nA;
            }
        }

        // => Generate Metric <= //

        for (size_t C = 0, dA = 0; C < auxiliary_atoms.size(); C++) {
            int C2 = auxiliary_atoms[C];
            int nC = auxiliary_->nshell_on_center(C2);
            int oC = auxiliary_->shell_on_center(C2,0);
            for (int A = oC; A < oC + nC; A++) {
                int nA = auxiliary_->shell(A).nfunction();
                for (size_t D = 0, dB = 0; D < auxiliary_atoms.size(); D++) {
                    int D2 = auxiliary_atoms[D];
                    int nD = auxiliary_->nshell_on_center(D2);
                    int oD = auxiliary_->shell_on_center(D2,0);
                    for (int B = oD; B < oD + nD; B++) {
                        int nB = auxiliary_->shell(B).nfunction();
                        if (B > A) { dB += nB; continue; }
                        if (algorithm) {
                            Jints2[thread]->compute_shell(A,B);
                        } else {
                            Jints1[thread]->compute_shell(A,0,B,0);
                        }
                        for (int a = 0, index = 0; a < nA; a++) {
                            for (int b = 0; b < nB; b++) {
                                Jp[a + dA][b + dB] = Jp[b + dB][a + dA] = Jbuffer[index++];
                            }
                        }
                        dB += nB;
                    }
                }
                dA += nA;
            }
        }

        // => "Bump" the metric <= //

        bump(J,bump_atoms,auxiliary_atoms, false);

        // => Invert Metric <= //

        J->power(-1.0,condition_);

        // => "Bump" the inverse metric <= //

        bump(J,bump_atoms,auxiliary_atoms, true);

        // => Apply Metric <= //

        C_DGEMM('N','N',naux,npq,naux,1.0,Jp[0],naux,Ap[0],npq,0.0,Bp[0],npq);

        Bpq_[pair] = Bpq;
    }
}
void FastDFJK::bump(std::shared_ptr<Matrix> J, const std::vector<double>& bump_atoms, const std::vector<int>& auxiliary_atoms, bool bump_diagonal)
{
    double** Jp = J->pointer();
    for (size_t C = 0, dA = 0; C < auxiliary_atoms.size(); C++) {
        int C2 = auxiliary_atoms[C];
        int nC = auxiliary_->nshell_on_center(C2);
        int oC = auxiliary_->shell_on_center(C2,0);
        for (int A = oC; A < oC + nC; A++) {
            int nA = auxiliary_->shell(A).nfunction();
            for (size_t D = 0, dB = 0; D < auxiliary_atoms.size(); D++) {

                double scale = ((!bump_diagonal && C == D) ? 1.0 : bump_atoms[C] * bump_atoms[D]);

                int D2 = auxiliary_atoms[D];
                int nD = auxiliary_->nshell_on_center(D2);
                int oD = auxiliary_->shell_on_center(D2,0);
                for (int B = oD; B < oD + nD; B++) {
                    int nB = auxiliary_->shell(B).nfunction();
                    for (int a = 0; a < nA; a++) {
                        for (int b = 0; b < nB; b++) {
                            Jp[a + dA][b + dB] *= scale;
                        }
                    }
                    dB += nB;
                }
            }
            dA += nA;
        }
    }
}

void FastDFJK::build_J(std::shared_ptr<Matrix> Z,
                       const std::vector<std::shared_ptr<Matrix> >& D,
                       const std::vector<std::shared_ptr<Matrix> >& J)
{
    // => Sizing <= //

    int natom = auxiliary_->molecule()->natom();
    int nso  = primary_->nbf();
    int naux = auxiliary_->nbf();

    int max_naux_per_atom = 0;
    int max_nso_per_atom = 0;
    for (int A = 0; A < natom; A++) {
        int P1 = primary_->shell_on_center(A,0);
        int P2 = (A == natom - 1 ? primary_->nshell() : primary_->shell_on_center(A+1,0));
        int nP = (A == natom - 1 ? nso - primary_->shell(P1).function_index() : primary_->shell(P2).function_index() - primary_->shell(P1).function_index());
        max_nso_per_atom = (max_nso_per_atom >= nP ? max_nso_per_atom : nP);
        int A1 = auxiliary_->shell_on_center(A,0);
        int A2 = (A == natom - 1 ? auxiliary_->nshell() : auxiliary_->shell_on_center(A+1,0));
        int nA = (A == natom - 1 ? naux - auxiliary_->shell(A1).function_index() : auxiliary_->shell(A2).function_index() - auxiliary_->shell(A1).function_index());
        max_naux_per_atom = (max_naux_per_atom >= nA ? max_naux_per_atom : nA);
    }
    size_t max_atom = 0;
    for (size_t pair = 0L; pair < auxiliary_atoms_.size(); pair++) {
        max_atom = (max_atom >= auxiliary_atoms_[pair].size() ? max_atom : auxiliary_atoms_[pair].size());
    }
    int max_nso2 = max_nso_per_atom * max_nso_per_atom;
    int max_naux = max_atom * max_naux_per_atom;

    int nthread = 1;
    #ifdef _OPENMP
        nthread = Process::environment.get_n_threads();
    #endif

    // => Temporaries <= //

    std::vector<std::shared_ptr<Vector> > vs;
    std::vector<std::shared_ptr<Vector> > cs;
    std::vector<std::shared_ptr<Vector> > c;
    for (int thread = 0; thread < nthread; thread++) {
        vs.push_back(std::shared_ptr<Vector>(new Vector("vs", max_nso2)));
        cs.push_back(std::shared_ptr<Vector>(new Vector("cs", max_naux)));
        c.push_back( std::shared_ptr<Vector>(new Vector("c", naux)));
    }
    std::shared_ptr<Vector> d(new Vector("d",naux));

    // ==> Master Loop over J Tasks <= //

    for (size_t ind = 0; ind < D.size(); ind++) {

        double** Zp = Z->pointer();
        double** Dp = D[ind]->pointer();
        double** Jp = J[ind]->pointer();

        // => B_rs^A D_rs -> c_A <= //

        for (int thread = 0; thread < nthread; thread++) {
            c[thread]->zero();
        }

        #pragma omp parallel for num_threads(nthread) schedule(dynamic)
        for (size_t pair = 0L; pair < atom_pairs_.size(); pair++) {

            int thread = 0;
            #ifdef _OPENMP
                thread = omp_get_thread_num();
            #endif

            const std::vector<std::pair<int,int> >& shell_pairs = shell_pairs_[pair];
            const std::vector<int>& auxiliary_atoms = auxiliary_atoms_[pair];
            std::shared_ptr<Matrix> B = Bpq_[pair];

            int nauxs = B->rowspi()[0];
            int nso2s = B->colspi()[0];

            double** Bp = B->pointer();
            double* vsp = vs[thread]->pointer();
            double* csp = cs[thread]->pointer();
            double* cp  = c[thread]->pointer();

            for (size_t PQ = 0, dPQ = 0; PQ < shell_pairs.size(); PQ++) {
                int P = shell_pairs[PQ].first;
                int Q = shell_pairs[PQ].second;
                int nP = primary_->shell(P).nfunction();
                int nQ = primary_->shell(Q).nfunction();
                int oP = primary_->shell(P).function_index();
                int oQ = primary_->shell(Q).function_index();
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        vsp[p * nQ + q + dPQ] = Dp[p + oP][q + oQ];
                        if (P != Q) {
                            vsp[p * nQ + q + dPQ] += Dp[q + oQ][p + oP];
                        }
                    }
                }
                dPQ += nP * nQ;
            }

            C_DGEMV('N',nauxs,nso2s,1.0,Bp[0],nso2s,vsp,1,0.0,csp,1);

            for (size_t C2 = 0, dA = 0; C2 < auxiliary_atoms.size(); C2++) {
                int C = auxiliary_atoms[C2];
                int nC = auxiliary_->nshell_on_center(C);
                int oC = auxiliary_->shell_on_center(C,0);
                for (int A = 0; A < nC; A++) {
                    int nA = auxiliary_->shell(A + oC).nfunction();
                    int oA = auxiliary_->shell(A + oC).function_index();
                    for (int a = 0; a < nA; a++) {
                        cp[a + oA] += csp[a + dA];
                    }
                    dA += nA;
                }
            }
        }

        for (int thread = 1; thread < nthread; thread++) {
            c[0]->add(c[thread]);
        }

        // => Z^BA c_A -> d_B <= //

        double* c2p = c[0]->pointer();
        double* dp = d->pointer();
        C_DGEMV('N',naux,naux,1.0,Zp[0],naux,c2p,1,0.0,dp,1);

        // => B_pq^B d_B -> J_pq <= //

        #pragma omp parallel for num_threads(nthread) schedule(dynamic)
        for (size_t pair = 0L; pair < atom_pairs_.size(); pair++) {

            int thread = 0;
            #ifdef _OPENMP
                thread = omp_get_thread_num();
            #endif

            const std::vector<std::pair<int,int> >& shell_pairs = shell_pairs_[pair];
            const std::vector<int>& auxiliary_atoms = auxiliary_atoms_[pair];
            std::shared_ptr<Matrix> B = Bpq_[pair];

            int nauxs = B->rowspi()[0];
            int nso2s = B->colspi()[0];

            double** Bp = B->pointer();
            double* vsp = vs[thread]->pointer();
            double* csp = cs[thread]->pointer();

            for (size_t C2 = 0, dA = 0; C2 < auxiliary_atoms.size(); C2++) {
                int C = auxiliary_atoms[C2];
                int nC = auxiliary_->nshell_on_center(C);
                int oC = auxiliary_->shell_on_center(C,0);
                for (int A = 0; A < nC; A++) {
                    int nA = auxiliary_->shell(A + oC).nfunction();
                    int oA = auxiliary_->shell(A + oC).function_index();
                    for (int a = 0; a < nA; a++) {
                        csp[a + dA] = dp[a + oA];
                    }
                    dA += nA;
                }
            }

            C_DGEMV('T',nauxs,nso2s,1.0,Bp[0],nso2s,csp,1,0.0,vsp,1);

            for (size_t PQ = 0, dPQ = 0; PQ < shell_pairs.size(); PQ++) {
                int P = shell_pairs[PQ].first;
                int Q = shell_pairs[PQ].second;
                int nP = primary_->shell(P).nfunction();
                int nQ = primary_->shell(Q).nfunction();
                int oP = primary_->shell(P).function_index();
                int oQ = primary_->shell(Q).function_index();
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        Jp[p + oP][q + oQ] = vsp[p * nQ + q + dPQ];
                        if (P != Q) {
                            Jp[q + oQ][p + oP] = vsp[p * nQ + q + dPQ];
                        }
                    }
                }
                dPQ += nP * nQ;
            }
        }
    }
}
void FastDFJK::build_K(std::shared_ptr<Matrix> /*Z*/,
                       const std::vector<std::shared_ptr<Matrix> >& /*D*/,
                       const std::vector<std::shared_ptr<Matrix> >& /*K*/)
{
    throw PSIEXCEPTION("K: Not implemented yet");
}
}
