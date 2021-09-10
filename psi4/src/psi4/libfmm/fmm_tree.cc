#include "psi4/pragma.h"

#include "psi4/libfmm/multipoles_helper.h"
#include "psi4/libfmm/fmm_tree.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector3.h"
#include "psi4/libmints/gshell.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/onebody.h"
#include "psi4/libmints/multipoles.h"
#include "psi4/libmints/overlap.h"
#include "psi4/libmints/twobody.h"
#include "psi4/libqt/qt.h"
#include "psi4/libpsi4util/process.h"

#include <functional>
#include <memory>
#include <tuple>
#include <vector>
#include <cmath>
#include <unordered_map>
#include <utility>
#include <csignal>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace psi {

int num_digits(long n) {
    if (n == 0) return 1;
    return (int) std::log10(std::abs(n)) + 1;
}

ShellPair::ShellPair(std::shared_ptr<BasisSet>& basisset, std::pair<int, int> pair_index, 
                     std::shared_ptr<HarmonicCoefficients>& mpole_coefs) {     
    basisset_ = basisset;
    pair_index_ = pair_index;

    const GaussianShell& Pshell = basisset_->shell(pair_index.first);
    const GaussianShell& Qshell = basisset_->shell(pair_index.second);

    Vector3 pcenter = Pshell.center();
    Vector3 qcenter = Qshell.center();

    center_ = Vector3(0.0, 0.0, 0.0);
    exp_ = INFINITY;

    int nprim_p = Pshell.nprimitive();
    int nprim_q = Qshell.nprimitive();
    for (int pp = 0; pp < nprim_p; pp++) {
        double pcoef = Pshell.coef(pp);
        double pexp = Pshell.exp(pp);
        for (int qp = 0; qp < nprim_q; qp++) {
            double qcoef = Qshell.coef(qp);
            double qexp = Qshell.exp(qp);

            const double pq_exp = std::abs(pexp + qexp);
            Vector3 pq_center = (pexp * pcenter + qexp * qcenter) / pq_exp;

            center_ += pq_center;
            exp_ = std::min(exp_, pq_exp);
        }
    }
    center_ /= (nprim_p * nprim_q);
    extent_ = ERFCI10 * std::sqrt(2.0 / exp_);

    mpole_coefs_ = mpole_coefs;
}

void ShellPair::calculate_mpoles(Vector3 box_center, std::shared_ptr<OneBodyAOInt> s_ints,
                            std::shared_ptr<OneBodyAOInt> mpole_ints, int lmax) {
    
    // number of total multipoles to compute, -1 since the overlap is not computed
    int nmpoles = (lmax + 1) * (lmax + 2) * (lmax + 3) / 6 - 1;

    int P = pair_index_.first;
    int Q = pair_index_.second;

    // Calculate the overlap integrals (Order 0 multipole integrals)
    s_ints->compute_shell(P, Q);
    const double* sbuffer = s_ints->buffer();

    // Calculate the multipole integrals
    mpole_ints->set_origin(box_center);
    mpole_ints->compute_shell(P, Q);
    const double* mbuffer = mpole_ints->buffer();

    const GaussianShell& Pshell = basisset_->shell(P);
    const GaussianShell& Qshell = basisset_->shell(Q);

    int p_start = Pshell.start();
    int num_p = Pshell.nfunction();

    int q_start = Qshell.start();
    int num_q = Qshell.nfunction();

    for (int p = p_start; p < p_start + num_p; p++) {
        int dp = p - p_start;
        for (int q = q_start; q < q_start + num_q; q++) {
            int dq = q - q_start;

            std::shared_ptr<RealSolidHarmonics> pq_mpoles = std::make_shared<RealSolidHarmonics>(lmax, box_center, Regular);

            pq_mpoles->add(0, 0, sbuffer[dp * num_q + dq]);

            int running_index = 0;
            for (int l = 1; l <= lmax; l++) {
                int l_ncart = ncart(l);
                for (int m = -l; m <= l; m++) {
                    int mu = m_addr(m);
                    std::unordered_map<int, double>& mpole_terms= mpole_coefs_->get_terms(l, mu);

                    int powdex = 0;
                    for (int ii = 0; ii <= l; ii++) {
                        int a = l - ii;
                        for (int jj = 0; jj <= ii; jj++) {
                            int b = ii - jj;
                            int c = jj;
                            int ind = a * l_ncart * l_ncart + b * l_ncart + c;

                            if (mpole_terms.count(ind)) {
                                double coef = mpole_terms[ind];
                                int abcindex = powdex + running_index;
                                pq_mpoles->add(l, mu, pow(-1.0, (double) l+1) * coef * mbuffer[abcindex * num_p * num_q + dp * num_q + dq]);
                            }
                            powdex += 1;
                        } // end jj
                    } // end ii
                } // end m loop
                running_index += l_ncart;
            } // end l
            mpoles_.push_back(pq_mpoles);
        } // end q
    } // end p

}

CFMMBox::CFMMBox(std::shared_ptr<CFMMBox> parent, std::vector<std::shared_ptr<ShellPair>> shell_pairs, 
              Vector3 origin, double length, int level, int lmax, int ws) {
    parent_ = parent;

    shell_pairs_ = shell_pairs;
    origin_ = origin;
    center_ = origin_ + 0.5 * Vector3(length, length, length);
    length_ = length;
    level_ = level;
    lmax_ = lmax;
    ws_ = ws;

    mpoles_ = std::make_shared<RealSolidHarmonics>(lmax, center_, Regular);
    Vff_ = std::make_shared<RealSolidHarmonics>(lmax, center_, Irregular);

    nthread_ = 1;
#ifdef _OPENMP
    nthread_ = Process::environment.get_n_threads();
#endif

}

void CFMMBox::make_children() {

    int nchild = (level_ > 0) ? 16 : 8;
    std::vector<std::vector<std::shared_ptr<ShellPair>>> child_shell_pair_buffer(nchild);

    // Fill order (ws,z,y,x) (0)000 (0)001 (0)010 (0)011 (0)100 (0)101 (0)110 (0)111
    // (1)000 (1)001 (1)010 (1)011 (1)100 (1)101 (1)110 (1)111
    for (std::shared_ptr<ShellPair> shell_pair : shell_pairs_) {
        Vector3 sp_center = shell_pair->get_center();
        double x = sp_center[0];
        double y = sp_center[1];
        double z = sp_center[2];
        double extent = shell_pair->get_extent();
        int ws = std::max(2, 2 * (int)std::ceil(extent / length_));

        int xbit = (x < center_[0]) ? 0 : 1;
        int ybit = (y < center_[1]) ? 0 : 1;
        int zbit = (z < center_[2]) ? 0 : 1;
        int rbit = (level_ == 0 || ws < 2 * ws_) ? 0 : 1;

        int boxind = 8 * rbit + 4 * zbit + 2 * ybit + 1 * xbit;
        child_shell_pair_buffer[boxind].push_back(shell_pair);
    }

    // Make the children
    for (int boxind = 0; boxind < nchild; boxind++) {
        int xbit = boxind % 2;
        int ybit = (boxind / 2) % 2;
        int zbit = (boxind / 4) % 2;
        int rbit = (boxind / 8) % 2;
        Vector3 new_origin = origin_ + Vector3(xbit * 0.5 * length_, ybit * 0.5 * length_, zbit * 0.5 * length_);
        int child_ws = 2 * ws_ - 2 + 2 * rbit;
        children_.push_back(std::make_shared<CFMMBox>(this->get(), child_shell_pair_buffer[boxind], new_origin, 
                                                          0.5 * length_, level_ + 1, lmax_, child_ws));
    }
}

void CFMMBox::compute_far_field() {

    // Creates a temporary parent shared pointer
    std::shared_ptr<CFMMBox> parent = parent_.lock();

    // Parent is not a nullpointer
    if (parent) {
        // Siblings of this box (Technically near fields include self in this implementation)
        for (std::shared_ptr<CFMMBox> sibling : parent->children_) {
            Vector3 Rab = center_ - sibling->center_;
            double rab = std::sqrt(Rab.dot(Rab));

            int ref_ws = (ws_ + sibling->ws_) / 2;
            if (rab <= ref_ws * length_) {
                near_field_.push_back(sibling);
            } else {
                local_far_field_.push_back(sibling);
            }
        }

        // Parent's near field (Cousins)
        for (std::shared_ptr<CFMMBox> uncle : parent->near_field_) {
            if (uncle.get() == parent.get()) continue;
            for (std::shared_ptr<CFMMBox> cousin : uncle->children_) {
                Vector3 Rab = center_ - cousin->center_;
                double rab = std::sqrt(Rab.dot(Rab));

                int ref_ws = (ws_ + cousin->ws_) / 2;
                if (rab <= ref_ws * length_) {
                    near_field_.push_back(cousin);
                } else {
                    local_far_field_.push_back(cousin);
                }
            }
        }

        // Calculate the far field vector
        for (std::shared_ptr<CFMMBox> box : local_far_field_) {
            // The far field effect the boxes have on this particular box
            std::shared_ptr<RealSolidHarmonics> far_field = box->mpoles_->far_field_vector(center_);
            Vff_->add(far_field);
        }

        // Add the parent's far field contribution
        Vff_->add(parent->Vff_->translate(center_));
    }

    else {
        near_field_.push_back(this->get());
    }

}

void CFMMBox::compute_mpoles(std::shared_ptr<BasisSet>& basisset, std::vector<SharedMatrix>& D) {

    std::shared_ptr<IntegralFactory> int_factory = std::make_shared<IntegralFactory>(basisset);

    std::vector<std::shared_ptr<OneBodyAOInt>> mpints(nthread_);
    std::vector<std::shared_ptr<OneBodyAOInt>> sints(nthread_);

    for (int thread = 0; thread < nthread_; thread++) {
        mpints[thread] = std::shared_ptr<OneBodyAOInt>(int_factory->ao_multipoles(lmax_));
        sints[thread] = std::shared_ptr<OneBodyAOInt>(int_factory->ao_overlap());
        mpints[thread]->set_origin(center_);
    }

    // Compute multipoles for all function pairs in the shell pair
#pragma omp parallel for
    for (int ind = 0; ind < shell_pairs_.size(); ind++) {
        std::shared_ptr<ShellPair> sp = shell_pairs_[ind];
        
        int thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif
        sp->calculate_mpoles(center_, sints[thread], mpints[thread], lmax_);
    }

    // Contract the multipoles with the density matrix to get box multipoles
    for (const auto& sp : shell_pairs_) {
        std::vector<std::shared_ptr<RealSolidHarmonics>>& sp_mpoles = sp->get_mpoles();

        std::pair<int, int> PQ = sp->get_shell_pair_index();
        int P = PQ.first;
        int Q = PQ.second;

        double prefactor = (P == Q) ? 2.0 : 4.0;

        const GaussianShell& Pshell = basisset->shell(P);
        const GaussianShell& Qshell = basisset->shell(Q);

        int p_start = Pshell.start();
        int num_p = Pshell.nfunction();

        int q_start = Qshell.start();
        int num_q = Qshell.nfunction();

        for (int N = 0; N < D.size(); N++) {
            for (int p = p_start; p < p_start + num_p; p++) {
                int dp = p - p_start;
                for (int q = q_start; q < q_start + num_q; q++) {
                    int dq = q - q_start;
                    std::shared_ptr<RealSolidHarmonics> basis_mpole = sp_mpoles[dp * num_q + dq]->copy();
                    
                    basis_mpole->scale(prefactor * D[N]->get(p, q));
                    mpoles_->add(basis_mpole);
                } // end q
            } // end p
        } // end N
    }

}

void CFMMBox::compute_mpoles_from_children() {

    for (std::shared_ptr<CFMMBox> child : children_) {
        std::shared_ptr<RealSolidHarmonics> child_mpoles = child->mpoles_->translate(center_);
        mpoles_->add(child_mpoles);
    }

}

CFMMTree::CFMMTree(std::shared_ptr<Molecule> molecule, std::shared_ptr<BasisSet> basisset, std::vector<SharedMatrix>& D, 
                std::vector<SharedMatrix>& J, const std::vector<std::pair<int, int>>& shell_pairs, int nlevels, int lmax) {
    molecule_ = molecule;
    basisset_ = basisset;
    D_ = D;
    J_ = J;
    nlevels_ = nlevels;
    lmax_ = lmax;
    int num_boxes = (nlevels_ == 1) ? 1 : (0.5 * std::pow(16, nlevels_) + 7) / 15;
    tree_.resize(num_boxes);
    Options& options = Process::environment.options;
    mpole_coefs_ = std::make_shared<HarmonicCoefficients>(options.get_int("CFMM_MAX_MPOLE_ORDER"), Regular);

    std::shared_ptr<IntegralFactory> factory = std::make_shared<IntegralFactory>(basisset_);
    std::shared_ptr<TwoBodyAOInt> eri = std::shared_ptr<TwoBodyAOInt>(factory->eri());
    ints_.push_back(eri);

    nthread_ = 1;
#ifdef _OPENMP
    nthread_ = Process::environment.get_n_threads();
#endif

    for (int thread = 1; thread < nthread_; thread++) {
        ints_.push_back(std::shared_ptr<TwoBodyAOInt>(eri->clone()));
    }

    for (const auto& pair : shell_pairs) {
        shell_pairs_.push_back(std::make_shared<ShellPair>(basisset_, pair, mpole_coefs_));
    }

    // sort_shell_pairs();
    make_root_node();
    make_children();

    int print = options.get_int("PRINT");
    if (print >= 2) print_out();
}

void CFMMTree::sort_shell_pairs() {

    // Number of digits of each phase of the sort
    // Resolution for floats is 0.001 bohr
    // Sort by z, y, x, and then radial extents
    int zdig = 1, ydig = 1, xdig = 1, rdig = 1;

    for (const auto& shell_pair : shell_pairs_) {
        Vector3 center = shell_pair->get_center();
        double x = center[0];
        double y = center[1];
        double z = center[2];
        double r = shell_pair->get_extent();

        zdig = std::max(zdig, num_digits(z * 1000));
        ydig = std::max(ydig, num_digits(y * 1000));
        xdig = std::max(xdig, num_digits(x * 1000));
        rdig = std::max(rdig, num_digits(r * 1000));
    }

    std::vector<std::vector<std::shared_ptr<ShellPair>>> buckets(19);

    // Sort by x coordinates
    long curr10 = 10;

    for (int iter = 0; iter < xdig; iter++) {
        for (int ind = 0; ind < shell_pairs_.size(); ind++) {
            int dig = ((long)(shell_pairs_[ind]->get_center()[0] * 1000) / curr10) % 10;
            buckets[dig + 9].push_back(shell_pairs_[ind]);
        }
        shell_pairs_.clear();
        for (int b = 0; b < 19; b++) {
            int nelem = buckets[b].size();
            for (int idx = 0; idx < nelem; idx++) {
                shell_pairs_.push_back(buckets[b].back());
                buckets[b].pop_back();
            }
        }
        curr10 *= 10;
    }

    // Sort by y coordinates
    curr10 = 10;

    for (int iter = 0; iter < ydig; iter++) {
        for (int ind = 0; ind < shell_pairs_.size(); ind++) {
            int dig = ((long)(shell_pairs_[ind]->get_center()[1] * 1000) / curr10) % 10;
            buckets[dig + 9].push_back(shell_pairs_[ind]);
        }
        shell_pairs_.clear();
        for (int b = 0; b < 19; b++) {
            int nelem = buckets[b].size();
            for (int idx = 0; idx < nelem; idx++) {
                shell_pairs_.push_back(buckets[b].back());
                buckets[b].pop_back();
            }
        }
        curr10 *= 10;
    }

    // Sort by z coordinates
    curr10 = 10;

    for (int iter = 0; iter < zdig; iter++) {
        for (int ind = 0; ind < shell_pairs_.size(); ind++) {
            int dig = ((long)(shell_pairs_[ind]->get_center()[2] * 1000) / curr10) % 10;
            buckets[dig + 9].push_back(shell_pairs_[ind]);
        }
        shell_pairs_.clear();
        for (int b = 0; b < 19; b++) {
            int nelem = buckets[b].size();
            for (int idx = 0; idx < nelem; idx++) {
                shell_pairs_.push_back(buckets[b].back());
                buckets[b].pop_back();
            }
        }
        curr10 *= 10;
    }

    // Sort by radial extents
    curr10 = 10;
    for (int iter = 0; iter < rdig; iter++) {
        for (int ind = 0; ind < shell_pairs_.size(); ind++) {
            int dig = ((long)(shell_pairs_[ind]->get_extent() * 1000) / curr10) % 10;
            buckets[dig + 9].push_back(shell_pairs_[ind]);
        }
        shell_pairs_.clear();
        for (int b = 0; b < 19; b++) {
            int nelem = buckets[b].size();
            for (int idx = 0; idx < nelem; idx++) {
                shell_pairs_.push_back(buckets[b].back());
                buckets[b].pop_back();
            }
        }
        curr10 *= 10;
    }
}

void CFMMTree::make_root_node() {
    double min_dim = molecule_->x(0);
    double max_dim = molecule_->x(0);

    for (int atom = 0; atom < molecule_->natom(); atom++) {
        double x = molecule_->x(atom);
        double y = molecule_->y(atom);
        double z = molecule_->z(atom);
        min_dim = std::min(x, min_dim);
        min_dim = std::min(y, min_dim);
        min_dim = std::min(z, min_dim);
        max_dim = std::max(x, max_dim);
        max_dim = std::max(y, max_dim);
        max_dim = std::max(z, max_dim);
    }

    max_dim += 0.1; // Add a small buffer to the box

    Vector3 origin = Vector3(min_dim, min_dim, min_dim);
    double length = (max_dim - min_dim);

    tree_[0] = std::make_shared<CFMMBox>(nullptr, shell_pairs_, origin, length, 0, lmax_, 2);
}

void CFMMTree::make_children() {

    timer_on("CFMMTree::make_children");

    for (int level = 0; level <= nlevels_ - 2; level += 1) {
        int start, end;
        if (level == 0) {
            start = 0;
            end = 1;
        } else {
            start = (0.5 * std::pow(16, level) + 7) / 15;
            end = (0.5 * std::pow(16, level+1) + 7) / 15;
        }

#pragma omp parallel for
        for (int bi = start; bi < end; bi++) {
            tree_[bi]->make_children();
            auto children = tree_[bi]->get_children();

            for (int ci = 0; ci < children.size(); ci++) {
                int ti = (level == 0) ? ci + 1 : bi * 16 - 7 + ci;
                tree_[ti] = children[ci];
            }
        }
    }

    timer_off("CFMMTree::make_children");
}

void CFMMTree::calculate_multipoles() {
    timer_on("CFMMTree::calculate_multipoles");

    // Compute bottom layer mpoles (parallel in call)
    int start, end;
    if (nlevels_ == 1) {
        start = 0;
        end = 1;
    } else {
        start = (0.5 * std::pow(16, nlevels_-1) + 7) / 15;
        end = (0.5 * std::pow(16, nlevels_) + 7) / 15;
    }

    for (int bi = start; bi < end; bi += 1) {
        if (tree_[bi]->get_nsp() == 0) continue;
        tree_[bi]->compute_mpoles(basisset_, D_);
    }

    // Translate parents (parallel here)
    for (int level = nlevels_ - 2; level >= 0; level -= 1) {
        int start, end;
        if (level == 0) {
            start = 0;
            end = 1;
        } else {
            start = (0.5 * std::pow(16, level) + 7) / 15;
            end = (0.5 * std::pow(16, level+1) + 7) / 15;
        }

#pragma omp parallel for
        for (int bi = start; bi < end; bi++) {
            if (tree_[bi]->get_nsp() == 0) continue;
            tree_[bi]->compute_mpoles_from_children();
        }
    }

    timer_off("CFMMTree::calculate_multipoles");
}

void CFMMTree::compute_far_field() {

    timer_on("CFMMTree::compute_far_field");

    for (int level = 0; level <= nlevels_ - 1; level += 1) {
        int start, end;
        if (level == 0) {
            start = 0;
            end = 1;
        } else {
            start = (0.5 * std::pow(16, level) + 7) / 15;
            end = (0.5 * std::pow(16, level+1) + 7) / 15;
        }

#pragma omp parallel for
        for (int bi = start; bi < end; bi++) {
            if (tree_[bi]->get_nsp() == 0) continue;
            tree_[bi]->compute_far_field();
        }
    }

    timer_off("CFMMTree::compute_far_field");

}

void CFMMTree::build_nf_J() {

    timer_on("CFMMTree::build_nf_J");

    int nshell = basisset_->nshell();

    // Starting and ending box indices
    int start = (nlevels_ == 1) ? 0 : (0.5 * std::pow(16, nlevels_-1) + 7) / 15;
    int end = (nlevels_ == 1) ? 1 : (0.5 * std::pow(16, nlevels_) + 7) / 15;
    
    // A map of the function (num_p * num_q) offsets per shell-pair
    std::unordered_map<int, int> offsets;

    // Maximum space (pfunc * qfunc) to allocate per box-task
    size_t max_alloc = 0;
    
    for (int bi = start; bi < end; bi++) {
        auto& PQshells = tree_[bi]->get_shell_pairs();
        int PQoff = 0;
        for (int PQind = 0; PQind < PQshells.size(); PQind++) {
            std::pair<int, int> PQ = PQshells[PQind]->get_shell_pair_index();
            int P = PQ.first;
            int Q = PQ.second;
            offsets[P * nshell + Q] = PQoff;
            int Pfunc = basisset_->shell(P).nfunction();
            int Qfunc = basisset_->shell(Q).nfunction();
            PQoff += Pfunc * Qfunc;
        }
        max_alloc = std::max((size_t)PQoff, max_alloc);
    }

    // Make intermediate buffers (for threading purposes and take advantage of 8-fold perm sym)
    std::vector<std::vector<std::vector<double>>> JT;
    for (int thread = 0; thread < nthread_; thread++) {
        std::vector<std::vector<double>> J2;
        for (size_t N = 0; N < D_.size(); N++) {
            std::vector<double> temp(2 * max_alloc);
            J2.push_back(temp);
        }
        JT.push_back(J2);
    }

#pragma omp parallel for
    for (int bi = start; bi < end; bi++) {
        std::shared_ptr<CFMMBox> curr = tree_[bi];
        auto& PQshells = curr->get_shell_pairs();
        auto& nf_boxes = curr->near_field_boxes();

        int thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif
        
        for (int nfi = 0; nfi < nf_boxes.size(); nfi++) {
            std::shared_ptr<CFMMBox> neighbor = nf_boxes[nfi];
            auto& RSshells = neighbor->get_shell_pairs();
            
            bool touched = false;
            for (int PQind = 0; PQind < PQshells.size(); PQind++) {
                std::pair<int, int> PQ = PQshells[PQind]->get_shell_pair_index();
                int P = PQ.first;
                int Q = PQ.second;
            
                for (int RSind = 0; RSind < RSshells.size(); RSind++) {
                    std::pair<int, int> RS = RSshells[RSind]->get_shell_pair_index();
                    int R = RS.first;
                    int S = RS.second;
                    
                    if (R * nshell + S > P * nshell + Q) continue;
                    if (!ints_[thread]->shell_significant(P, Q, R, S)) continue;

                    double prefactor = 1.0;
                    if (P != Q) prefactor *= 2;
                    if (R != S) prefactor *= 2;
                    if (P == R && Q == S) prefactor *= 0.5;

                    const GaussianShell& Pshell = basisset_->shell(P);
                    const GaussianShell& Qshell = basisset_->shell(Q);
                    const GaussianShell& Rshell = basisset_->shell(R);
                    const GaussianShell& Sshell = basisset_->shell(S);

                    int p_start = Pshell.start();
                    int num_p = Pshell.nfunction();

                    int q_start = Qshell.start();
                    int num_q = Qshell.nfunction();

                    int r_start = Rshell.start();
                    int num_r = Rshell.nfunction();

                    int s_start = Sshell.start();
                    int num_s = Sshell.nfunction();

                    int PQoff = offsets[P * nshell + Q];
                    int RSoff = offsets[R * nshell + S];

                    if (ints_[thread]->compute_shell(P, Q, R, S) == 0) continue;
                    const double* pqrs = ints_[thread]->buffer();

                    for (int N = 0; N < D_.size(); N++) {
                        double** Jp = J_[N]->pointer();
                        double** Dp = D_[N]->pointer();
                        double* JTp = JT[thread][N].data();
                        const double* pqrs2 = pqrs;
                        
                        if (!touched) {
                            ::memset((void*)(&JTp[0L * max_alloc]), '\0', max_alloc * sizeof(double));
                            ::memset((void*)(&JTp[1L * max_alloc]), '\0', max_alloc * sizeof(double));
                        }

                        double* J1p = &JTp[0L * max_alloc];
                        double* J2p = &JTp[1L * max_alloc];

                        for (int p = p_start; p < p_start + num_p; p++) {
                            int dp = p - p_start;
                            for (int q = q_start; q < q_start + num_q; q++) {
                                int dq = q - q_start;
                                for (int r = r_start; r < r_start + num_r; r++) {
                                    int dr = r - r_start;
                                    for (int s = s_start; s < s_start + num_s; s++) {
                                        int ds = s - s_start;

                                        int pq = PQoff + dp * num_q + dq;
                                        int rs = RSoff + dr * num_s + ds;
                                        
                                        J1p[pq] += prefactor * (*pqrs2) * Dp[r][s];
                                        J2p[rs] += prefactor * (*pqrs2) * Dp[p][q];
                                        pqrs2++;
                                    } // end s
                                } // end r
                            } // end q
                        } // end p

                    } // end N
                    touched = true;
                } // end RSind
            } // end PQind
            if (!touched) continue;
            
            // = > Stripeout < = //
            for (int N = 0; N < D_.size(); N++) {
                double** Jp = J_[N]->pointer();
                double** Dp = D_[N]->pointer();
                double* JTp = JT[thread][N].data();
                
                double* J1p = &JTp[0L * max_alloc];
                double* J2p = &JTp[1L * max_alloc];
                
                for (int PQind = 0; PQind < PQshells.size(); PQind++) {
                    std::pair<int, int> PQ = PQshells[PQind]->get_shell_pair_index();
                    int P = PQ.first;
                    int Q = PQ.second;

                    int PQoff = offsets[P * nshell + Q];
                
                    const GaussianShell& Pshell = basisset_->shell(P);
                    const GaussianShell& Qshell = basisset_->shell(Q);
                
                    int p_start = Pshell.start();
                    int num_p = Pshell.nfunction();

                    int q_start = Qshell.start();
                    int num_q = Qshell.nfunction();
                
                    for (int p = p_start; p < p_start + num_p; p++) {
                        int dp = p - p_start;
                        for (int q = q_start; q < q_start + num_q; q++) {
                            int dq = q - q_start;
                            int pq = PQoff + dp * num_q + dq;
#pragma omp atomic
                            Jp[p][q] += J1p[pq];
                        }
                    }
                }
            
                for (int RSind = 0; RSind < RSshells.size(); RSind++) {
                    std::pair<int, int> RS = RSshells[RSind]->get_shell_pair_index();
                    int R = RS.first;
                    int S = RS.second;

                    int RSoff = offsets[R * nshell + S];
                
                    const GaussianShell& Rshell = basisset_->shell(R);
                    const GaussianShell& Sshell = basisset_->shell(S);
                
                    int r_start = Rshell.start();
                    int num_r = Rshell.nfunction();

                    int s_start = Sshell.start();
                    int num_s = Sshell.nfunction();
                
                    for (int r = r_start; r < r_start + num_r; r++) {
                        int dr = r - r_start;
                        for (int s = s_start; s < s_start + num_s; s++) {
                            int ds = s - s_start;
                            int rs = RSoff + dr * num_s + ds;
#pragma omp atomic
                            Jp[r][s] += J2p[rs];
                        }
                    }
                }
            }
            
            // => End Stripeout <= //
            
        } // end nfi
    } // end bi

    timer_off("CFMMTree::build_nf_J");
}

void CFMMTree::build_ff_J() {

    timer_on("CFMMTree::build_ff_J");

    int start = (nlevels_ == 1) ? 0 : (0.5 * std::pow(16, nlevels_-1) + 7) / 15;
    int end = (nlevels_ == 1) ? 1 : (0.5 * std::pow(16, nlevels_) + 7) / 15;

#pragma omp parallel for
    for (int bi = start; bi < end; bi += 1) {
        std::shared_ptr<CFMMBox> curr = tree_[bi];
        std::shared_ptr<RealSolidHarmonics>& Vff = curr->far_field_vector();
        auto& PQshells = curr->get_shell_pairs();

        for (int PQind = 0; PQind < PQshells.size(); PQind++) {
            std::shared_ptr<ShellPair> PQshell = PQshells[PQind];
            std::pair<int, int> PQ = PQshell->get_shell_pair_index();
            std::vector<std::shared_ptr<RealSolidHarmonics>>& PQ_mpoles = PQshell->get_mpoles();

            int P = PQ.first;
            int Q = PQ.second;

            double prefactor = (P == Q) ? 0.5 : 1.0;

            const GaussianShell& Pshell = basisset_->shell(P);
            const GaussianShell& Qshell = basisset_->shell(Q);

            int p_start = Pshell.start();
            int num_p = Pshell.nfunction();

            int q_start = Qshell.start();
            int num_q = Qshell.nfunction();

            for (int p = p_start; p < p_start + num_p; p++) {
                int dp = p - p_start;
                for (int q = q_start; q < q_start + num_q; q++) {
                    int dq = q - q_start;
                    for (int N = 0; N < J_.size(); N++) {
                        double** Jp = J_[N]->pointer();
                        // Far field multipole contributions
                        for (int l = 0; l <= lmax_; l++) {
                            for (int mu = 0; mu < 2*l+1; mu++) {
                                Jp[p][q] += prefactor * Vff->get(l, mu) * PQ_mpoles[dp * num_q + dq]->get(l, mu);
                            } // end mu
                        } // end l
                    } // end N
                } // end q
            } // end p

        } // end PQind
    } // end bi

    timer_off("CFMMTree::build_ff_J");

}

void CFMMTree::build_J() {
    timer_on("CFMMTree::build_J");

    // Zero the J matrix
    for (int ind = 0; ind < D_.size(); ind++) {
        J_[ind]->zero();
    }

    // Compute multipoles and far field
    calculate_multipoles();
    compute_far_field();

    // Compute near field J and far field J
    build_nf_J();
    build_ff_J();

    // Hermitivitize J matrix afterwards
    for (int ind = 0; ind < D_.size(); ind++) {
        J_[ind]->hermitivitize();
    }

    timer_off("CFMMTree::build_J");
}

void CFMMTree::print_out() {
    for (int bi = 0; bi < tree_.size(); bi++) {
        std::shared_ptr<CFMMBox> box = tree_[bi];
        auto sp = box->get_shell_pairs();
        int nshells = sp.size();
        int level = box->get_level();
        int ws = box->get_ws();
        if (nshells > 0) {
            outfile->Printf("  BOX INDEX: %d, LEVEL: %d, WS: %d, NSHELLS: %d\n", bi, level, ws, nshells);
            for (int si = 0; si < sp.size(); si++) {
                Vector3 center = sp[si]->get_center();
                outfile->Printf("  SHELL: %d, x: %8.5f, y: %8.5f, z: %8.5f\n\n", si, center[0], center[1], center[2]);
            }
        }
    }
}

} // end namespace psi