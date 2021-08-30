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

std::shared_ptr<CFMMBox> CFMMBox::get() {
    return shared_from_this();
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

void CFMMBox::set_nf_lff() {

    // Creates a temporary parent shared pointer
    std::shared_ptr<CFMMBox> parent = parent_.lock();

    // Parent is not a nullpointer
    if (parent) {
        // Siblings of this box (Technically near fields include self in this implementation)
        for (std::shared_ptr<CFMMBox> sibling : parent->children_) {
            Vector3 Rab = center_ - sibling->center_;
            double dist = Rab.norm();

            int ref_ws = (ws_ + sibling->ws_) / 2;
            if (dist <= ref_ws * length_ * std::sqrt(3.0)) {
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
                double dist = Rab.norm();

                int ref_ws = (ws_ + cousin->ws_) / 2;
                if (dist <= ref_ws * length_ * std::sqrt(3.0)) {
                    near_field_.push_back(cousin);
                } else {
                    local_far_field_.push_back(cousin);
                }
            }
        }
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

    // Compute multipoles for all basis sets in the basis pair
#pragma omp parallel for
    for (int ind = 0; ind < shell_pairs_.size(); ind++) {
        std::shared_ptr<ShellPair> sp = shell_pairs_[ind];
        int thread = 0;

#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif
        sp->calculate_mpoles(center_, sints[thread], mpints[thread], lmax_);
        std::vector<std::shared_ptr<RealSolidHarmonics>>& sp_mpoles = sp->get_mpoles();

        std::pair<int, int> PQ = sp->get_shell_pair_index();
        int P = PQ.first;
        int Q = PQ.second;

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
                    basis_mpole->scale(2.0 * D[N]->get(p, q));
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

void CFMMBox::compute_far_field_vector() {

    for (std::shared_ptr<CFMMBox> box : local_far_field_) {
        // The far field effect the boxes have on this particular box
        std::shared_ptr<RealSolidHarmonics> far_field = box->mpoles_->far_field_vector(center_);
        Vff_->add(far_field);
    }

    // If parent is not null, add the parent's far field
    // Creates a temporary parent shared pointer
    std::shared_ptr<CFMMBox> parent = parent_.lock();

    // If parent is not null, add the parent's far field
    if (parent) {
        Vff_->add(parent->Vff_->translate(center_));
    }

}

void CFMMBox::compute_nf_J(std::shared_ptr<BasisSet> basisset, std::vector<SharedMatrix>& D, std::vector<SharedMatrix>& J) {

    std::vector<std::shared_ptr<TwoBodyAOInt>> ints;
    std::shared_ptr<IntegralFactory> factory = std::make_shared<IntegralFactory>(basisset);
    std::shared_ptr<TwoBodyAOInt> eri = std::shared_ptr<TwoBodyAOInt>(std::shared_ptr<TwoBodyAOInt>(factory->eri()));
    ints.push_back(eri);

    for (int thread = 1; thread < nthread_; thread++) {
        ints.push_back(std::shared_ptr<TwoBodyAOInt>(eri->clone()));
    }

#pragma omp parallel for
    for (int indA = 0; indA < shell_pairs_.size(); indA++) {
        std::shared_ptr<ShellPair> spA = shell_pairs_[indA];
        std::pair<int, int> PQ = spA->get_shell_pair_index();
        int P = PQ.first;
        int Q = PQ.second;

        int thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif

        const GaussianShell& Pshell = basisset->shell(P);
        const GaussianShell& Qshell = basisset->shell(Q);

        int p_start = Pshell.start();
        int num_p = Pshell.nfunction();

        int q_start = Qshell.start();
        int num_q = Qshell.nfunction();

        for (int N = 0; N < D.size(); N++) {
            double** Jp = J[N]->pointer();
            double** Dp = D[N]->pointer();

            for (int nf = 0; nf < near_field_.size(); nf++) {
                std::shared_ptr<CFMMBox> nfbox = near_field_[nf];
                for (int indB = 0; indB < nfbox->shell_pairs_.size(); indB++) {
                    std::shared_ptr<ShellPair> spB = nfbox->shell_pairs_[indB];
                    std::pair<int, int> RS = spB->get_shell_pair_index();
                    int R = RS.first;
                    int S = RS.second;

                    const GaussianShell& Rshell = basisset->shell(R);
                    const GaussianShell& Sshell = basisset->shell(S);

                    int r_start = Rshell.start();
                    int num_r = Rshell.nfunction();

                    int s_start = Sshell.start();
                    int num_s = Sshell.nfunction();

                    ints[thread]->compute_shell(P, Q, R, S);
                    const double* pqrs = ints[thread]->buffer();

                    for (int p = p_start; p < p_start + num_p; p++) {
                        int dp = p - p_start;
                        for (int q = q_start; q < q_start + num_q; q++) {
                            int dq = q - q_start;
                            for (int r = r_start; r < r_start + num_r; r++) {
                                int dr = r - r_start;
                                for (int s = s_start; s < s_start + num_s; s++) {
                                    int ds = s - s_start;
                                    Jp[p][q] += pqrs[dp * num_q * num_r * num_s + dq * num_r * num_s + dr * num_s + ds] * Dp[r][s];
                                }
                            }
                        }
                    }
                }
            }
        }
    }

}

void CFMMBox::compute_ff_J(std::shared_ptr<BasisSet> basisset, std::vector<SharedMatrix>& D, std::vector<SharedMatrix>& J) {

#pragma omp parallel for
    for (int ind = 0; ind < shell_pairs_.size(); ind++) {
        std::shared_ptr<ShellPair> sp = shell_pairs_[ind];
        std::vector<std::shared_ptr<RealSolidHarmonics>>& sp_mpoles = sp->get_mpoles();

        int thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif

        std::pair<int, int> PQ = sp->get_shell_pair_index();
        int P = PQ.first;
        int Q = PQ.second;

        const GaussianShell& Pshell = basisset->shell(P);
        const GaussianShell& Qshell = basisset->shell(Q);

        int p_start = Pshell.start();
        int num_p = Pshell.nfunction();

        int q_start = Qshell.start();
        int num_q = Qshell.nfunction();

        for (int N = 0; N < D.size(); N++) {
            double** Jp = J[N]->pointer();
            double** Dp = D[N]->pointer();

            for (int p = p_start; p < p_start + num_p; p++) {
                int dp = p - p_start;
                for (int q = q_start; q < q_start + num_q; q++) {
                    int dq = q - q_start;
                    double val = 0.0;
                    for (int l = 0; l <= lmax_; l++) {
                        for (int mu = 0; mu < 2*l+1; mu++) {
                            val += Vff_->get(l, mu) * sp_mpoles[dp * num_q + dq]->get(l, mu);
                        }
                    }
                    Jp[p][q] += 0.5 * val;
                }
            }
        }
    }
}

void CFMMBox::compute_J(std::shared_ptr<BasisSet> basisset, std::vector<SharedMatrix>& D, std::vector<SharedMatrix>& J) {
    compute_nf_J(basisset, D, J);
    compute_ff_J(basisset, D, J);
}

CFMMTree::CFMMTree(std::shared_ptr<Molecule> molecule, std::shared_ptr<BasisSet> basisset, std::vector<SharedMatrix>& D, 
                std::vector<SharedMatrix>& J, const std::vector<std::pair<int, int>>& shell_pairs, int nlevels, int lmax) {
    molecule_ = molecule;
    basisset_ = basisset;
    D_ = D;
    J_ = J;
    nlevels_ = nlevels;
    lmax_ = lmax;
    // int num_boxes = (0.5 * std::pow(16, nlevels_) + 7) / 15;
    Options& options = Process::environment.options;
    mpole_coefs_ = std::make_shared<HarmonicCoefficients>(options.get_int("CFMM_MAX_MPOLE_ORDER"), Regular);

    for (const auto& pair : shell_pairs) {
        shell_pairs_.push_back(std::make_shared<ShellPair>(basisset_, pair, mpole_coefs_));
        if (pair.first != pair.second) {
            std::pair<int, int> reverse_pair = std::make_pair(pair.second, pair.first);
            shell_pairs_.push_back(std::make_shared<ShellPair>(basisset_, reverse_pair, mpole_coefs_));
        }
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

    tree_.push_back(std::make_shared<CFMMBox>(nullptr, shell_pairs_, origin, length, 0, lmax_, 2));
}

void CFMMTree::make_children() {

    int bi = 0;
    int end = (0.5 * std::pow(16, nlevels_ - 1) + 7) / 15;

    while (bi < end) {
        tree_[bi]->make_children();
        for (std::shared_ptr<CFMMBox> child : tree_[bi]->get_children()) {
            tree_.push_back(child);
        }
        bi += 1;
    }
}

void CFMMTree::calculate_multipoles() {
    timer_on("CFMMTree::calculate_multipoles");

    for (int bi = tree_.size() - 1; bi >= 0; bi -= 1) {
        std::shared_ptr<CFMMBox> box = tree_[bi];
        int level = box->get_level();
        if (level == nlevels_ - 1) box->compute_mpoles(basisset_, D_);
        else box->compute_mpoles_from_children();
    }

    timer_off("CFMMTree::calculate_multipoles");
}

void CFMMTree::set_nf_lff() {
    for (int bi = 0; bi < tree_.size(); bi++) {
        tree_[bi]->set_nf_lff();
    }
}

void CFMMTree::compute_far_field() {

    timer_on("CFMMTree::compute_far_field");
    for (int bi = 0; bi < tree_.size(); bi++) {
        tree_[bi]->compute_far_field_vector();
    }
    timer_off("CFMMTree::compute_far_field");
}

void CFMMTree::build_J() {
    // Zero the J matrix

    timer_on("CFMMTree::build_J");

    for (int ind = 0; ind < D_.size(); ind++) {
        J_[ind]->zero();
    }

    calculate_multipoles();
    set_nf_lff();
    compute_far_field();

    for (int bi = tree_.size() - 1; bi >= 0; bi -= 1) {
        std::shared_ptr<CFMMBox> box = tree_[bi];
        int level = box->get_level();
        if (level == nlevels_ - 1) box->compute_J(basisset_, D_, J_);
        else break;
    }

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