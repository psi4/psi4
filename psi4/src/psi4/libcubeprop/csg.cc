/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2024 The Psi4 Developers.
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

#include "csg.h"

#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "psi4/psi4-dec.h"

#include "psi4/libciomr/libciomr.h"
#include "psi4/libfilesystem/path.h"
#include "psi4/libfock/cubature.h"
#include "psi4/libfock/points.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/potential.h"
#include "psi4/libmints/vector.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/libqt/qt.h"

namespace psi {

CubicScalarGrid::CubicScalarGrid(std::shared_ptr<BasisSet> primary, Options& options)
    : primary_(primary), mol_(primary->molecule()), options_(options) {
    filepath_ = "";
    npoints_ = 0L;
    x_ = nullptr;
    y_ = nullptr;
    z_ = nullptr;
    w_ = nullptr;

    build_grid();  // Defaults from Options
}
CubicScalarGrid::~CubicScalarGrid() {
    if (x_) delete[] x_;
    if (y_) delete[] y_;
    if (z_) delete[] z_;
    if (w_) delete[] w_;
}
void CubicScalarGrid::build_grid(const std::string filepath, int* N, double* D, double* O) {
    filepath_ = filepath;

    for (int k = 0; k < 3; k++) {
        N_[k] = N[k];
        O_[k] = O[k];
        D_[k] = D[k];
    }

    populate_grid();
}
void CubicScalarGrid::build_grid(std::shared_ptr<CubicScalarGrid> other) {
    filepath_ = other->filepath();

    for (int k = 0; k < 3; k++) {
        N_[k] = other->N()[k];
        O_[k] = other->O()[k];
        D_[k] = other->D()[k];
    }

    populate_grid();
}
void CubicScalarGrid::build_grid() {
    filepath_ = ".";

    double L[3];
    if (options_["CUBIC_GRID_OVERAGE"].size() != 3) {
        L[0] = 4.0;
        L[1] = 4.0;
        L[2] = 4.0;
    } else {
        L[0] = options_["CUBIC_GRID_OVERAGE"][0].to_double();
        L[1] = options_["CUBIC_GRID_OVERAGE"][1].to_double();
        L[2] = options_["CUBIC_GRID_OVERAGE"][2].to_double();
    }

    double D[3];
    if (options_["CUBIC_GRID_SPACING"].size() != 3) {
        D[0] = 0.2;
        D[1] = 0.2;
        D[2] = 0.2;
    } else {
        D[0] = options_["CUBIC_GRID_SPACING"][0].to_double();
        D[1] = options_["CUBIC_GRID_SPACING"][1].to_double();
        D[2] = options_["CUBIC_GRID_SPACING"][2].to_double();
    }

    double Xmin[3];
    double Xmax[3];
    double Xdel[3];
    int N[3];
    double O[3];

    for (int k = 0; k < 3; k++) {
        Xmin[k] = Xmax[k] = mol_->xyz(0, k);
        for (int A = 0; A < mol_->natom(); A++) {
            Xmin[k] = (Xmin[k] > mol_->xyz(A, k) ? mol_->xyz(A, k) : Xmin[k]);
            Xmax[k] = (Xmax[k] < mol_->xyz(A, k) ? mol_->xyz(A, k) : Xmax[k]);
        }
        Xdel[k] = Xmax[k] - Xmin[k];
        N[k] = (int)((Xmax[k] - Xmin[k] + 2.0 * L[k]) / (D[k]));
        if (D[k] * (double)N[k] < (Xmax[k] - Xmin[k] + 2.0 * L[k])) N[k]++;
        O[k] = Xmin[k] - (D[k] * (double)N[k] - (Xmax[k] - Xmin[k])) / 2.0;
        N_[k] = N[k];
        O_[k] = O[k];
        D_[k] = D[k];
    }

    populate_grid();
}
void CubicScalarGrid::populate_grid() {
    if (x_) delete[] x_;
    if (y_) delete[] y_;
    if (z_) delete[] z_;
    if (w_) delete[] w_;

    npoints_ = (N_[0] + 1L) * (N_[1] + 1L) * (N_[2] + 1L);
    x_ = new double[npoints_];
    y_ = new double[npoints_];
    z_ = new double[npoints_];
    w_ = new double[npoints_];

    double epsilon = options_.get_double("CUBIC_BASIS_TOLERANCE");
    extents_ = std::make_shared<BasisExtents>(primary_, epsilon);

    int max_points = options_.get_int("CUBIC_BLOCK_MAX_POINTS");
    nxyz_ = std::llround(pow((double)max_points, 1.0 / 3.0));

    blocks_.clear();
    size_t offset = 0L;
    for (int istart = 0L; istart <= N_[0]; istart += nxyz_) {
        int ni = (istart + nxyz_ > N_[0] ? (N_[0] + 1) - istart : nxyz_);
        for (int jstart = 0L; jstart <= N_[1]; jstart += nxyz_) {
            int nj = (jstart + nxyz_ > N_[1] ? (N_[1] + 1) - jstart : nxyz_);
            for (int kstart = 0L; kstart <= N_[2]; kstart += nxyz_) {
                int nk = (kstart + nxyz_ > N_[2] ? (N_[2] + 1) - kstart : nxyz_);

                double* xp = &x_[offset];
                double* yp = &y_[offset];
                double* zp = &z_[offset];
                double* wp = &w_[offset];

                size_t block_size = 0L;
                for (int i = istart; i < istart + ni; i++) {
                    for (int j = jstart; j < jstart + nj; j++) {
                        for (int k = kstart; k < kstart + nk; k++) {
                            x_[offset] = O_[0] + i * D_[0];
                            y_[offset] = O_[1] + j * D_[1];
                            z_[offset] = O_[2] + k * D_[2];
                            w_[offset] = D_[0] * D_[1] * D_[2];
                            offset++;
                            block_size++;
                        }
                    }
                }
                blocks_.push_back(std::make_shared<BlockOPoints>(0, block_size, xp, yp, zp, wp, extents_));
            }
        }
    }

    int max_functions = 0L;
    for (int ind = 0; ind < blocks_.size(); ind++) {
        max_functions = (max_functions >= blocks_[ind]->functions_local_to_global().size()
                             ? max_functions
                             : blocks_[ind]->functions_local_to_global().size());
    }

    points_ = std::make_shared<RKSFunctions>(primary_, max_points, max_functions);
    points_->set_ansatz(0);
}
void CubicScalarGrid::print_header() {
    outfile->Printf("  ==> CubicScalarGrid <==\n\n");

    outfile->Printf("    Filepath     = %s\n", filepath_.c_str());
    outfile->Printf("    Total Points = %16zu\n", npoints_);
    outfile->Printf("    XYZ Blocking = %16zu\n", nxyz_);
    outfile->Printf("    X Points     = %16zu\n", N_[0] + 1L);
    outfile->Printf("    Y Points     = %16zu\n", N_[1] + 1L);
    outfile->Printf("    Z Points     = %16zu\n", N_[2] + 1L);
    outfile->Printf("    X Spacing    = %16.3E\n", D_[0]);
    outfile->Printf("    Y Spacing    = %16.3E\n", D_[1]);
    outfile->Printf("    Z Spacing    = %16.3E\n", D_[2]);
    outfile->Printf("    X Minimum    = %16.3E\n", O_[0]);
    outfile->Printf("    Y Minimum    = %16.3E\n", O_[1]);
    outfile->Printf("    Z Minimum    = %16.3E\n", O_[2]);
    outfile->Printf("    X Maximum    = %16.3E\n", O_[0] + D_[0] * N_[0]);
    outfile->Printf("    Y Maximum    = %16.3E\n", O_[1] + D_[1] * N_[1]);
    outfile->Printf("    Z Maximum    = %16.3E\n", O_[2] + D_[2] * N_[2]);
    outfile->Printf("\n");

    primary_->print();
}
void CubicScalarGrid::write_gen_file(double* v, const std::string& name, const std::string& type,
                                     const std::string& comment) {
    if (type == "CUBE") {
        write_cube_file(v, name, comment);
    } else {
        throw PSIEXCEPTION("CubicScalarGrid: Unrecognized output file type");
    }
}
void CubicScalarGrid::write_cube_file(double* v, const std::string& name, const std::string& comment) {
    // => Reorder the grid <= //

    auto v2 = std::vector<double>(npoints_);
    size_t offset = 0L;
    for (int istart = 0L; istart <= N_[0]; istart += nxyz_) {
        int ni = (istart + nxyz_ > N_[0] ? (N_[0] + 1) - istart : nxyz_);
        for (int jstart = 0L; jstart <= N_[1]; jstart += nxyz_) {
            int nj = (jstart + nxyz_ > N_[1] ? (N_[1] + 1) - jstart : nxyz_);
            for (int kstart = 0L; kstart <= N_[2]; kstart += nxyz_) {
                int nk = (kstart + nxyz_ > N_[2] ? (N_[2] + 1) - kstart : nxyz_);
                for (int i = istart; i < istart + ni; i++) {
                    for (int j = jstart; j < jstart + nj; j++) {
                        for (int k = kstart; k < kstart + nk; k++) {
                            size_t index = i * (N_[1] + 1L) * (N_[2] + 1L) + j * (N_[2] + 1L) + k;
                            v2[index] = v[offset];
                            offset++;
                        }
                    }
                }
            }
        }
    }

    // => Drop the grid out <= //

    std::stringstream ss;
    ss << filepath_ << "/" << name << ".cube";

    // Is filepath a valid directory?
    if (filesystem::path(filepath_).make_absolute().is_directory() == false) {
        printf("Filepath \"%s\" is not valid.  Please create this directory.\n", filepath_.c_str());
        outfile->Printf("Filepath \"%s\" is not valid.  Please create this directory.\n", filepath_.c_str());
        exit(Failure);
    }

    FILE* fh = fopen(ss.str().c_str(), "w");
    // Two comment lines
    fprintf(fh, "Psi4 Gaussian Cube File.\n");
    fprintf(fh, "Property: %s%s\n", name.c_str(), comment.c_str());

    // Number of atoms plus origin of data
    fprintf(fh, "%6d %10.6f %10.6f %10.6f\n", mol_->natom(), O_[0], O_[1], O_[2]);

    // Number of points along axis, displacement along x,y,z
    fprintf(fh, "%6d %10.6f %10.6f %10.6f\n", N_[0] + 1, D_[0], 0.0, 0.0);
    fprintf(fh, "%6d %10.6f %10.6f %10.6f\n", N_[1] + 1, 0.0, D_[1], 0.0);
    fprintf(fh, "%6d %10.6f %10.6f %10.6f\n", N_[2] + 1, 0.0, 0.0, D_[2]);

    // Atoms of molecule (Z, Q?, x, y, z)
    for (int A = 0; A < mol_->natom(); A++) {
        fprintf(fh, "%3d %10.6f %10.6f %10.6f %10.6f\n", (int)mol_->true_atomic_number(A), 0.0, mol_->x(A), mol_->y(A),
                mol_->z(A));
    }

    // Data, striped (x, y, z)
    for (size_t ind = 0; ind < npoints_; ind++) {
        fprintf(fh, "%12.5E ", v2[ind]);
        if (ind % 6 == 5) fprintf(fh, "\n");
    }

    fclose(fh);
}
void CubicScalarGrid::add_density(double* v, std::shared_ptr<Matrix> D) {
    points_->set_pointers(D);
    std::shared_ptr<Vector> rho = points_->point_value("RHO_A");
    double* rhop = rho->pointer();

    size_t offset = 0L;
    for (int ind = 0; ind < blocks_.size(); ind++) {
        points_->compute_points(blocks_[ind]);
        size_t npoints = blocks_[ind]->npoints();
        C_DAXPY(npoints, 0.5, rhop, 1, &v[offset], 1);
        offset += npoints;
    }
}
void CubicScalarGrid::add_esp(double* v, std::shared_ptr<Matrix> D, const std::vector<double>& nuc_weights) {
    // => Auxiliary Basis Set <= //

    if (!auxiliary_) {
        throw PSIEXCEPTION("Auxiliary basis is required for ESP computations.");
    }

    double cutoff = options_.get_double("INTS_TOLERANCE");
    double condition = options_.get_double("DF_FITTING_CONDITION");

    // => Sizing <= //

    int nbf = primary_->nbf();
    int naux = auxiliary_->nbf();
    int maxP = auxiliary_->max_function_per_shell();

    int nthreads = 1;
#ifdef _OPENMP
    nthreads = Process::environment.get_n_threads();
#endif

    // => Density Fitting (TODO: Could be sped up) <= //

    std::shared_ptr<IntegralFactory> Ifact =
        std::make_shared<IntegralFactory>(auxiliary_, BasisSet::zero_ao_basis_set(), primary_, primary_);
    std::vector<std::shared_ptr<TwoBodyAOInt>> ints;
    for (int thread = 0; thread < nthreads; thread++) {
        ints.push_back(std::shared_ptr<TwoBodyAOInt>(Ifact->eri()));
    }

    const std::vector<std::pair<int, int>>& pairs = ints[0]->shell_pairs();

    auto c = std::make_shared<Vector>("c", naux);
    double* cp = c->pointer();

    auto Amn = std::make_shared<Matrix>("Amn", maxP, nbf * nbf);
    double** Amnp = Amn->pointer();

    double** Dp = D->pointer();

    for (int P = 0; P < auxiliary_->nshell(); P++) {
        int nP = auxiliary_->shell(P).nfunction();
        int oP = auxiliary_->shell(P).function_index();

        Amn->zero();

// Integrals
#pragma omp parallel for schedule(dynamic)
        for (int task = 0; task < pairs.size(); task++) {
            int thread = 0;
#ifdef _OPENMP
            thread = omp_get_thread_num();
#endif

            int M = pairs[task].first;
            int N = pairs[task].second;

            ints[thread]->compute_shell(P, 0, M, N);
            const double* buffer = ints[thread]->buffer();

            int nM = primary_->shell(M).nfunction();
            int oM = primary_->shell(M).function_index();
            int nN = primary_->shell(N).nfunction();
            int oN = primary_->shell(N).function_index();

            int index = 0;
            for (int p = 0; p < nP; p++) {
                for (int m = 0; m < nM; m++) {
                    for (int n = 0; n < nN; n++) {
                        Amnp[p][(m + oM) * nbf + (n + oN)] = Amnp[p][(n + oN) * nbf + (m + oM)] = buffer[index++];
                    }
                }
            }
        }

        // Contraction
        C_DGEMV('N', nP, nbf * nbf, 1.0, Amnp[0], nbf * nbf, Dp[0], 1, 0.0, cp + oP, 1);
    }

    Amn.reset();
    Ifact.reset();
    ints.clear();

    auto J = std::make_shared<Matrix>("J", naux, naux);
    double** Jp = J->pointer();

    std::shared_ptr<IntegralFactory> Jfact = std::make_shared<IntegralFactory>(
        auxiliary_, BasisSet::zero_ao_basis_set(), auxiliary_, BasisSet::zero_ao_basis_set());

    std::shared_ptr<TwoBodyAOInt> Jints(Jfact->eri());

    for (int P = 0; P < auxiliary_->nshell(); P++) {
        int nP = auxiliary_->shell(P).nfunction();
        int oP = auxiliary_->shell(P).function_index();

        for (int Q = 0; Q <= P; Q++) {
            int nQ = auxiliary_->shell(Q).nfunction();
            int oQ = auxiliary_->shell(Q).function_index();

            Jints->compute_shell(P, 0, Q, 0);
            const double* Jbuffer = Jints->buffer();

            int index = 0;
            for (int p = 0; p < nP; p++) {
                for (int q = 0; q < nQ; q++) {
                    Jp[p + oP][q + oQ] = Jp[q + oQ][p + oP] = Jbuffer[index++];
                }
            }
        }
    }

    Jfact.reset();
    Jints.reset();

    J->power(-1.0, condition);

    auto d = std::make_shared<Vector>("d", naux);
    double* dp = d->pointer();

    C_DGEMV('N', naux, naux, 1.0, Jp[0], naux, cp, 1, 0.0, dp, 1);

    // c->print();
    // d->print();

    J.reset();

    // => Electronic Part <= //

    std::shared_ptr<IntegralFactory> Vfact = std::make_shared<IntegralFactory>(
        auxiliary_, BasisSet::zero_ao_basis_set(), auxiliary_, BasisSet::zero_ao_basis_set());
    std::vector<std::shared_ptr<Matrix>> VtempT;
    std::vector<std::shared_ptr<PotentialInt>> VintT;
    for (int thread = 0; thread < nthreads; thread++) {
        VtempT.push_back(std::make_shared<Matrix>("Vtemp", naux, 1));
        VintT.push_back(std::shared_ptr<PotentialInt>(static_cast<PotentialInt*>(Vfact->ao_potential().release())));
    }

#pragma omp parallel for schedule(dynamic)
    for (int P = 0; P < npoints_; P++) {
        // Thread info
        int thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif
        double** VtempTp = VtempT[thread]->pointer();

        // Integrals
        VtempT[thread]->zero();
        VintT[thread]->set_charge_field({{1.0, {x_[P], y_[P], z_[P]}}});
        VintT[thread]->compute(VtempT[thread]);

        // Contraction
        v[P] += C_DDOT(naux, dp, 1, VtempTp[0], 1);  // Potential integrals are negative definite already
    }

    // => Nuclear Part <= //

    for (int A = 0; A < mol_->natom(); A++) {
        double Z = mol_->Z(A) * (nuc_weights.size() ? nuc_weights[A] : 1.0);
        double x = mol_->x(A);
        double y = mol_->y(A);
        double z = mol_->z(A);
        for (int P = 0; P < npoints_; P++) {
            double R = sqrt((x - x_[P]) * (x - x_[P]) + (y - y_[P]) * (y - y_[P]) + (z - z_[P]) * (z - z_[P]));
            v[P] += (R >= 1.0E-15 ? Z / R : 0.0);
        }
    }
}
void CubicScalarGrid::add_basis_functions(double** v, const std::vector<int>& indices) {
    std::shared_ptr<Matrix> phi = points_->basis_value("PHI");
    double** phip = phi->pointer();

    size_t offset = 0L;
    for (int ind = 0; ind < blocks_.size(); ind++) {
        points_->compute_functions(blocks_[ind]);

        size_t npoints = blocks_[ind]->npoints();
        const std::vector<int>& function_map = blocks_[ind]->functions_local_to_global();
        int nlocal = function_map.size();
        int nglobal = points_->max_functions();

        for (int ind1 = 0; ind1 < indices.size(); ind1++) {
            for (int ind2 = 0; ind2 < function_map.size(); ind2++) {
                if (indices[ind1] == function_map[ind2]) {
                    C_DAXPY(npoints, 1.0, &phip[0][ind2], nglobal, &v[ind1][offset], 1);
                }
            }
        }

        offset += npoints;
    }
}
void CubicScalarGrid::add_orbitals(double** v, std::shared_ptr<Matrix> C) {
    int na = C->colspi()[0];

    points_->set_Cs(C);
    std::shared_ptr<Matrix> psi = points_->orbital_value("PSI_A");
    double** psip = psi->pointer();

    size_t offset = 0L;
    for (int ind = 0; ind < blocks_.size(); ind++) {
        points_->compute_orbitals(blocks_[ind]);

        size_t npoints = blocks_[ind]->npoints();
        for (int a = 0; a < na; a++) {
            C_DAXPY(npoints, 1.0, psip[a], 1, &v[a][offset], 1);
        }

        offset += npoints;
    }
}
void CubicScalarGrid::add_LOL(double* v, std::shared_ptr<Matrix> D) {
    points_->set_ansatz(2);
    points_->set_pointers(D);

    std::shared_ptr<Vector> rho = points_->point_value("RHO_A");
    std::shared_ptr<Vector> tau = points_->point_value("TAU_A");
    double* rhop = rho->pointer();
    double* taup = tau->pointer();

    double C = 3.0 / 5.0 * pow(6.0 * M_PI * M_PI, 2.0 / 3.0);

    size_t offset = 0L;
    for (int ind = 0; ind < blocks_.size(); ind++) {
        points_->compute_points(blocks_[ind]);
        size_t npoints = blocks_[ind]->npoints();
        for (int P = 0; P < npoints; P++) {
            double tau_LSDA = C * pow(0.5 * rhop[P], 5.0 / 3.0);
            double tau_EX = taup[P];
            double t = tau_LSDA / tau_EX;
            double v2 = (std::fabs(tau_EX / tau_LSDA) < 1.0E-15 ? 1.0 : t / (1.0 + t));
            v[P + offset] += v2;
        }
        offset += npoints;
    }

    points_->set_ansatz(0);
}
void CubicScalarGrid::add_ELF(double* v, std::shared_ptr<Matrix> D) {
    points_->set_ansatz(2);
    points_->set_pointers(D);

    std::shared_ptr<Vector> rho = points_->point_value("RHO_A");
    std::shared_ptr<Vector> gam = points_->point_value("GAMMA_AA");
    std::shared_ptr<Vector> tau = points_->point_value("TAU_A");
    double* rhop = rho->pointer();
    double* gamp = gam->pointer();
    double* taup = tau->pointer();

    double C = 3.0 / 5.0 * pow(6.0 * M_PI * M_PI, 2.0 / 3.0);

    size_t offset = 0L;
    for (int ind = 0; ind < blocks_.size(); ind++) {
        points_->compute_points(blocks_[ind]);
        size_t npoints = blocks_[ind]->npoints();
        for (int P = 0; P < npoints; P++) {
            double tau_LSDA = C * pow(0.5 * rhop[P], 5.0 / 3.0);
            double tau_EX = taup[P];
            double D_EX = tau_EX - 0.25 * gamp[P] / rhop[P];
            double D_LSDA = tau_LSDA;
            double B = D_EX / D_LSDA;
            double v2 = (std::fabs(D_LSDA / D_EX) < 1.0E-15 ? 0.0 : 1.0 / (1.0 + B * B));
            v[P + offset] += v2;
        }
        offset += npoints;
    }

    points_->set_ansatz(0);
}
void CubicScalarGrid::compute_density(std::shared_ptr<Matrix> D, const std::string& name, const std::string& type) {
    auto v = std::vector<double>(npoints_, 0);
    add_density(v.data(), D);
    // Get adaptive isocountour range
    std::pair<double, double> isocontour_range = compute_isocontour_range(v.data(), 1.0);
    double density_percent = 100.0 * options_.get_double("CUBEPROP_ISOCONTOUR_THRESHOLD");
    std::stringstream comment;
    comment << " [e/a0^3]. Isocontour range for " << density_percent << "% of the density: (" << isocontour_range.first
            << "," << isocontour_range.second << ")." << ecp_header();
    // Write to disk
    write_gen_file(v.data(), name, type, comment.str());
}
void CubicScalarGrid::compute_esp(std::shared_ptr<Matrix> D, const std::vector<double>& w, const std::string& name,
                                  const std::string& type) {
    auto v = std::vector<double>(npoints_, 0);
    add_esp(v.data(), D, w);
    write_gen_file(v.data(), name, type, " [Eh/e]");
}
void CubicScalarGrid::compute_basis_functions(const std::vector<int>& indices, const std::string& name,
                                              const std::string& type) {
    auto v = Matrix(indices.size(), npoints_);
    auto vp = v.pointer();
    add_basis_functions(vp, indices);
    for (int k = 0; k < indices.size(); k++) {
        std::stringstream ss;
        ss << name << "_" << (indices[k] + 1);
        write_gen_file(vp[k], ss.str(), type);
    }
}
void CubicScalarGrid::compute_orbitals(std::shared_ptr<Matrix> C, const std::vector<int>& indices,
                                       const std::vector<std::string>& labels, const std::string& name,
                                       const std::string& type) {
    auto C2 = std::make_shared<Matrix>(primary_->nbf(), indices.size());
    double** Cp = C->pointer();
    double** C2p = C2->pointer();
    for (int k = 0; k < indices.size(); k++) {
        C_DCOPY(primary_->nbf(), &Cp[0][indices[k]], C->colspi()[0], &C2p[0][k], C2->colspi()[0]);
    }
    auto v = Matrix(indices.size(), npoints_);
    auto vp = v.pointer();
    add_orbitals(vp, C2);
    for (int k = 0; k < indices.size(); k++) {
        // Get adaptive isocountour range
        std::pair<double, double> isocontour_range = compute_isocontour_range(vp[k], 2.0);
        double density_percent = 100.0 * options_.get_double("CUBEPROP_ISOCONTOUR_THRESHOLD");
        std::stringstream comment;
        comment << ". Isocontour range for " << density_percent << "% of the density: (" << isocontour_range.first
                << "," << isocontour_range.second << ")";
        // Write to disk
        std::stringstream ss;
        ss << name << "_" << (indices[k] + 1) << "_" << labels[k];
        write_gen_file(vp[k], ss.str(), type, comment.str());
    }
}
void CubicScalarGrid::compute_difference(std::shared_ptr<Matrix> C, const std::vector<int>& indices,
                                         const std::string& label, const bool square, const std::string& type) {
    auto C2 = std::make_shared<Matrix>(primary_->nbf(), indices.size());
    double** Cp = C->pointer();
    double** C2p = C2->pointer();
    for (int k = 0; k < indices.size(); k++) {
        C_DCOPY(primary_->nbf(), &Cp[0][indices[k]], C->colspi()[0], &C2p[0][k], C2->colspi()[0]);
    }
    auto v_t = std::make_shared<Matrix>(indices.size(), npoints_);
    double** v_tp = v_t->pointer();
    auto v = std::make_shared<Vector>(npoints_);
    double* vp = v->pointer();
    add_orbitals(&v_tp[0], C2);
    for (int i = 0; i < npoints_; i++) {
         if (square) {
             v->set(0, i, (v_t->get(0,i) - v_t->get(1,i))*(v_t->get(0,i) + v_t->get(1,i)));
         } else {
             v->set(0, i, (v_t->get(0,i) - v_t->get(1,i)));
         }
    }
    std::pair<double, double> isocontour_range = compute_isocontour_range(&vp[0], 2.0);
    double density_percent = 100.0 * options_.get_double("CUBEPROP_ISOCONTOUR_THRESHOLD");
    std::stringstream comment;
    comment << ". Isocontour range for " << density_percent << "% of the density: (" << isocontour_range.first
            << "," << isocontour_range.second << ")";
    // Write to disk
    write_gen_file(&vp[0], label, type, comment.str());
}
void CubicScalarGrid::compute_LOL(std::shared_ptr<Matrix> D, const std::string& name, const std::string& type) {
    auto v = std::vector<double>(npoints_, 0);
    add_LOL(v.data(), D);
    write_gen_file(v.data(), name, type);
}
void CubicScalarGrid::compute_ELF(std::shared_ptr<Matrix> D, const std::string& name, const std::string& type) {
    auto v = std::vector<double>(npoints_, 0);
    add_ELF(v.data(), D);
    write_gen_file(v.data(), name, type);
}
std::pair<double, double> CubicScalarGrid::compute_isocontour_range(double* v2, double exponent) {
    double cumulative_threshold = options_.get_double("CUBEPROP_ISOCONTOUR_THRESHOLD");

    // Store the points with their weights and compute the sum of weights
    double sum_weight = 0.0;
    std::vector<std::pair<double, double>> sorted_points(npoints_);
    for (size_t ind = 0; ind < npoints_; ind++) {
        double value = v2[ind];
        double weight = std::pow(std::fabs(value), exponent);
        sum_weight += weight;
        sorted_points[ind] = std::make_pair(weight, value);
    }

    // Sort the points
    std::sort(sorted_points.rbegin(), sorted_points.rend());

    // Determine the positive and negative bounds
    double sum = 0.0;
    double negative_isocontour = 0.0;
    double positive_isocontour = 0.0;
    for (const std::pair<double, double>& p : sorted_points) {
        double value = p.second;
        if (value >= 0.0) {
            positive_isocontour = value;
        } else {
            negative_isocontour = value;
        }
        sum += p.first / sum_weight;
        if (sum > cumulative_threshold) break;
    }
    return std::make_pair(positive_isocontour, negative_isocontour);
}
std::string CubicScalarGrid::ecp_header() {
    std::stringstream ecp_head;
    ecp_head.str("");
    if (primary_->has_ECP()) {
        ecp_head << " Total core: " << primary_->n_ecp_core() << " [e] from 1-indexed atoms (";
        std::stringstream ecp_atoms;
        std::stringstream ecp_ncore;
        for (int A = 0; A < mol_->natom(); A++) {
            if (primary_->n_ecp_core(mol_->label(A)) > 0) {
                ecp_atoms << A+1 << "[" << mol_->symbol(A) << "], ";
                ecp_ncore << primary_->n_ecp_core(mol_->label(A)) << ", ";
            }
        }
        ecp_head << ecp_atoms.str().substr(0,ecp_atoms.str().length()-2) << ") electrons ("
                 << ecp_ncore.str().substr(0,ecp_ncore.str().length()-2) << ").";
    }
    return ecp_head.str();
}
}  // namespace psi
