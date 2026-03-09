/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2025 The Psi4 Developers.
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

#include "psi4/libciomr/libciomr.h"
#include "psi4/libmints/cdsalclist.h"
#include "psi4/libmints/potential.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/sobasis.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "typedefs.h"

#include <libint2/engine.h>

using namespace psi;

PotentialInt::PotentialInt(std::vector<SphericalTransform> &st, std::shared_ptr<BasisSet> bs1,
                           std::shared_ptr<BasisSet> bs2, int deriv)
    : OneBodyAOInt(st, bs1, bs2, deriv) {
    // When computing potential integrals with different bra and ket basis, the atom definition is taken from the bra basis.

    int max_am = std::max(basis1()->max_am(), basis2()->max_am());
    int max_nprim = std::max(basis1()->max_nprimitive(), basis2()->max_nprimitive());

    // Setup the initial field of partial charges
    std::vector<std::pair<double, std::array<double, 3>>> params;
    for (int A = 0; A < bs1_->molecule()->natom(); A++) {
        params.push_back({
            (double)bs1_->molecule()->Z(A),
            {bs1_->molecule()->x(A), bs1_->molecule()->y(A), bs1_->molecule()->z(A)}});
    }

    engine0_ = std::make_unique<libint2::Engine>(libint2::Operator::nuclear, max_nprim, max_am, 0);
    engine0_->set_params(params);

    if (deriv == 1) {
        const auto nresults = 3 * (2 + bs1_->molecule()->natom());

        set_chunks(nresults);

        engine1_ = std::make_unique<libint2::Engine>(libint2::Operator::nuclear, max_nprim, max_am, 1);
        engine1_->set_params(params);

    } else if (deriv == 2) {
        constexpr auto nopers = libint2::operator_traits<libint2::Operator::nuclear>::nopers;
        const auto nresults = nopers * libint2::num_geometrical_derivatives(bs1_->molecule()->natom(), 2);

        set_chunks(nresults);

        engine1_ = std::make_unique<libint2::Engine>(libint2::Operator::nuclear, max_nprim, max_am, 1);
        engine1_->set_params(params);
        engine2_ = std::make_unique<libint2::Engine>(libint2::Operator::nuclear, max_nprim, max_am, 2);
        engine2_->set_params(params);
    } else if (deriv > 2) {
        throw PSIEXCEPTION("PotentialInt only supports derivatives <= 2.");
    }

    buffer_ = nullptr;
    buffers_.resize(nchunk_);
}

PotentialInt::~PotentialInt() {}



void PotentialInt::set_charge_field(const std::vector<std::pair<double, std::array<double, 3>>>& Zxyz) {
    engine0_->set_params(Zxyz);
    if (engine1_) engine1_->set_params(Zxyz);
    if (engine2_) engine2_->set_params(Zxyz);
    Zxyz_ = Zxyz;
}


PotentialSOInt::PotentialSOInt(const std::shared_ptr<OneBodyAOInt> &aoint, const std::shared_ptr<IntegralFactory> &fact)
    : OneBodySOInt(aoint, fact) {
    natom_ = ob_->basis1()->molecule()->natom();
}

PotentialSOInt::PotentialSOInt(const std::shared_ptr<OneBodyAOInt> &aoint, const IntegralFactory *fact)
    : OneBodySOInt(aoint, fact) {
    natom_ = ob_->basis1()->molecule()->natom();
}

void PotentialSOInt::compute_deriv1(std::vector<SharedMatrix> result, const CdSalcList &cdsalcs) {

    // Do some checks:
    if (deriv_ < 1)
        throw SanityCheckError("OneBodySOInt::compute_deriv1: integral object not created to handle derivatives.",
                               __FILE__, __LINE__);

    if (result.size() != cdsalcs.ncd())
        throw SanityCheckError("OneBodySOInt::compute_deriv1: result vector size does not match SALC size.", __FILE__,
                               __LINE__);

    int ns1 = b1_->nshell();
    int ns2 = b2_->nshell();

    // Loop over unique SO shells.
    for (int ish = 0; ish < ns1; ++ish) {
        const SOTransform &t1 = b1_->sotrans(ish);
        int nao1 = b1_->naofunction(ish);

        for (int jsh = 0; jsh < ns2; ++jsh) {
            const SOTransform &t2 = b2_->sotrans(jsh);
            int nao2 = b2_->naofunction(jsh);

            int nao12 = nao1 * nao2;

            // loop through the AO shells that make up this SO shell
            for (int i = 0; i < t1.naoshell; ++i) {
                const SOTransformShell &s1 = t1.aoshell[i];

                for (int j = 0; j < t2.naoshell; ++j) {
                    const SOTransformShell &s2 = t2.aoshell[j];

                    ob_->compute_shell_deriv1(s1.aoshell, s2.aoshell);
                    const auto &buffers = ob_->buffers(); 
                    size_t nchunks = buffers.size();
                    int icenter = b1_->basis()->shell(s1.aoshell).ncenter();
                    int jcenter = b2_->basis()->shell(s2.aoshell).ncenter();

                    for (int itr = 0; itr < s1.nfunc(); ++itr) {
                        const SOTransformFunction &ifunc = s1.func[itr];
                        // SO transform coefficient
                        double icoef = ifunc.coef;
                        // AO function offset in a linear array
                        int iaofunc = ifunc.aofunc;
                        // SO function offset in a linear array
                        int isofunc = b1_->function_offset_within_shell(ish, ifunc.irrep) + ifunc.sofunc;
                        // AO function offset in a linear array
                        int iaooff = iaofunc;
                        // Relative position of the SO function within its irrep
                        int irel = b1_->function_within_irrep(ish, isofunc);
                        int iirrep = ifunc.irrep;

                        for (int jtr = 0; jtr < s2.nfunc(); ++jtr) {
                            const SOTransformFunction &jfunc = s2.func[jtr];
                            double jcoef = jfunc.coef * icoef;
                            int jaofunc = jfunc.aofunc;
                            int jsofunc = b2_->function_offset_within_shell(jsh, jfunc.irrep) + jfunc.sofunc;
                            int jaooff = iaooff * nao2 + jaofunc;
                            int jrel = b2_->function_within_irrep(jsh, jsofunc);
                            int jirrep = jfunc.irrep;

                            for (int chunk = 0; chunk < nchunks; ++chunk) {
                                int atom = (chunk < 1 ? icenter : (chunk < 2 ? jcenter : chunk-2) );
                                const CdSalcWRTAtom &cdsalc1 = cdsalcs.atom_salc(atom);

                                double jcoef_aobuf = jcoef * buffers[3*chunk + 0][jaooff];
                                for (int nx = 0; nx < cdsalc1.nx(); ++nx) {
                                    const CdSalcWRTAtom::Component element = cdsalc1.x(nx);
                                    double temp = jcoef_aobuf * element.coef;
                                    if ((iirrep ^ jirrep) == element.irrep && std::fabs(temp) > 1.0e-10) {
                                        result[element.salc]->add(iirrep, irel, jrel, temp);
                                    }
                                }

                                jcoef_aobuf = jcoef * buffers[3*chunk+1][jaooff];
                                for (int ny = 0; ny < cdsalc1.ny(); ++ny) {
                                    const CdSalcWRTAtom::Component element = cdsalc1.y(ny);
                                    double temp = jcoef_aobuf * element.coef;
                                    if ((iirrep ^ jirrep) == element.irrep && std::fabs(temp) > 1.0e-10) {
                                        result[element.salc]->add(iirrep, irel, jrel, temp);
                                    }
                                }

                                jcoef_aobuf = jcoef * buffers[3*chunk+2][jaooff];
                                for (int nz = 0; nz < cdsalc1.nz(); ++nz) {
                                    const CdSalcWRTAtom::Component element = cdsalc1.z(nz);
                                    double temp = jcoef_aobuf * element.coef;
                                    if ((iirrep ^ jirrep) == element.irrep && std::fabs(temp) > 1.0e-10) {
                                        result[element.salc]->add(iirrep, irel, jrel, temp);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}


