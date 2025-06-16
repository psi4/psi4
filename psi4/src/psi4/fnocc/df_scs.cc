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

#include "psi4/psi4-dec.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libqt/qt.h"
#include "psi4/libciomr/libciomr.h"
#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_wtime() 0.0
#endif

#include "blas.h"
#include "ccsd.h"
#include "psi4/libmints/basisset.h"
#include "psi4/lib3index/3index.h"

using namespace psi;

namespace psi {
namespace fnocc {

std::tuple<double, double, SharedMatrix, SharedMatrix> DFCoupledCluster::ComputePair(const std::string& name) {
    long int v = nvirt;
    long int o = ndoccact;
    long int rs = nmo;

    double ssenergy = 0.0;
    double osenergy = 0.0;

    // df (ia|bj) formerly E2klcd
    F_DGEMM('n', 't', o * v, o * v, nQ, 1.0, Qov, o * v, Qov, o * v, 0.0, integrals, o * v);

    if (t2_on_disk) {
        auto psio = std::make_shared<PSIO>();
        psio->open(PSIF_DCC_T2, PSIO_OPEN_OLD);
        psio->read_entry(PSIF_DCC_T2, "t2", (char*)&tempv[0], o * o * v * v * sizeof(double));
        psio->close(PSIF_DCC_T2, 1);
        tb = tempv;
    }

    auto matAA = std::make_shared<Matrix>(name + " Alpha-Alpha Pair Energies", o, o);
    auto matAB = std::make_shared<Matrix>(name + " Alpha-Beta Pair Energies", o, o);

    // The sum over i, j gives the "Coulomb-like" part of the (i, j) energy and the "Exchnge-like"
    // part of the (j, i) pair energy. For this reason, we need sum over both i, j and also add to the
    // i, j and j, i pair energies.
    for (long int i = 0; i < o; i++) {
        for (long int j = 0; j < o; j++) {
            double pair_os = 0;
            double pair_ss = 0;
            for (long int a = 0; a < v; a++) {
                for (long int b = 0; b < v; b++) {
                    long int ijab = a * v * o * o + b * o * o + i * o + j;
                    long int iajb = i * v * v * o + a * v * o + j * v + b;
                    long int jaib = iajb + (i - j) * v * (1 - v * o);

                    pair_os += integrals[iajb] * (tb[ijab] + t1[a * o + i] * t1[b * o + j]);
                    pair_ss += integrals[iajb] * (tb[ijab] - tb[b * o * o * v + a * o * o + i * o + j]);
                    pair_ss += integrals[iajb] *
                                (t1[a * o + i] * t1[b * o + j] - t1[b * o + i] * t1[a * o + j]);
                }
            }
            osenergy += pair_os;
            ssenergy += pair_ss;
            matAA->add(i, j, pair_ss);
            matAB->add(i, j, pair_os);
            if (i != j) {
                matAA->add(j, i, pair_ss);
                matAB->add(j, i, pair_os);
            }
        }
    }

    return std::make_tuple(osenergy, ssenergy, matAA, matAB);
}

std::tuple<SharedMatrix, SharedMatrix> spin_adapt(SharedMatrix ss, SharedMatrix os) {
    auto triplet = ss->clone();
    triplet->scale(1.5);
    auto singlet = os->clone();
    singlet->scale(2);
    singlet->axpy(-0.5, ss);
    for (int i = 0; i < singlet->nrow(); i++) {
        singlet->set(i, i, singlet->get(i, i) / 2);
    }
    return std::make_tuple(singlet, triplet);
}

void DFCoupledCluster::SCS_CCSD() {
    SharedMatrix CCA, CCB;
    std::tie(eccsd_os, eccsd_ss, CCA, CCB) = ComputePair("CC");
    eccsd = eccsd_os + eccsd_ss;
    CCA->add(delta_pair_energies_ss);
    CCB->add(delta_pair_energies_os);
    set_array_variable("CC ALPHA-ALPHA PAIR ENERGIES", CCA);
    set_array_variable("CC ALPHA-BETA PAIR ENERGIES", CCB);
    set_array_variable("CCSD ALPHA-ALPHA PAIR ENERGIES", CCA);
    set_array_variable("CCSD ALPHA-BETA PAIR ENERGIES", CCB);
    SharedMatrix singlet, triplet;
    std::tie(singlet, triplet) = spin_adapt(CCA, CCB);
    singlet->set_name("CC SINGLET PAIR ENERGIES");
    set_array_variable("CC SINGLET PAIR ENERGIES", singlet);
    set_array_variable("CCSD SINGLET PAIR ENERGIES", singlet);
    triplet->set_name("CC TRIPLET PAIR ENERGIES");
    set_array_variable("CC TRIPLET PAIR ENERGIES", triplet);
    set_array_variable("CCSD TRIPLET PAIR ENERGIES", triplet);
}

void DFCoupledCluster::SCS_MP2() {
    SharedMatrix MPA, MPB;
    std::tie(emp2_os, emp2_ss, MPA, MPB) = ComputePair("MP2");
    emp2 = emp2_os + emp2_ss;
    if (has_array_variable("MP2 ALPHA-ALPHA PAIR ENERGIES")) {
        // We just computed the truncated MP2 pair energies.
        delta_pair_energies_ss = array_variable("MP2 ALPHA-ALPHA PAIR ENERGIES")->clone();
        delta_pair_energies_os = array_variable("MP2 ALPHA-BETA PAIR ENERGIES")->clone();
        delta_pair_energies_ss->subtract(MPA);
        delta_pair_energies_os->subtract(MPB);
    } else {
        // We just computed the true MP2 pair energies.
        set_array_variable("MP2 ALPHA-ALPHA PAIR ENERGIES", MPA);
        set_array_variable("MP2 ALPHA-BETA PAIR ENERGIES", MPB);
        delta_pair_energies_ss = std::make_shared<Matrix>("Same-Spin Pair Energies", MPA->rowspi(), MPA->colspi());
        delta_pair_energies_os = std::make_shared<Matrix>("Opposite-Spin Pair Energies", MPB->rowspi(), MPB->colspi());
        SharedMatrix singlet, triplet;
        std::tie(singlet, triplet) = spin_adapt(MPA, MPB);
        singlet->set_name("MP2 SINGLET PAIR ENERGIES");
        set_array_variable("MP2 SINGLET PAIR ENERGIES", singlet);
        triplet->set_name("MP2 TRIPLET PAIR ENERGIES");
        set_array_variable("MP2 TRIPLET PAIR ENERGIES", triplet);
    }
}
}
}
