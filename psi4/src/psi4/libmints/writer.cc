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
#include <cstdio>
#include <utility>
#include <algorithm>
#include <vector>
#include <libint2/config.h>
#include "psi4/libmints/writer.h"
#include "psi4/psi4-dec.h"
#include "psi4/physconst.h"
#include "psi4/masses.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/pointgrp.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/petitelist.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libqt/qt.h"

using namespace psi;
;


MoldenWriter::MoldenWriter(std::shared_ptr<Wavefunction> wavefunction) : wavefunction_(wavefunction) {
    outfile->Printf("\tConstructing a MoldenWriter and then calling write instead of using `wfn.write_molden(name)`\n");
    outfile->Printf("\tis both buggy and deprecated, and as soon as 1.5 it will stop working.\n\n");
}

void MoldenWriter::write(const std::string &filename, std::shared_ptr<Matrix> Ca, std::shared_ptr<Matrix> Cb,
                         std::shared_ptr<Vector> Ea, std::shared_ptr<Vector> Eb, std::shared_ptr<Vector> OccA,
                         std::shared_ptr<Vector> OccB, bool dovirtual) {
    auto mode = std::ostream::trunc;
    auto printer = std::make_shared<PsiOutStream>(filename, mode);

    int atom;

    printer->Printf("[Molden Format]\n");

    // Get the molecule for ease
    BasisSet &basisset = *wavefunction_->basisset().get();
    Molecule &mol = *basisset.molecule().get();

    //    basisset.print_detail();

    // Print the molecule to molden
    printer->Printf("[Atoms] (AU)\n");
    for (atom = 0; atom < mol.natom(); ++atom) {
        Vector3 coord = mol.xyz(atom);
        printer->Printf("%-2s  %2d  %3d   %20.10f %20.10f %20.10f\n", mol.symbol(atom).c_str(), atom + 1,
                        static_cast<int>(mol.Z(atom)), coord[0], coord[1], coord[2]);
    }

    // Dump the basis set using code adapted from psi2molden
    printer->Printf("[GTO]\n");

    // For each atom
    for (atom = 0; atom < mol.natom(); ++atom) {
        printer->Printf("  %d 0\n", atom + 1);

        // Go through all the shells on this center
        for (int shell = 0; shell < basisset.nshell_on_center(atom); ++shell) {
            int overall_shell = basisset.shell_on_center(atom, shell);

            const GaussianShell &gs = basisset.shell(overall_shell);

            printer->Printf(" %c%5d  1.00\n", gs.amchar(), gs.nprimitive());

            for (int prim = 0; prim < gs.nprimitive(); ++prim) {
                printer->Printf("%20.10f %20.10f\n", gs.exp(prim), gs.original_coef(prim));
            }
        }

        // An empty line separates atoms
        printer->Printf("\n");
    }

    // Convert Ca & Cb
    auto pl = std::make_shared<PetiteList>(wavefunction_->basisset(), wavefunction_->integral());
    // get the "aotoso" transformation matrix, ao by so
    SharedMatrix aotoso = pl->aotoso();
    // need dimensions
    const Dimension aos = pl->AO_basisdim();
    const Dimension sos = pl->SO_basisdim();
    const Dimension nmo = Ca->colspi();

    auto Ca_ao_mo = std::make_shared<Matrix>("Ca AO x MO", aos, nmo);
    auto Cb_ao_mo = std::make_shared<Matrix>("Cb AO x MO", aos, nmo);

    // do the half transform
    Ca_ao_mo->gemm(false, false, 1.0, aotoso, Ca, 0.0);
    Cb_ao_mo->gemm(false, false, 1.0, aotoso, Cb, 0.0);

    //    aotoso->print();
    //    Ca_ao_mo->print();
    //    Cb_ao_mo->print();

    // The order Molden expects
    //     P: x, y, z
    //    5D: D 0, D+1, D-1, D+2, D-2
    //    6D: xx, yy, zz, xy, xz, yz
    //
    //    7F: F 0, F+1, F-1, F+2, F-2, F+3, F-3
    //   10F: xxx, yyy, zzz, xyy, xxy, xxz, xzz, yzz, yyz, xyz
    //
    //    9G: G 0, G+1, G-1, G+2, G-2, G+3, G-3, G+4, G-4
    //   15G: xxxx yyyy zzzz xxxy xxxz yyyx yyyz zzzx zzzy,
    //        xxyy xxzz yyzz xxyz yyxz zzxy
    // Since Molden doesn't handle higher than g we'll just leave them be.
    int molden_cartesian_order[][15] = {
        {2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},      // p
        {0, 3, 4, 1, 5, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0},      // d
        {0, 4, 5, 3, 9, 6, 1, 8, 7, 2, 0, 0, 0, 0, 0},      // f
        {0, 3, 4, 9, 12, 10, 5, 13, 14, 7, 1, 6, 11, 8, 2}  // g
    };

    int nirrep = Ca_ao_mo->nirrep();
    Dimension countpi(nirrep);
    Dimension zeropi(nirrep);
    Dimension ncartpi(nirrep);

    for (int i = 0; i < basisset.nshell(); i++) {
        int am = basisset.shell(i).am();

        int ncart = basisset.shell(i).nfunction();
        if ((am == 1 && basisset.has_puream()) || (am > 1 && am < 5 && basisset.shell(i).is_cartesian())) {
            for (int h = 0; h < nirrep; ++h) ncartpi[h] = ncart;

            Slice row_slice(countpi, countpi + ncartpi);
            Slice acol_slice(zeropi, zeropi + Ca_ao_mo->colspi());
            Slice bcol_slice(zeropi, zeropi + Cb_ao_mo->colspi());
            SharedMatrix temp_a = Ca_ao_mo->get_block(row_slice, acol_slice);
            SharedMatrix temp_b = Cb_ao_mo->get_block(row_slice, bcol_slice);

            for (int j = 0; j < ncart; j++) {
                for (int h = 0; h < Ca_ao_mo->nirrep(); ++h) {
                    for (int k = 0; k < Ca_ao_mo->coldim(h); ++k) {
                        // outfile->Printf( "am %d\n, from %d to %d\n", am, j, countpi[h] +
                        // molden_cartesian_order[am-1][j]);
                        Ca_ao_mo->set(h, countpi[h] + molden_cartesian_order[am - 1][j], k, temp_a->get(h, j, k));
                        Cb_ao_mo->set(h, countpi[h] + molden_cartesian_order[am - 1][j], k, temp_b->get(h, j, k));
                    }
                }
            }
        }

        for (int h = 0; h < nirrep; ++h) countpi[h] += ncart;
    }

    if (basisset.has_puream()) {
        // Tell Molden to use spherical.  5d implies 5d and 7f.
        printer->Printf("[5D]\n[9G]\n\n");
    }
    CharacterTable ct = mol.point_group()->char_table();

    // Dump MO's to the molden file
    printer->Printf("[MO]\n");

    std::vector<std::pair<double, std::pair<int, int> > > mos;

    // Number of MOs to write
    std::vector<int> nmoh(wavefunction_->nirrep());
    for (int h = 0; h < wavefunction_->nirrep(); ++h) {
        if (dovirtual)
            nmoh[h] = wavefunction_->nmopi()[h];
        else
            nmoh[h] = wavefunction_->doccpi()[h] + wavefunction_->soccpi()[h];
    }

    // do alpha's
    bool SameOcc = true;
    for (int h = 0; h < wavefunction_->nirrep(); ++h) {
        for (int n = 0; n < nmoh[h]; ++n) {
            mos.push_back(std::make_pair(Ea->get(h, n), std::make_pair(h, n)));
            if (std::fabs(OccA->get(h, n) - OccB->get(h, n)) > 1e-10) SameOcc = false;
        }
    }
    std::sort(mos.begin(), mos.end());

    for (int i = 0; i < (int)mos.size(); ++i) {
        int h = mos[i].second.first;
        int n = mos[i].second.second;

        printer->Printf(" Sym= %s\n", ct.gamma(h).symbol());
        printer->Printf(" Ene= %24.10e\n", Ea->get(h, n));
        printer->Printf(" Spin= Alpha\n");
        if (Ca == Cb && Ea == Eb && SameOcc)
            printer->Printf(" Occup= %24.10e\n", OccA->get(h, n) + OccB->get(h, n));
        else
            printer->Printf(" Occup= %24.10e\n", OccA->get(h, n));
        for (int so = 0; so < wavefunction_->nso(); ++so)
            printer->Printf("%3d %24.10e\n", so + 1, Ca_ao_mo->get(h, so, n));
    }

    // do beta's
    mos.clear();
    if (Ca != Cb || Ea != Eb || !SameOcc) {
        for (int h = 0; h < wavefunction_->nirrep(); ++h) {
            for (int n = 0; n < nmoh[h]; ++n) {
                mos.push_back(std::make_pair(Eb->get(h, n), std::make_pair(h, n)));
            }
        }
        std::sort(mos.begin(), mos.end());

        for (int i = 0; i < (int)mos.size(); ++i) {
            int h = mos[i].second.first;
            int n = mos[i].second.second;

            printer->Printf(" Sym= %s\n", ct.gamma(h).symbol());
            printer->Printf(" Ene= %24.10e\n", Eb->get(h, n));
            printer->Printf(" Spin= Beta\n");
            printer->Printf(" Occup= %24.10e\n", OccB->get(h, n));
            for (int so = 0; so < wavefunction_->nso(); ++so)
                printer->Printf("%3d %24.10e\n", so + 1, Cb_ao_mo->get(h, so, n));
        }
    }
}

FCHKWriter::FCHKWriter(std::shared_ptr<Wavefunction> wavefunction) : wavefunction_(wavefunction) {}

void FCHKWriter::write_number(const char *label, double value) { fprintf(chk_, "%-43sR%27.15e\n", label, value); }

void FCHKWriter::write_number(const char *label, int value) { fprintf(chk_, "%-43sI%17d\n", label, value); }

void FCHKWriter::write_sym_matrix(const char *label, const SharedMatrix &mat) {
    int dim = mat->rowdim();
    fprintf(chk_, "%-43s%-3s N=%12d\n", label, "R", (dim * dim + dim) / 2);

    int count = 0;
    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j <= i; ++j) {
            fprintf(chk_, "%16.8e", mat->get(i, j));
            if (count % 5 == 4) fprintf(chk_, "\n");
            ++count;
        }
    }
    if (count % 5) fprintf(chk_, "\n");
}

void FCHKWriter::write_matrix(const char *label, const SharedVector &mat) {
    int dim = mat->dim();
    fprintf(chk_, "%-43s%-3s N=%12d\n", label, "R", dim);

    int count = 0;
    for (int i = 0; i < dim; ++i) {
        fprintf(chk_, "%16.8e", mat->get(count));
        if (count % 5 == 4) fprintf(chk_, "\n");
        ++count;
    }
    if (count % 5) fprintf(chk_, "\n");
}

void FCHKWriter::write_matrix(const char *label, const SharedMatrix &mat) {
    int rowdim = mat->rowdim();
    int coldim = mat->coldim();
    fprintf(chk_, "%-43s%-3s N=%12d\n", label, "R", rowdim * coldim);

    int count = 0;
    for (int i = 0; i < rowdim; ++i) {
        for (int j = 0; j < coldim; ++j) {
            fprintf(chk_, "%16.8e", mat->get(i, j));
            if (count % 5 == 4) fprintf(chk_, "\n");
            ++count;
        }
    }
    if (count % 5) fprintf(chk_, "\n");
}

void FCHKWriter::write_matrix(const char *label, const std::vector<double> &mat) {
    int dim = mat.size();
    fprintf(chk_, "%-43s%-3s N=%12d\n", label, "R", dim);

    int count = 0;
    for (int i = 0; i < dim; ++i) {
        fprintf(chk_, "%16.8e", mat[count]);
        if (count % 5 == 4) fprintf(chk_, "\n");
        ++count;
    }
    if (count % 5) fprintf(chk_, "\n");
}

void FCHKWriter::write_matrix(const char *label, const std::vector<int> &mat) {
    int dim = mat.size();
    fprintf(chk_, "%-43s%-3s N=%12d\n", label, "I", dim);

    int count = 0;
    for (int i = 0; i < dim; ++i) {
        fprintf(chk_, "%12d", mat[count]);
        if (count % 6 == 5) fprintf(chk_, "\n");
        ++count;
    }
    if (count % 6) fprintf(chk_, "\n");
}

void FCHKWriter::set_postscf_density_label(const std::string &label) {
    postscf_density_label_ = ("Total" + label);
    spin_postscf_density_label_ = ("Spin" + label);
}

void FCHKWriter::write(const std::string &filename) {
    chk_ = fopen(filename.c_str(), "w");
    std::shared_ptr<BasisSet> basis = wavefunction_->basisset();
    int maxam = basis->max_am();
    std::shared_ptr<Molecule> mol = wavefunction_->molecule();
    std::shared_ptr<Wavefunction> refwf = wavefunction_->reference_wavefunction();

    // Orbitals
    Ca_ao = wavefunction_->Ca_subset("AO");
    Cb_ao = wavefunction_->Cb_subset("AO");

    const std::string &name = wavefunction_->name();
    const std::string &basisname = basis->name();
    int nbf = basis->nbf();
    int nmo = Ca_ao->ncol();
    int nalpha = wavefunction_->nalpha();
    int nbeta = wavefunction_->nbeta();
    int natoms = mol->natom();
    int nprimitive = basis->nprimitive();

    // SCF density matrices
    SharedMatrix Da_ao;
    SharedMatrix Db_ao;
    // Post-Hartree-Fock density matrices
    SharedMatrix DPHFa_ao;
    SharedMatrix DPHFb_ao;

    // Post Hartree Fock?
    bool pHF = (refwf != NULL);

    if (pHF) {
        Da_ao = refwf->Da_subset("AO");
        Db_ao = refwf->Db_subset("AO");
        DPHFa_ao = wavefunction_->Da_subset("AO");
        DPHFb_ao = wavefunction_->Db_subset("AO");
    } else {
        // SCF level of theory
        Da_ao = wavefunction_->Da_subset("AO");
        Db_ao = wavefunction_->Db_subset("AO");
    }

    // Total and spin density
    SharedMatrix Dspin_ao;
    SharedMatrix DPHFtot_ao;
    SharedMatrix DPHFspin_ao;

    Dtot_ao = Da_ao->clone();
    Dtot_ao->add(Db_ao);
    Dspin_ao = Da_ao->clone();
    Dspin_ao->subtract(Db_ao);

    if (pHF) {
        DPHFtot_ao = DPHFa_ao->clone();
        DPHFtot_ao->add(DPHFb_ao);
        DPHFspin_ao = DPHFa_ao->clone();
        DPHFspin_ao->subtract(DPHFb_ao);
    }

    std::vector<double> coords;
    std::vector<double> nuc_charges;
    std::vector<int> atomic_numbers;
    std::vector<int> int_atomic_weights;
    std::vector<double> atomic_weights;
    double to_bohr = mol->units() == Molecule::Angstrom ? 1.0 / pc_bohr2angstroms : 1.0;
    for (int atom = 0; atom < natoms; ++atom) {
        int intZ = static_cast<int>(mol->Z(atom));
        atomic_weights.push_back(mol->mass(atom));
        int_atomic_weights.push_back(mol->mass_number(atom));
        nuc_charges.push_back(mol->Z(atom));
        atomic_numbers.push_back(intZ > 0 ? mol->true_atomic_number(atom) : intZ); // care about ECP & ghosts!
        const Vector3 &xyz = mol->xyz(atom);
        coords.push_back(xyz[0]);
        coords.push_back(xyz[1]);
        coords.push_back(xyz[2]);
    }
    // For Cartesian functions we need to add a basis function normalization
    // constant of
    //      _______________________________
    //     / (2lx-1)!! (2ly-1)!! (2lz-1)!!
    //    /  -----------------------------
    //  \/             (2l-1)!!
    //
    // which is omitted in the CCA standard, adopted by Psi4.
    //
    // We also need to order basis functions to the Gaussian / GAMESS convention.
    // * When psi4_SHGAUSS_ORDERING=gaussian (usual case), spherical
    // harmonics are already defined appropriately, with the exception of
    // the fact that p functions are ordered 0 +1 -1, which is Z X Y, but
    // the FCHK format calls for X Y Z; this is a simple reordering operation.
    // * When psi4_SHGAUSS_ORDERING=standard (unusual but configurable), all
    // spherical harmonics need reordering.

    const double pureP_from_gss[3][3] = {
        //           0    1    2
        // Psi4:     Z    X    Y
        // Psi4:     0   +1   -1
        // Expected:+1   -1    0
        // Expected: X    Y    Z
        /* 0 */ {0.0, 1.0, 0.0},
        /* 1 */ {0.0, 0.0, 1.0},
        /* 2 */ {1.0, 0.0, 0.0},
    };

    const double pureP_from_sss[3][3] = {
        //            0    1    2
        // Psi4:     -1    0   +1
        // Expected: +1   -1    0
        /* 0 */ {0.0, 0.0, 1.0},
        /* 1 */ {1.0, 0.0, 0.0},
        /* 2 */ {0.0, 1.0, 0.0}
    };
    const double pureD_from_sss[5][5] = {
        //            0    1    2    3    4
        // Psi4:     -2   -1    0   +1   +2
        // Expected:  0   +1   -1   +2   -2
        /* 0 */ {0.0, 0.0, 1.0, 0.0, 0.0},
        /* 1 */ {0.0, 0.0, 0.0, 1.0, 0.0},
        /* 2 */ {0.0, 1.0, 0.0, 0.0, 0.0},
        /* 3 */ {0.0, 0.0, 0.0, 0.0, 1.0},
        /* 4 */ {1.0, 0.0, 0.0, 0.0, 0.0}
    };
    const double pureF_from_sss[7][7] = {
        //            0    1    2    3    4    5    6
        // Psi4:     -3   -2   -1    0   +1   +2   +3
        // Expected:  0   +1   -1   +2   -2   +3   -3
        /* 0 */ {0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0},
        /* 1 */ {0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0},
        /* 2 */ {0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0},
        /* 3 */ {0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0},
        /* 4 */ {0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        /* 5 */ {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0},
        /* 6 */ {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
    };
    const double pureG_from_sss[9][9] = {
        //            0    1    2    3    4    5    6    7    8
        // Psi4:     -4   -3   -2   -1    0   +1   +2   +3   +4
        // Expected:  0   +1   -1   +2   -2   +3   -3   +4   -4
        /* 0 */ {0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0},
        /* 1 */ {0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0},
        /* 2 */ {0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        /* 3 */ {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0},
        /* 4 */ {0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        /* 5 */ {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0},
        /* 6 */ {0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        /* 7 */ {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0},
        /* 8 */ {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
    };

    double pf1, pf2, pf3, pf4;
    pf1 = 1.0;              // aa
    pf2 = sqrt(1.0 / 3.0);  // ab
    const double cartD[6][6] = {
        //             0    1    2    3    4    5
        // Psi4:      XX   XY   XZ   YY   YZ   ZZ
        // Expected:  XX   YY   ZZ   XY   XZ   YZ
        /* 0 */ {pf1, 0.0, 0.0, 0.0, 0.0, 0.0},
        /* 1 */ {0.0, 0.0, 0.0, pf1, 0.0, 0.0},
        /* 2 */ {0.0, 0.0, 0.0, 0.0, 0.0, pf1},
        /* 3 */ {0.0, pf2, 0.0, 0.0, 0.0, 0.0},
        /* 4 */ {0.0, 0.0, pf2, 0.0, 0.0, 0.0},
        /* 5 */ {0.0, 0.0, 0.0, 0.0, pf2, 0.0},
    };

    pf1 = 1.0;               // aaa
    pf2 = sqrt(1.0 / 5.0);   // aab
    pf3 = sqrt(1.0 / 15.0);  // abc
    const double cartF[10][10] = {
        //            0    1    2    3    4    5    6    7    8    9
        // Psi4:     XXX  XXY  XXZ  XYY  XYZ  XZZ  YYY  YYZ  YZZ  ZZZ
        // Expected: XXX  YYY  ZZZ  XYY  XXY  XXZ  XZZ  YZZ  YYZ  XYZ
        /* 0 */ {pf1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        /* 1 */ {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, pf1, 0.0, 0.0, 0.0},
        /* 2 */ {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, pf1},
        /* 3 */ {0.0, 0.0, 0.0, pf2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        /* 4 */ {0.0, pf2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        /* 5 */ {0.0, 0.0, pf2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        /* 6 */ {0.0, 0.0, 0.0, 0.0, 0.0, pf2, 0.0, 0.0, 0.0, 0.0},
        /* 7 */ {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, pf2, 0.0},
        /* 8 */ {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, pf2, 0.0, 0.0},
        /* 9 */ {0.0, 0.0, 0.0, 0.0, pf3, 0.0, 0.0, 0.0, 0.0, 0.0},
    };
    pf1 = 1.0;               // aaaa
    pf2 = sqrt(1.0 / 7.0);   // aaab
    pf3 = sqrt(3.0 / 35.0);  // aabb
    pf4 = sqrt(1.0 / 35.0);  // abcc
    const double cartG[15][15] = {
        //             0    1    2    3    4    5    6    7    8    9   10   11   12   13   14
        // Psi4:     XXXX XXXY XXXZ XXYY XXYZ XXZZ XYYY XYYZ XYZZ XZZZ YYYY YYYZ YYZZ YZZZ ZZZZ
        // Expected: XXXX YYYY ZZZZ XXXY XXXZ XYYY YYYZ XZZZ YZZZ XXYY XXZZ YYZZ XXYZ XYYZ XYZZ
        // This is the ordering that I would expect to generate...
        /*  0 */ {pf1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        /*  1 */ {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, pf1, 0.0, 0.0, 0.0, 0.0},
        /*  2 */ {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, pf1},
        /*  3 */ {0.0, pf2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        /*  4 */ {0.0, 0.0, pf2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        /*  5 */ {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, pf2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        /*  6 */ {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, pf2, 0.0, 0.0, 0.0},
        /*  7 */ {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, pf2, 0.0, 0.0, 0.0, 0.0, 0.0},
        /*  8 */ {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, pf2, 0.0},
        /*  9 */ {0.0, 0.0, 0.0, pf3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        /* 10 */ {0.0, 0.0, 0.0, 0.0, 0.0, pf3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        /* 11 */ {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, pf3, 0.0, 0.0},
        /* 12 */ {0.0, 0.0, 0.0, 0.0, pf4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        /* 13 */ {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, pf4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        /* 14 */ {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, pf4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        // ... however, I need to use the following ordering to match a FCHK file that I found for
        // QZ water. We will go with the former, because that's the order expected by GDMA.
        // I doubt anybody will ever use a Cartesian basis set with G functions, so this is
        // most likely a non-issue.  I'll leave this here in case anybody has problems in future (ACS).
        ///*  2 */  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, pf1 },
        ///*  8 */  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, pf2, 0.0 },
        ///* 11 */  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, pf3, 0.0, 0.0 },
        ///*  6 */  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, pf2, 0.0, 0.0, 0.0 },
        ///*  1 */  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, pf1, 0.0, 0.0, 0.0, 0.0 },
        ///*  7 */  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, pf2, 0.0, 0.0, 0.0, 0.0, 0.0 },
        ///* 14 */  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, pf4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
        ///* 13 */  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, pf4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
        ///*  5 */  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, pf2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
        ///* 10 */  { 0.0, 0.0, 0.0, 0.0, 0.0, pf3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
        ///* 12 */  { 0.0, 0.0, 0.0, 0.0, pf4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
        ///*  9 */  { 0.0, 0.0, 0.0, pf3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
        ///*  4 */  { 0.0, 0.0, pf2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
        ///*  3 */  { 0.0, pf2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
        ///*  0 */  { pf1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
    };

    auto transmat = std::make_shared<Matrix>("Reorder", nbf, nbf);
    transmat->identity();
    int offset = 0;
    for (int nshell = 0; nshell < basis->nshell(); ++nshell) {
        const GaussianShell &shell = basis->shell(nshell);
        int am = shell.am();
        int nfunc = shell.nfunction();
        if (basis->has_puream()) {
#if psi4_SHGSHELL_ORDERING == LIBINT_SHGSHELL_ORDERING_STANDARD
            if (am == 1) {
                for (int row = 0; row < 3; ++row)
                    for (int col = 0; col < 3; ++col) transmat->set(offset + row, offset + col, pureP_from_sss[row][col]);
            } else if (am == 2) {
                for (int row = 0; row < 5; ++row)
                    for (int col = 0; col < 5; ++col) transmat->set(offset + row, offset + col, pureD_from_sss[row][col]);
            } else if (am == 3) {
                for (int row = 0; row < 7; ++row)
                    for (int col = 0; col < 7; ++col) transmat->set(offset + row, offset + col, pureF_from_sss[row][col]);
            } else if (am == 4) {
                for (int row = 0; row < 9; ++row)
                    for (int col = 0; col < 9; ++col) transmat->set(offset + row, offset + col, pureG_from_sss[row][col]);
            } else if (am >= 5) {
                throw PSIEXCEPTION("The Psi4 FCHK writer only supports up to G shell (l=4) spherical functions");
            }
#elif psi4_SHGSHELL_ORDERING == LIBINT_SHGSHELL_ORDERING_GAUSSIAN
            // Spherical harmonics - everything is fine, apart from P orbitals
            if (am == 1) {
                for (int row = 0; row < 3; ++row)
                    for (int col = 0; col < 3; ++col) transmat->set(offset + row, offset + col, pureP_from_gss[row][col]);
            }
#else
#  error "unknown value of macro psi4_SHGSHELL_ORDERING"
#endif
        } else {
            // Cartesians - S and P orbitals are fine, but higher terms need reordering
            if (am == 2) {
                for (int row = 0; row < 6; ++row)
                    for (int col = 0; col < 6; ++col) transmat->set(offset + row, offset + col, cartD[row][col]);
            } else if (am == 3) {
                for (int row = 0; row < 10; ++row)
                    for (int col = 0; col < 10; ++col) transmat->set(offset + row, offset + col, cartF[row][col]);
            } else if (am == 4) {
                for (int row = 0; row < 15; ++row)
                    for (int col = 0; col < 15; ++col) transmat->set(offset + row, offset + col, cartG[row][col]);
            } else if (am >= 5) {
                throw PSIEXCEPTION("The Psi4 FCHK writer only supports up to G shell (l=4) cartesian functions");
            }
        }
        offset += nfunc;
    }

    SharedMatrix reorderedDt(Dtot_ao->clone());
    reorderedDt->back_transform(Dtot_ao, transmat);
    SharedMatrix reorderedDs(Dspin_ao->clone());
    reorderedDs->back_transform(Dspin_ao, transmat);

    SharedMatrix reorderedDPHFt;
    SharedMatrix reorderedDPHFs;
    if (pHF) {
        reorderedDPHFt = (DPHFtot_ao->clone());
        reorderedDPHFt->back_transform(DPHFtot_ao, transmat);
        reorderedDPHFs = (DPHFspin_ao->clone());
        reorderedDPHFs->back_transform(DPHFspin_ao, transmat);
    }

    auto reorderedCa = std::make_shared<Matrix>("Reordered Ca", Ca_ao->ncol(), Ca_ao->nrow());
    auto reorderedCb = std::make_shared<Matrix>("Reordered Cb", Cb_ao->ncol(), Cb_ao->nrow());
    reorderedCa->gemm(true, true, 1.0, Ca_ao, transmat, 0.0);
    reorderedCb->gemm(true, true, 1.0, Cb_ao, transmat, 0.0);
    for (int i = 0; i < reorderedDt->nrow(); ++i)
        for (int j = 0; j < reorderedDt->ncol(); ++j)
            if (std::fabs(reorderedDt->get(i, j)) < 1E-12) reorderedDt->set(i, j, 0.0);
    for (int i = 0; i < reorderedCa->nrow(); ++i)
        for (int j = 0; j < reorderedCa->ncol(); ++j)
            if (std::fabs(reorderedCa->get(i, j)) < 1E-12) reorderedCa->set(i, j, 0.0);
    for (int i = 0; i < reorderedCb->nrow(); ++i)
        for (int j = 0; j < reorderedCb->ncol(); ++j)
            if (std::fabs(reorderedCb->get(i, j)) < 1E-12) reorderedCb->set(i, j, 0.0);
    std::vector<double> shell_coords;
    std::vector<double> coefficients;
    std::vector<double> exponents;
    std::vector<int> prim_per_shell;
    std::vector<int> shell_to_atom;
    std::vector<int> shell_am;
    int nshell = basis->nshell();
    for (int shell = 0; shell < nshell; ++shell) {
        const GaussianShell &s = basis->shell(shell);
        int am = s.am();
        int nprim = s.nprimitive();
        prim_per_shell.push_back(nprim);
        int center = basis->shell_to_center(shell);
        shell_to_atom.push_back(center + 1);
        const Vector3 &xyz = mol->xyz(center);
        shell_coords.push_back(xyz[0]);
        shell_coords.push_back(xyz[1]);
        shell_coords.push_back(xyz[2]);
        int shell_type_fac = s.am() > 1 && s.is_pure() ? -1 : 1;
        shell_am.push_back(shell_type_fac * am);
        double normfac = 1.0;
        if (am > 1)
            // Undo the angular momentum normalization, applied by the ERD code
            normfac = sqrt(df[2 * am] / pow(2.0, 2.0 * am));
        for (int prim = 0; prim < nprim; ++prim) {
            exponents.push_back(s.exp(prim));
            double normfac2 = normfac / pow(s.exp(prim), 0.5 * ((double)am + 1.5));
            coefficients.push_back(normfac2 * s.erd_coef(prim));
        }
    }

    fprintf(chk_, "Generated by Psi4\n");
    std::string jobtype = wavefunction_->gradient() ? "FORCE" : "SP";
    fprintf(chk_, "%10s%30s%30s\n", jobtype.c_str(), name.c_str(), basisname.c_str());
    write_number("Number of atoms", natoms);
    write_number("Charge", mol->molecular_charge());
    write_number("Multiplicity", mol->multiplicity());
    write_number("Number of electrons", nalpha + nbeta);
    write_number("Number of alpha electrons", nalpha);
    write_number("Number of beta electrons", nbeta);
    write_number("Number of basis functions", nbf);
    write_number("Number of independent functions", nmo);
    write_matrix("Atomic numbers", atomic_numbers);
    write_matrix("Nuclear charges", nuc_charges);
    write_matrix("Current cartesian coordinates", coords);
    write_matrix("Integer atomic weights", int_atomic_weights);
    write_matrix("Real atomic weights", atomic_weights);
    write_number("Number of primitive shells", nprimitive);
    write_number("Number of contracted shells", nshell);
    write_number("Pure/Cartesian d shells", !basis->has_puream());
    write_number("Pure/Cartesian f shells", !basis->has_puream());
    write_number("Highest angular momentum", maxam);
    write_number("Largest degree of contraction", basis->max_nprimitive());
    write_matrix("Shell types", shell_am);
    write_matrix("Number of primitives per shell", prim_per_shell);
    write_matrix("Shell to atom map", shell_to_atom);
    write_matrix("Primitive exponents", exponents);
    write_matrix("Contraction coefficients", coefficients);
    write_matrix("Coordinates of each shell", shell_coords);
    write_number("Total Energy", wavefunction_->energy());

    // This is the format expected in the formatted checkpoint file
    write_matrix("Alpha Orbital Energies", wavefunction_->epsilon_a_subset("AO"));
    write_matrix("Alpha MO coefficients", reorderedCa);
    write_sym_matrix("Total SCF Density", reorderedDt);
    // labels for the densities are figured out at the python level
    if (pHF) write_sym_matrix(postscf_density_label_.c_str(), reorderedDPHFt);
    if (!wavefunction_->same_a_b_orbs() || !wavefunction_->same_a_b_dens()) {
        // These are only printed out if the orbitals or density is not spin-restricted
        write_matrix("Beta Orbital Energies", wavefunction_->epsilon_b_subset("AO"));
        write_matrix("Beta MO coefficients", reorderedCb);
        write_sym_matrix("Spin SCF Density", reorderedDs);
        if (pHF) write_sym_matrix(spin_postscf_density_label_.c_str(), reorderedDPHFs);
    }

    SharedMatrix gradient = wavefunction_->gradient();
    if (gradient) write_matrix("Cartesian Gradient", gradient);

    fclose(chk_);
    chk_ = 0;
}

MOWriter::MOWriter(std::shared_ptr<Wavefunction> wavefunction)
    : wavefunction_(wavefunction), restricted_(wavefunction->same_a_b_orbs()) {}

void MOWriter::write() {
    // Get the molecule for ease
    BasisSet &basisset = *wavefunction_->basisset().get();
    Molecule &mol = *basisset.molecule().get();

    // Convert Ca & Cb
    // make copies
    Matrix Ca(wavefunction_->Ca());
    Matrix Cb(wavefunction_->Cb());
    Vector &Ea = *wavefunction_->epsilon_a().get();
    Vector &Eb = *wavefunction_->epsilon_b().get();

    auto pl = std::make_shared<PetiteList>(wavefunction_->basisset(), wavefunction_->integral());

    // get the "aotoso" transformation matrix, ao by so
    SharedMatrix aotoso = pl->aotoso();
    // need dimensions
    const Dimension aos = pl->AO_basisdim();
    const Dimension sos = pl->SO_basisdim();
    const Dimension mos = wavefunction_->nmopi();

    auto Ca_ao_mo = std::make_shared<Matrix>("Ca AO x MO", aos, mos);
    auto Cb_ao_mo = std::make_shared<Matrix>("Cb AO x MO", aos, mos);

    // do the half transform
    Ca_ao_mo->gemm(false, false, 1.0, aotoso, Ca, 0.0);
    Cb_ao_mo->gemm(false, false, 1.0, aotoso, Cb, 0.0);

    int nirrep = Ca_ao_mo->nirrep();

    // order orbitals in terms of energy:
    int minorb;
    nmo = mos.sum();

    map = new int[nmo];
    auto *skip = new bool[nmo];
    for (int orb = 0; orb < nmo; orb++) skip[orb] = false;
    for (int orb = 0; orb < nmo; orb++) {
        int count = 0;
        double minen = 1.0e9;
        for (int h = 0; h < nirrep; h++) {
            for (int n = 0; n < wavefunction_->nmopi()[h]; n++) {
                if (skip[count]) {
                    count++;
                    continue;
                }
                if (Ea.get(h, n) <= minen) {
                    minen = Ea.get(h, n);
                    minorb = count;
                }

                count++;
            }
        }
        map[orb] = minorb;
        skip[minorb] = true;
    }

    // reorder orbitals:
    nso = wavefunction_->nso();
    eps = new double[nmo];
    sym = new int[nmo];
    occ = new int[nmo];
    Ca_pointer = new double[nmo * nso];
    for (int i = 0; i < nmo * nso; i++) Ca_pointer[i] = 0.0;

    int count = 0;
    int extra = restricted_ ? 1 : 0;
    for (int h = 0; h < nirrep; h++) {
        double **Ca_old = Ca_ao_mo->pointer(h);
        for (int n = 0; n < wavefunction_->nmopi()[h]; n++) {
            occ[count] = n < (wavefunction_->doccpi()[h] + wavefunction_->soccpi()[h]) ? 1 : 0;
            occ[count] += n < wavefunction_->doccpi()[h] ? extra : 0;
            eps[count] = Ea.get(h, n);
            sym[count] = h;
            for (int mu = 0; mu < nso; mu++) {
                Ca_pointer[mu * nmo + count] = Ca_old[mu][n];
            }
            count++;
        }
    }

    // dump to output file
    outfile->Printf("\n");
    if (restricted_)
        outfile->Printf("  ==> Molecular Orbitals <==\n");
    else
        outfile->Printf("  ==> Alpha-Spin Molecular Orbitals <==\n");
    outfile->Printf("\n");

    write_mos(mol);

    // now for beta spin
    if (!restricted_) {
        // order orbitals in terms of energy
        for (int orb = 0; orb < nmo; orb++) skip[orb] = false;
        for (int orb = 0; orb < nmo; orb++) {
            int count = 0;
            double minen = 1.0e9;
            for (int h = 0; h < nirrep; h++) {
                for (int n = 0; n < wavefunction_->nmopi()[h]; n++) {
                    if (skip[count]) {
                        count++;
                        continue;
                    }

                    if (Eb.get(h, n) <= minen) {
                        minen = Eb.get(h, n);
                        minorb = count;
                    }
                    count++;
                }
            }
            map[orb] = minorb;
            skip[minorb] = true;
        }

        // reorder orbitals:
        for (int i = 0; i < nmo * nso; i++) Ca_pointer[i] = 0.0;
        count = 0;
        for (int h = 0; h < nirrep; h++) {
            double **Ca_old = Cb_ao_mo->pointer(h);
            for (int n = 0; n < wavefunction_->nmopi()[h]; n++) {
                occ[count] = n < wavefunction_->doccpi()[h] ? 1 : 0;
                eps[count] = Eb.get(h, n);
                sym[count] = h;
                for (int mu = 0; mu < nso; mu++) {
                    Ca_pointer[mu * nmo + count] = Ca_old[mu][n];
                }
                count++;
            }
        }

        // dump to output file
        outfile->Printf("\n");
        outfile->Printf("  ==> Beta-Spin Molecular Orbitals <==\n");
        outfile->Printf("\n");
        write_mos(mol);
    }

    delete[] skip;
    delete occ;
    delete sym;
    delete eps;
    delete Ca_pointer;
}

void MOWriter::write_mos(Molecule &mol) {
    CharacterTable ct = mol.point_group()->char_table();

    BasisSet &basisset = *wavefunction_->basisset().get();

    std::vector<std::string> ao_labels;
    for (int s = 0; s < basisset.nshell(); s++) {
        GaussianShell shell = basisset.shell(s);
        int center = shell.ncenter() + 1;
        int am = shell.am();
        char amchar = shell.amchar();
        std::string basename = mol.symbol(shell.ncenter()) + std::to_string(center) + " ";
        basename += char(amchar);

        if (shell.is_pure()) {
            ao_labels.push_back(basename + "0");
            for (int i = 0; i < am; i++) {
                ao_labels.push_back(basename + "+" + std::to_string(i + 1));
                ao_labels.push_back(basename + "-" + std::to_string(i + 1));
            }
            continue;
        }

        for (int j = 0; j < am + 1; j++) {
            int lx = am - j;
            for (int lz = 0; lz < j + 1; lz++) {
                int ly = j - lz;
                std::string x = "";
                for (int count = 0; count < lx; count++) x += "x";
                std::string y = "";
                for (int count = 0; count < ly; count++) y += "y";
                std::string z = "";
                for (int count = 0; count < lz; count++) z += "z";
                std::string name = basename + x + y + z;
                ao_labels.push_back(basename + x + y + z);
            }
        }
    }

    // print mos (5 columns)
    int ncols = 5;
    int ncolsleft = nmo % ncols;
    int nrows = (nmo - ncolsleft) / ncols;

    // print the full rows:
    int count = 0;
    for (int i = 0; i < nrows; i++) {
        // print blank space
        outfile->Printf("                ");
        // print mo number
        for (int j = 0; j < ncols; j++) {
            outfile->Printf("%13d", count + j + 1);
        }
        outfile->Printf("\n");
        outfile->Printf("\n");
        // print orbitals
        for (int mu = 0; mu < nso; mu++) {
            // print ao labels
            outfile->Printf(" %-4d %-10s", mu + 1, ao_labels[mu].c_str());
            for (int j = 0; j < ncols; j++) {
                outfile->Printf("%13.7lf", Ca_pointer[mu * nmo + map[count + j]]);
            }
            outfile->Printf("\n");
        }
        outfile->Printf("\n");
        // print energy
        outfile->Printf("            Ene ");
        for (int j = 0; j < ncols; j++) {
            outfile->Printf("%13.7lf", eps[map[count + j]]);
        }
        outfile->Printf("\n");
        // print symmetry
        outfile->Printf("            Sym ");
        for (int j = 0; j < ncols; j++) {
            outfile->Printf("%13s", ct.gamma(sym[map[count + j]]).symbol());
        }
        outfile->Printf("\n");
        // print occupancy
        outfile->Printf("            Occ ");
        for (int j = 0; j < ncols; j++) {
            outfile->Printf("%13d", occ[map[count + j]]);
        }
        outfile->Printf("\n");
        outfile->Printf("\n");
        outfile->Printf("\n");
        count += ncols;
    }

    // print the partial rows:
    if (ncolsleft > 0) {
        // print blank space
        outfile->Printf("               ");
        // print mo number
        for (int j = 0; j < ncolsleft; j++) {
            outfile->Printf("%13d", count + j + 1);
        }
        outfile->Printf("\n");
        outfile->Printf("\n");
        // print orbitals
        for (int mu = 0; mu < nso; mu++) {
            // print ao labels
            outfile->Printf(" %-4d %-10s", mu + 1, ao_labels[mu].c_str());
            for (int j = 0; j < ncolsleft; j++) {
                outfile->Printf("%13.7lf", Ca_pointer[mu * nmo + map[count + j]]);
            }
            outfile->Printf("\n");
        }
        outfile->Printf("\n");
        // print energy
        outfile->Printf("            Ene ");
        for (int j = 0; j < ncolsleft; j++) {
            outfile->Printf("%13.7lf", eps[map[count + j]]);
        }
        outfile->Printf("\n");
        outfile->Printf("            Sym ");
        for (int j = 0; j < ncolsleft; j++) {
            outfile->Printf("%13s", ct.gamma(sym[map[count + j]]).symbol());
        }
        outfile->Printf("\n");
        // print occupancy
        outfile->Printf("            Occ ");
        for (int j = 0; j < ncolsleft; j++) {
            outfile->Printf("%13d", occ[map[count + j]]);
        }
        outfile->Printf("\n");
        outfile->Printf("\n");
    }
}
