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
#include <utility>
#include <algorithm>
#include "psi4/libmints/writer.h"
#include "psi4/libmints/view.h"
#include "psi4/psi4-dec.h"
#include "psi4/physconst.h"
#include "psi4/masses.h"
#include "psi4/libparallel/ParallelPrinter.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/pointgrp.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/petitelist.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/mintshelper.h"


using namespace std;
using namespace psi;
;

GradientWriter::GradientWriter(std::shared_ptr<Molecule> mol, const Matrix& grad)
    : molecule_(mol), gradient_(grad)
{
}

void GradientWriter::write(const std::string &filename)
{
   std::shared_ptr<OutFile> printer(new OutFile(filename,APPEND));
   int i;


    printer->Printf("%-59.59s %-10.10s%-9.9s\n",
            molecule_->name().c_str(),
            "(wfn)",
            "(dertype)");

    printer->Printf("%5d%20.10lf\n", molecule_->natom(), Process::environment.globals["CURRENT ENERGY"]);

    for (i=0; i<molecule_->natom(); ++i) {
        printer->Printf("%20.10lf%20.10lf%20.10lf%20.10lf\n",
                double(molecule_->Z(i)), molecule_->x(i), molecule_->y(i), molecule_->z(i));
    }

    for (i=0; i<molecule_->natom(); ++i) {
        printer->Printf("                    %20.10lf%20.10lf%20.10lf\n",
                gradient_(i, 0), gradient_(i, 1), gradient_(i, 2));
    }
}

MoldenWriter::MoldenWriter(std::shared_ptr<Wavefunction> wavefunction)
    : wavefunction_(wavefunction)
{

}
void MoldenWriter::write(const std::string &filename, std::shared_ptr<Matrix> Ca, std::shared_ptr<Matrix> Cb, std::shared_ptr<Vector> Ea, std::shared_ptr<Vector> Eb, std::shared_ptr<Vector> OccA, std::shared_ptr<Vector> OccB, bool dovirtual)
{
    std::shared_ptr<OutFile> printer(new OutFile(filename,APPEND));

    int atom;

    printer->Printf("[Molden Format]\n");

    // Get the molecule for ease
    BasisSet& basisset = *wavefunction_->basisset().get();
    Molecule& mol = *basisset.molecule().get();

    //    basisset.print_detail();

    // Print the molecule to molden
    printer->Printf("[Atoms] (AU)\n");
    for (atom=0; atom<mol.natom(); ++atom) {
        Vector3 coord = mol.xyz(atom);
        printer->Printf("%-2s  %2d  %3d   %20.12f %20.12f %20.12f\n",
                mol.symbol(atom).c_str(), atom+1, static_cast<int>(mol.Z(atom)), coord[0], coord[1], coord[2]);
    }

    // Dump the basis set using code adapted from psi2molden
    printer->Printf("[GTO]\n");

    // For each atom
    for (atom=0; atom<mol.natom(); ++atom) {
        printer->Printf("  %d 0\n", atom+1);

        // Go through all the shells on this center
        for (int shell=0; shell < basisset.nshell_on_center(atom); ++shell) {
            int overall_shell = basisset.shell_on_center(atom, shell);

            const GaussianShell& gs = basisset.shell(overall_shell);

            printer->Printf(" %c%5d  1.00\n", gs.amchar(), gs.nprimitive());

            for (int prim=0; prim<gs.nprimitive(); ++prim) {
                printer->Printf("%20.10f %20.10f\n", gs.exp(prim), gs.original_coef(prim));
            }
        }

        // An empty line separates atoms
        printer->Printf("\n");
    }

    // Convert Ca & Cb
    std::shared_ptr<PetiteList> pl(new PetiteList(wavefunction_->basisset(), wavefunction_->integral()));
    // get the "aotoso" transformation matrix, ao by so
    SharedMatrix aotoso = pl->aotoso();
    // need dimensions
    const Dimension aos = pl->AO_basisdim();
    const Dimension sos = pl->SO_basisdim();
    const Dimension nmo = Ca->colspi();

    SharedMatrix Ca_ao_mo(new Matrix("Ca AO x MO", aos, nmo));
    SharedMatrix Cb_ao_mo(new Matrix("Cb AO x MO", aos, nmo));

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
        { 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },    // p
        { 0, 3, 4, 1, 5, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0 },    // d
        { 0, 4, 5, 3, 9, 6, 1, 8, 7, 2, 0, 0, 0, 0, 0 },    // f
        { 0, 3, 4, 9, 12, 10, 5, 13, 14, 7, 1, 6, 11, 8, 2} // g
    };

    int nirrep = Ca_ao_mo->nirrep();
    Dimension countpi(nirrep);
    Dimension zeropi(nirrep);
    Dimension ncartpi(nirrep);

    for(int i = 0; i < basisset.nshell(); i++) {
        int am = basisset.shell(i).am();

        int ncart = basisset.shell(i).nfunction();
        if((am == 1 && basisset.has_puream()) || (am > 1 && am < 5 && basisset.shell(i).is_cartesian())) {
            for (int h=0; h<nirrep; ++h)
                ncartpi[h] = ncart;

            View block_a(Ca_ao_mo, ncartpi, Ca_ao_mo->colspi(), countpi, zeropi);
            View block_b(Cb_ao_mo, ncartpi, Cb_ao_mo->colspi(), countpi, zeropi);

            SharedMatrix temp_a = block_a();
            SharedMatrix temp_b = block_b();

            for( int j =0; j < ncart; j++) {
                for (int h=0; h < Ca_ao_mo->nirrep(); ++h) {
                    for (int k=0; k<Ca_ao_mo->coldim(h); ++k) {
                        // outfile->Printf( "am %d\n, from %d to %d\n", am, j, countpi[h] + molden_cartesian_order[am-1][j]);
                        Ca_ao_mo->set(h, countpi[h] + molden_cartesian_order[am-1][j], k, temp_a->get(h, j, k));
                        Cb_ao_mo->set(h, countpi[h] + molden_cartesian_order[am-1][j], k, temp_b->get(h, j, k));
                    }
                }
            }
        }

        for (int h=0; h<nirrep; ++h)
            countpi[h] += ncart;
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
    int nmoh[wavefunction_->nirrep()];
    for (int h=0; h<wavefunction_->nirrep(); ++h) {
	if (dovirtual)
	    nmoh[h] = wavefunction_->nmopi()[h];
	else
	    nmoh[h] = wavefunction_->doccpi()[h]+wavefunction_->soccpi()[h];
    }

    // do alpha's
    bool SameOcc = true;
    for (int h=0; h<wavefunction_->nirrep(); ++h) {
        for (int n=0; n<nmoh[h]; ++n) {
            mos.push_back(make_pair(Ea->get(h, n), make_pair(h, n)));
            if(fabs(OccA->get(h,n) - OccB->get(h,n)) > 1e-10)
                SameOcc = false;
        }
    }
    std::sort(mos.begin(), mos.end());

    for (int i=0; i<(int)mos.size(); ++i) {
        int h = mos[i].second.first;
        int n = mos[i].second.second;

        printer->Printf(" Sym= %s\n", ct.gamma(h).symbol());
        printer->Printf(" Ene= %20.10f\n", Ea->get(h, n));
        printer->Printf(" Spin= Alpha\n");
        if(Ca == Cb && Ea == Eb && SameOcc)
            printer->Printf(" Occup= %7.4lf\n", OccA->get(h,n)+OccB->get(h,n));
        else
            printer->Printf(" Occup= %7.4lf\n", OccA->get(h,n));
        for (int so=0; so<wavefunction_->nso(); ++so)
            printer->Printf("%3d %20.12lf\n", so+1, Ca_ao_mo->get(h, so, n));
    }

    // do beta's
    mos.clear();
    if (Ca != Cb || Ea != Eb || !SameOcc) {
        for (int h=0; h<wavefunction_->nirrep(); ++h) {
            for (int n=0; n<nmoh[h]; ++n) {
                mos.push_back(make_pair(Eb->get(h, n), make_pair(h, n)));
            }
        }
        std::sort(mos.begin(), mos.end());

        for (int i=0; i<(int)mos.size(); ++i) {
            int h = mos[i].second.first;
            int n = mos[i].second.second;

            printer->Printf(" Sym= %s\n", ct.gamma(h).symbol());
            printer->Printf(" Ene= %20.10lf\n", Eb->get(h, n));
            printer->Printf(" Spin= Beta\n");
            printer->Printf(" Occup= %7.4lf\n", OccB->get(h,n));
            for (int so=0; so<wavefunction_->nso(); ++so)
                printer->Printf("%3d %20.12lf\n", so+1, Cb_ao_mo->get(h, so, n));
        }
    }


}

FCHKWriter::FCHKWriter(std::shared_ptr<Wavefunction> wavefunction)
    : wavefunction_(wavefunction)
{
}


void FCHKWriter::write_number(const char *label, double value)
{
    fprintf(chk_, "%-43sR%27.15e\n", label, value);
}


void FCHKWriter::write_number(const char *label, int value)
{
    fprintf(chk_, "%-43sI%17d\n", label, value);
}

void FCHKWriter::write_sym_matrix(const char *label, const SharedMatrix &mat)
{
    int dim = mat->rowdim();
    fprintf(chk_, "%-43s%-3s N=%12d\n", label, "R", (dim*dim+dim)/2);

    int count = 0;
    for(int i = 0; i < dim; ++i){
        for(int j = 0; j <= i; ++j){
            fprintf(chk_, "%16.8e", mat->get(i, j));
            if(count % 5 == 4)
                fprintf(chk_, "\n");
            ++count;
        }
    }
    if(count % 5)
        fprintf(chk_, "\n");
}

void FCHKWriter::write_matrix(const char *label, const SharedVector &mat)
{
    int dim = mat->dim();
    fprintf(chk_, "%-43s%-3s N=%12d\n", label, "R", dim);

    int count = 0;
    for(int i = 0; i < dim; ++i){
        fprintf(chk_, "%16.8e", mat->get(count));
        if(count % 5 == 4)
            fprintf(chk_, "\n");
        ++count;
    }
    if(count % 5)
        fprintf(chk_, "\n");
}


void FCHKWriter::write_matrix(const char *label, const SharedMatrix &mat)
{
    int rowdim = mat->rowdim();
    int coldim = mat->coldim();
    fprintf(chk_, "%-43s%-3s N=%12d\n", label, "R", rowdim*coldim);

    int count = 0;
    for(int i = 0; i < rowdim; ++i){
        for(int j = 0; j < coldim; ++j){
            fprintf(chk_, "%16.8e", mat->get(i, j));
            if(count % 5 == 4)
                fprintf(chk_, "\n");
            ++count;
        }
    }
    if(count % 5)
        fprintf(chk_, "\n");
}

void FCHKWriter::write_matrix(const char *label, const std::vector<double> &mat)
{
    int dim = mat.size();
    fprintf(chk_, "%-43s%-3s N=%12d\n", label, "R", dim);

    int count = 0;
    for(int i = 0; i < dim; ++i){
        fprintf(chk_, "%16.8e", mat[count]);
        if(count % 5 == 4)
            fprintf(chk_, "\n");
        ++count;
    }
    if(count % 5)
        fprintf(chk_, "\n");
}

void FCHKWriter::write_matrix(const char *label, const std::vector<int> &mat)
{
    int dim = mat.size();
    fprintf(chk_, "%-43s%-3s N=%12d\n", label, "I", dim);

    int count = 0;
    for(int i = 0; i < dim; ++i){
        fprintf(chk_, "%12d", mat[count]);
        if(count % 6 == 5)
            fprintf(chk_, "\n");
        ++count;
    }
    if(count % 6)
        fprintf(chk_, "\n");
}

void FCHKWriter::write(const std::string &filename)
{
    chk_ = fopen(filename.c_str(), "w");
    std::shared_ptr<BasisSet> basis = wavefunction_->basisset();
    int maxam = basis->max_am();
    if(maxam > 4)
        throw PSIEXCEPTION("The Psi4 FCHK writer only supports up to G functions");
    std::shared_ptr<Molecule> mol = wavefunction_->molecule();
    SharedMatrix Ca_ao = wavefunction_->Ca_subset("AO");
    SharedMatrix Cb_ao = wavefunction_->Cb_subset("AO");
    SharedMatrix Da_ao = wavefunction_->Da_subset("AO");
    SharedMatrix Db_ao = wavefunction_->Db_subset("AO");
    SharedMatrix Dtot_ao(Da_ao->clone());
    SharedMatrix Dspin_ao(Da_ao->clone());
    const std::string& name = wavefunction_->name();
    const std::string& basisname = basis->name();
    Dtot_ao->add(Db_ao);
    Dspin_ao->subtract(Db_ao);
    int nbf = basis->nbf();
    int nalpha = wavefunction_->nalpha();
    int nbeta = wavefunction_->nbeta();
    int natoms = mol->natom();
    int nprimitive = basis->nprimitive();

    std::vector<double> coords;
    std::vector<double> nuc_charges;
    std::vector<int> atomic_numbers;
    std::vector<int> int_atomic_weights;
    std::vector<double> atomic_weights;
    double to_bohr = mol->units() == Molecule::Angstrom ? 1.0 / pc_bohr2angstroms : 1.0;
    for(int atom = 0; atom < natoms; ++atom){
        double Z = mol->Z(atom);
        int intZ = static_cast<int>(Z);
        double mass = an2masses[intZ];
        int intmass = static_cast<int>(mass);
        atomic_weights.push_back(mass);
        int_atomic_weights.push_back(intmass);
        nuc_charges.push_back(Z);
        atomic_numbers.push_back(intZ);
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
    // which is omitted in the CCA standard, adopted by Psi4.  We also need to
    // order basis functions to the Gaussian / GAMESS convention.  Spherical
    // harmonics are already defined appropriately, with the exception of
    // the fact that p functions are ordered 0 +1 -1, which is Z X Y, but
    // the FCHK format calls for X Y Z; this is a simple reordering operation.

    const double pureP[3][3] = {
        //           0    1    2
        // Psi4:     Z    X    Y
        // Expected: X    Y    Z
          /* 0 */ { 0.0, 1.0, 0.0 },
          /* 1 */ { 0.0, 0.0, 1.0 },
          /* 2 */ { 1.0, 0.0, 0.0 },
    };
    double pf1, pf2, pf3, pf4;
    pf1 = 1.0;            // aa
    pf2 = sqrt(1.0/3.0);  // ab
    const double cartD[6][6] = {
        //             0    1    2    3    4    5
        // Psi4:      XX   XY   XZ   YY   YZ   ZZ
        // Expected:  XX   YY   ZZ   XY   XZ   YZ
         /* 0 */   { pf1, 0.0, 0.0, 0.0, 0.0, 0.0 },
         /* 1 */   { 0.0, 0.0, 0.0, pf1, 0.0, 0.0 },
         /* 2 */   { 0.0, 0.0, 0.0, 0.0, 0.0, pf1 },
         /* 3 */   { 0.0, pf2, 0.0, 0.0, 0.0, 0.0 },
         /* 4 */   { 0.0, 0.0, pf2, 0.0, 0.0, 0.0 },
         /* 5 */   { 0.0, 0.0, 0.0, 0.0, pf2, 0.0 },
    };


    pf1 = 1.0;            // aaa
    pf2 = sqrt(1.0/5.0);  // aab
    pf3 = sqrt(1.0/15.0); // abc
    const double cartF[10][10] = {
        //            0    1    2    3    4    5    6    7    8    9
        // Psi4:     XXX  XXY  XXZ  XYY  XYZ  XZZ  YYY  YYZ  YZZ  ZZZ
        // Expected: XXX  YYY  ZZZ  XYY  XXY  XXZ  XZZ  YZZ  YYZ  XYZ
         /* 0 */   { pf1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
         /* 1 */   { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, pf1, 0.0, 0.0, 0.0 },
         /* 2 */   { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, pf1 },
         /* 3 */   { 0.0, 0.0, 0.0, pf2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
         /* 4 */   { 0.0, pf2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
         /* 5 */   { 0.0, 0.0, pf2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
         /* 6 */   { 0.0, 0.0, 0.0, 0.0, 0.0, pf2, 0.0, 0.0, 0.0, 0.0 },
         /* 7 */   { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, pf2, 0.0 },
         /* 8 */   { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, pf2, 0.0, 0.0 },
         /* 9 */   { 0.0, 0.0, 0.0, 0.0, pf3, 0.0, 0.0, 0.0, 0.0, 0.0 },
    };
    pf1 = 1.0;             // aaaa
    pf2 = sqrt(1.0/7.0);   // aaab
    pf3 = sqrt(3.0/35.0);  // aabb
    pf4 = sqrt(1.0/35.0);  // abcc
    const double cartG[15][15] = {
        //             0    1    2    3    4    5    6    7    8    9   10   11   12   13   14
        // Psi4:     XXXX XXXY XXXZ XXYY XXYZ XXZZ XYYY XYYZ XYZZ XZZZ YYYY YYYZ YYZZ YZZZ ZZZZ
        // Expected: XXXX YYYY ZZZZ XXXY XXXZ XYYY YYYZ XZZZ YZZZ XXYY XXZZ YYZZ XXYZ XYYZ XYZZ
        // This is the ordering that I would expect to generate...
         /*  0 */  { pf1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
         /*  1 */  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, pf1, 0.0, 0.0, 0.0, 0.0 },
         /*  2 */  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, pf1 },
         /*  3 */  { 0.0, pf2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
         /*  4 */  { 0.0, 0.0, pf2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
         /*  5 */  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, pf2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
         /*  6 */  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, pf2, 0.0, 0.0, 0.0 },
         /*  7 */  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, pf2, 0.0, 0.0, 0.0, 0.0, 0.0 },
         /*  8 */  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, pf2, 0.0 },
         /*  9 */  { 0.0, 0.0, 0.0, pf3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
         /* 10 */  { 0.0, 0.0, 0.0, 0.0, 0.0, pf3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
         /* 11 */  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, pf3, 0.0, 0.0 },
         /* 12 */  { 0.0, 0.0, 0.0, 0.0, pf4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
         /* 13 */  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, pf4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
         /* 14 */  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, pf4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
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

    SharedMatrix transmat(new Matrix("Reorder", nbf, nbf));
    transmat->identity();
    int offset = 0;
    for(int nshell = 0; nshell < basis->nshell(); ++nshell){
        const GaussianShell &shell = basis->shell(nshell);
        int am = shell.am();
        int nfunc = shell.nfunction();
        if(basis->has_puream()){
            // Spherical harmonics - everything is fine, apart from P orbitals
            if(am == 1){
                for(int row = 0; row < 3; ++row)
                    for(int col = 0; col < 3; ++col)
                        transmat->set(offset+row, offset+col, pureP[row][col]);
            }
        }else{
            // Cartesians - P orbitals are fine, but higher terms need reordering
            if(am == 2){
                for(int row = 0; row < 6; ++row)
                    for(int col = 0; col < 6; ++col)
                        transmat->set(offset+row, offset+col, cartD[row][col]);
            }
            if(am == 3){
                for(int row = 0; row < 10; ++row)
                    for(int col = 0; col < 10; ++col)
                        transmat->set(offset+row, offset+col, cartF[row][col]);
            }
            if(am == 4){
                for(int row = 0; row < 15; ++row)
                    for(int col = 0; col < 15; ++col)
                        transmat->set(offset+row, offset+col, cartG[row][col]);
            }
        }
        offset += nfunc;
    }

    SharedMatrix reorderedDt(Dtot_ao->clone());
    SharedMatrix reorderedDs(Dtot_ao->clone());
    reorderedDt->back_transform(Dtot_ao, transmat);
    reorderedDs->back_transform(Dspin_ao, transmat);
    SharedMatrix reorderedCa(new Matrix("Reordered Ca", Ca_ao->ncol(), Ca_ao->nrow()));
    SharedMatrix reorderedCb(new Matrix("Reordered Cb", Cb_ao->ncol(), Cb_ao->nrow()));
    reorderedCa->gemm(true, true, 1.0, Ca_ao, transmat, 0.0);
    reorderedCb->gemm(true, true, 1.0, Cb_ao, transmat, 0.0);
    for(int i = 0; i < reorderedDt->nrow(); ++i)
        for(int j = 0; j < reorderedDt->ncol(); ++j)
            if(fabs(reorderedDt->get(i,j)) < 1E-12)
                reorderedDt->set(i,j,0.0);
    for(int i = 0; i < reorderedCa->nrow(); ++i)
        for(int j = 0; j < reorderedCa->ncol(); ++j)
            if(fabs(reorderedCa->get(i,j)) < 1E-12)
                reorderedCa->set(i,j,0.0);
    for(int i = 0; i < reorderedCb->nrow(); ++i)
        for(int j = 0; j < reorderedCb->ncol(); ++j)
            if(fabs(reorderedCb->get(i,j)) < 1E-12)
                reorderedCb->set(i,j,0.0);
    std::vector<double> shell_coords;
    std::vector<double> coefficients;
    std::vector<double> exponents;
    std::vector<int> prim_per_shell;
    std::vector<int> shell_to_atom;
    std::vector<int> shell_am;
    int nshell = basis->nshell();
    for(int shell = 0; shell < nshell; ++shell){
        const GaussianShell &s = basis->shell(shell);
        int am = s.am();
        int nprim = s.nprimitive();
        prim_per_shell.push_back(nprim);
        int center = basis->shell_to_center(shell);
        shell_to_atom.push_back(center+1);
        const Vector3 &xyz = mol->xyz(center);
        shell_coords.push_back(xyz[0]);
        shell_coords.push_back(xyz[1]);
        shell_coords.push_back(xyz[2]);
        int shell_type_fac = s.am()>1 && s.is_pure() ? -1 : 1;
        shell_am.push_back(shell_type_fac * am);
        double normfac = 1.0;
        if(am > 1)
            // Undo the angular momentum normalization, applied by the ERD code
            normfac = sqrt(df[2*am]/pow(2.0, 2.0*am));
        for(int prim = 0; prim < nprim; ++prim){
            exponents.push_back(s.exp(prim));
            double normfac2 = normfac / pow(s.exp(prim), 0.5*((double)am+1.5));
            coefficients.push_back(normfac2*s.erd_coef(prim));
        }
    }


    fprintf(chk_, "Generated by Psi4\n");
    std::string jobtype = wavefunction_->gradient() ? "FORCE" :"SP";
    fprintf(chk_, "%10s%30s%30s\n", jobtype.c_str(), name.c_str(), basisname.c_str());
    write_number("Number of atoms", natoms);
    write_number("Charge", mol->molecular_charge());
    write_number("Multiplicity", mol->multiplicity());
    write_number("Number of electrons", nalpha + nbeta);
    write_number("Number of alpha electrons", nalpha);
    write_number("Number of beta electrons", nbeta);
    write_number("Number of basis functions", nbf);
    write_number("Number of independent functions", nbf);
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
    write_number("Total Energy", wavefunction_->reference_energy());
    write_matrix("Alpha Orbital Energies", wavefunction_->epsilon_a_subset("AO"));
    write_matrix("Alpha MO coefficients", reorderedCa);
    write_matrix("Beta Orbital Energies", wavefunction_->epsilon_b_subset("AO"));
    write_matrix("Beta MO coefficients", reorderedCb);
    char* label = new char[256];
    std::string type = name == "DFT" ? "SCF" : name;
    sprintf(label, "Total %s Density", type.c_str());
    write_sym_matrix(label, reorderedDt);
    sprintf(label, "Spin %s Density", type.c_str());
    write_sym_matrix(label, reorderedDs);
    delete [] label;
    SharedMatrix gradient = wavefunction_->gradient();
    if(gradient)
        write_matrix("Cartesian Gradient", gradient);

    fclose(chk_);
    chk_ = 0;
}


NBOWriter::NBOWriter(std::shared_ptr<Wavefunction> wavefunction)
    : wavefunction_(wavefunction)
{


}

void NBOWriter::write(const std::string &filename)
{
    int pure_order[][7] = {
        { 1, 0, 0, 0, 0, 0, 0},      // s
        { 101, 102, 103, 0, 0, 0, 0}, // p
        // z2  xz   yz  x2-y2 xy
        { 255, 252, 253, 254, 251, 0, 0}, // d
        //z(z2-r2), x(z2-r2), y(z2-r2) z(x2-y2), xyz, x(x2-y2), y(x2-y2)
        { 351, 352, 353, 354, 355, 356, 357 } //f
    };

    MintsHelper helper(wavefunction_->basisset(), wavefunction_->options(), 0);
    SharedMatrix sotoao = helper.petite_list()->sotoao();
    std::shared_ptr<OutFile> printer(new OutFile(filename,APPEND));


    //Get the basis set and molecule from the wavefuntion
    BasisSet& basisset = *wavefunction_->basisset().get();
    Molecule& mol = *basisset.molecule().get();

    //NBO can only handle up to f functions
    if( basisset.max_am () > 3)
    {
        throw PSIEXCEPTION("NBO cannot handle angular momentum above f functions. \n");
    }
    //print $GENNBO section of file
    //BOHR indicates atomic units for the coordinates; now ANG but not sure about keyword
    //OPEN indicates that we'll provide separate alpha and beta matrices
    printer->Printf(" $GENNBO NATOMS = %d NBAS = %d BODM ", mol.natom(), basisset.nbf());

    //To make this more user-friendly in the case of RHF wavefunctions...
    bool open_shell = (wavefunction_->nalpha() != wavefunction_->nbeta());
    if(open_shell)
        printer->Printf(" OPEN $END\n");
    else
        printer->Printf(" $END\n");

    //print NBO section of file47; user can modify this to suit their needs
    printer->Printf(" $NBO       $END\n");

    //Now print out the molecule
    printer->Printf(" $COORD\n");
    printer->Printf(" GENNBO expects one comment line here. So, here's a comment line.\n");
    for( int i =0; i< mol.natom(); i++)
    {
        //the second mol.Z() should be modified when pseudopotentials are implemented
        printer->Printf( "%2d  %2d  %20.12f %20.12f %20.12f\n",
                static_cast<int>(mol.Z(i)), static_cast<int>(mol.Z(i)),
                mol.x(i)*pc_bohr2angstroms, mol.y(i)*pc_bohr2angstroms,
                mol.z(i)*pc_bohr2angstroms);
    }
    printer->Printf( " $END\n");


    //To form the BASIS and CONTRACT sections, we need some memory
    int nshells = basisset.nshell(); //Total number of shells
    int nprim = basisset.nprimitive(); //Total number of primitives
    Vector centers(basisset.nbf());
    Vector labels(basisset.nbf());
    Vector components(nshells); //Functions per shell
    Vector angmom(nshells); //Angular momentum of shell
    Vector nprimitives(nshells); //Primitives per shell
    Vector exponents(nprim); //Exponents of primitives
    //Coefficient matrix; first row is S, second P, third D, fourth F
    Matrix coefficient(4, nprim);
    coefficient.zero();
    int fnindex = 0;
    int primindex = 0;

    //Loop over all the shells
    for( int i =0; i < nshells; i++)
    {
        const GaussianShell& gshell = basisset.shell(i);
        int nfns = gshell.nfunction(); //get number of functions in shell
        components.set(0, i, nfns);
        int angm = gshell.am(); //get angular momentum of shell
        angmom.set(0, i, angm);
        for( int j = 0; j< nfns; j++)
        {
            centers.set (0, fnindex, gshell.ncenter());
            if(gshell.is_pure()) {
                //outfile->Printf( "fnindex %d pure_order[%d][%d] %d\n", fnindex, angm, j, pure_order[angm][j]);
                labels.set (0, fnindex, pure_order[angm][j]);
            }
            else
                labels.set (0, fnindex, angm*100+j+1);
            fnindex++;
        }
        int nshellprim = gshell.nprimitive();
        nprimitives.set (0, i, nshellprim);
        for( int k =0; k < nshellprim; k++)
        {
            exponents.set(0, primindex, gshell.exp(k));
            coefficient.set (0, angm, primindex, gshell.coef(k));
            primindex++;
        }
    }

    // now, we print out the basis section
    printer->Printf(" $BASIS\n");
    // CENTER section
    for(int i=0; i<basisset.nbf(); i++)
    {
        if((i+1)%10 == 1) {
            if(i==0)
                printer->Printf("\n  CENTER =");
            else
                printer->Printf("\n          ");
        }
        printer->Printf( " %6d", (int)centers.get(0, i)+1);
    }

    //The LABEL section
    for( int i =0; i < basisset.nbf(); i++)
    {
        if((i+1)%10 == 1) {
            if(i==0)
                printer->Printf("\n   LABEL =");
            else
                printer->Printf("\n          ");
        }
        printer->Printf( " %6d", (int)labels.get(0, i));
    }
    printer->Printf( "\n $END\n");

    //The CONTRACT heading
    printer->Printf( " $CONTRACT\n");
    printer->Printf( "  NSHELL = %6d\n", nshells);
    printer->Printf( "    NEXP = %6d\n", nprim);

    // List of the number of functions per shell
    for(int i=0; i < nshells; i++)
    {
        if((i+1)%10 == 1) {
            if(i==0)
                printer->Printf("\n   NCOMP =");
            else
                printer->Printf("\n          ");
        }
        printer->Printf( " %6d", (int)components.get(0, i));
    }
    // List the number of primitives per shell
    for(int i=0; i < nshells; i++)
    {
        if((i+1)%10 == 1) {
            if(i==0)
                printer->Printf("\n   NPRIM =");
            else
                printer->Printf("\n          ");
        }
        printer->Printf( " %6d", (int)nprimitives.get(0, i));
    }
    // location of the first exponent for each shell
    int ptr = 1;
    for( int i=0; i < nshells; i++)
    {
        if((i+1)%10 == 1) {
            if(i==0)
                printer->Printf("\n    NPTR =");
            else
                printer->Printf("\n          ");
        }
        printer->Printf( " %6d", ptr);
        ptr += nprimitives.get(0, i);
    }
    // exponents
    for( int i=0; i<nprim; i++)
    {
        if((i+1)%4 == 1) {
            if(i==0)
                printer->Printf("\n     EXP =");
            else
                printer->Printf("\n          ");
        }
        printer->Printf( "%15.6E", exponents.get(0, i));
    }
    // coefficients for s functions
    for(int i=0; i<nprim; i++)
    {
        if((i+1)%4 == 1) {
            if(i==0)
                printer->Printf("\n      CS =");
            else
                printer->Printf("\n          ");
        }
        printer->Printf( "%15.6E", coefficient.get (0, 0, i));
    }
    // coefficients for p functions
    for(int i=0; i<nprim; i++)
    {
        if((i+1)%4 == 1) {
            if(i==0)
                printer->Printf("\n      CP =");
            else
                printer->Printf("\n          ");
        }
        printer->Printf( "%15.6E", coefficient.get (0, 1, i));
    }
    // coefficients for d functions
    for(int i=0; i<nprim; i++)
    {
        if((i+1)%4 == 1) {
            if(i==0)
                printer->Printf("\n      CD =");
            else
                printer->Printf("\n          ");
        }
        printer->Printf( "%15.6E", coefficient.get (0, 2, i));
    }
    // coefficients for f functions
    for(int i=0; i<nprim; i++)
    {
        if((i+1)%4 == 1) {
            if(i==0)
                printer->Printf("\n      CF =");
            else
                printer->Printf("\n          ");
        }
        printer->Printf( "%15.6E", coefficient.get (0, 3, i));
    }
    printer->Printf( "\n $END");

    //Matrix transformation information we'll need
    int nbf = basisset.nbf ();

    //Now we need the overlap matrix in the AO basis
    SharedMatrix overlap = helper.ao_overlap();
    //Print overlap matrix
    printer->Printf( "\n $OVERLAP");
    for(int i=0; i<nbf; i++)
    {
        for(int j=0; j<nbf; j++)
        {
            if(((nbf*i+j+1)%5) == 1)
                printer->Printf("\n  ");
            printer->Printf( "%15.6E", overlap->get (0, i, j));
        }
    }
    printer->Printf( "\n $END");

    //Alpha Density Matrix
    SharedMatrix soalphadens = wavefunction_->Da();
    SharedMatrix alphadens(new Matrix(nbf, nbf));
    alphadens->remove_symmetry (soalphadens, sotoao);
    //Beta density
    SharedMatrix betadens(new Matrix(nbf, nbf));
    SharedMatrix sobetadens = wavefunction_->Db();
    betadens->remove_symmetry (sobetadens, sotoao);
    //Now print the density matrix
    printer->Printf( "\n $DENSITY");
    if(wavefunction_->same_a_b_dens ())
    {
        SharedMatrix density(new Matrix(nbf, nbf));
        density->copy (alphadens);
        density->add (betadens);
        for( int i=0; i<nbf; i++)
        {
            for(int j=0; j<nbf; j++)
            {
                if(((nbf*i+j+1)%5)==1)
                    printer->Printf( "\n  ");
                printer->Printf( "%15.6E", density->get(0, i, j));
            }
        }
    }
    else
    {
        int count = 0;
        for(int i =0; i < nbf; i++)
        {
            for(int j =0; j < nbf; j++)
            {
                count++;
                if(count%5 == 1)
                    printer->Printf( "\n  ");
                printer->Printf( "%15.6E", alphadens->get (0, i, j));
            }
        }
        for(int i=0; i<nbf; i++)
        {
            for(int j=0; j<nbf; j++)
            {
                count++;
                if(count%5 ==0)
                    printer->Printf( "\n  ");
                printer->Printf( "%15.6E", betadens->get (0, i, j));
            }
        }
    }
    printer->Printf( "\n $END");


    // alpha Fock matrix
    SharedMatrix alphasofock = wavefunction_->Fa();
    SharedMatrix alphafock(new Matrix(nbf, nbf));
    alphafock->remove_symmetry (alphasofock, sotoao);
    // print the Fock matrix
    printer->Printf( "\n $FOCK");
    if(wavefunction_->same_a_b_dens ())
    {
        for(int i = 0; i < nbf; i++)
        {
            for(int j = 0; j < nbf; j++)
            {
                if(((nbf*i+j+1)%5)==1)
                    printer->Printf( "\n  ");
                printer->Printf( "%15.6E", alphafock->get (0, i, j));
            }
        }
    }

    else
    {
        // beta Fock
        SharedMatrix betafock(new Matrix(nbf, nbf));
        SharedMatrix betasofock = wavefunction_->Fb();
        betafock->remove_symmetry(betasofock, sotoao);
        int count=0;
        for(int i =0; i < nbf; i++)
        {
            for(int j =0; j < nbf; j++)
            {
                count++;
                if(count%5 == 1)
                    printer->Printf( "\n  ");
                printer->Printf( "%15.6E", alphafock->get (0, i, j));
            }
        }
        for(int i=0; i<nbf; i++)
        {
            for(int j=0; j<nbf; j++)
            {
                count++;
                if(count%5 == 1)
                    printer->Printf( "\n  ");
                printer->Printf( "%15.6E", betafock->get (0, i, j));
            }
        }
    }
    printer->Printf( "\n $END");

    //Alpha AO->MO transformation
    SharedMatrix soalphac = wavefunction_->Ca();
    const Dimension aos = helper.petite_list()->AO_basisdim();
    const Dimension nmo = wavefunction_->Ca()->colspi();
    SharedMatrix alphac(new Matrix("Ca AO x MO", aos, nmo));
    alphac->gemm(true, false, 1.00, sotoao, soalphac, 0.00);

    printer->Printf( "\n $LCAOMO");
    if(wavefunction_->same_a_b_orbs ())
    {
        for(int i = 0; i < nbf; i++)
        {
            for(int j = 0; j < nbf; j++)
            {
                if(((nbf*i+j+1)%5)==1)
                    printer->Printf( "\n  ");
                printer->Printf( "%15.6E", alphac->get(0, i, j));
            }
        }
    }

    else
    {
        //Beta AO->MO transformation
        SharedMatrix betac(new Matrix(nbf, nbf));
        SharedMatrix sobetac = wavefunction_->Cb();
        betac->gemm(true, false, 1.00, sotoao, sobetac, 0.00);

        //Print the AO->MO coefficients
        int count = 0;
        for(int i=0; i<nbf; i++)
        {
            for(int j=0; j<nbf; j++)
            {
                count++;
                if(count%5==1)
                    printer->Printf( "\n  ");
                printer->Printf( "%15.6E", alphac->get (0, i, j));
            }
        }
        for(int i =0; i < nbf; i++)
        {
            for(int j =0; j < nbf; j++)
            {
                count++;
                if(count%5 ==1)
                    printer->Printf( "\n  ");
                printer->Printf( "%15.6E ", betac->get (0,   i, j));
            }
        }
    }
    printer->Printf( "\n $END\n");

}


MOWriter::MOWriter(std::shared_ptr<Wavefunction> wavefunction)
    : wavefunction_(wavefunction), restricted_(wavefunction->same_a_b_orbs())
{
}

void MOWriter::write()
{

    // Get the molecule for ease
    BasisSet& basisset = *wavefunction_->basisset().get();
    Molecule & mol = *basisset.molecule().get();

    // Convert Ca & Cb
    // make copies
    Matrix Ca(wavefunction_->Ca());
    Matrix Cb(wavefunction_->Cb());
    Vector& Ea = *wavefunction_->epsilon_a().get();
    Vector& Eb = *wavefunction_->epsilon_b().get();

    std::shared_ptr<PetiteList> pl(new PetiteList(wavefunction_->basisset(), wavefunction_->integral()));

    // get the "aotoso" transformation matrix, ao by so
    SharedMatrix aotoso = pl->aotoso();
    // need dimensions
    const Dimension aos = pl->AO_basisdim();
    const Dimension sos = pl->SO_basisdim();
    const Dimension mos = wavefunction_->nmopi();

    SharedMatrix Ca_ao_mo(new Matrix("Ca AO x MO", aos, mos));
    SharedMatrix Cb_ao_mo(new Matrix("Cb AO x MO", aos, mos));

    // do the half transform
    Ca_ao_mo->gemm(false, false, 1.0, aotoso, Ca, 0.0);
    Cb_ao_mo->gemm(false, false, 1.0, aotoso, Cb, 0.0);

    int nirrep = Ca_ao_mo->nirrep();

    // order orbitals in terms of energy:
    int minorb;
    nmo = mos.sum();

    map  = new int[nmo];
    bool * skip = new bool[nmo];
    for (int orb = 0; orb < nmo; orb++) skip[orb] = false;
    for (int orb = 0; orb < nmo; orb++) {

        int count = 0;
        double minen = 1.0e9;
        for (int h = 0; h < nirrep; h++) {
            for (int n = 0; n<wavefunction_->nmopi()[h]; n++) {

                if ( skip[count] ) {
                    count++;
                    continue;
                }
                if ( Ea.get(h,n) <= minen ) {
                    minen = Ea.get(h,n);
                    minorb = count;
                }

                count++;

            }
        }
        map[ orb ] = minorb;
        skip[ minorb ] = true;
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
        double ** Ca_old = Ca_ao_mo->pointer(h);
        for (int n = 0; n<wavefunction_->nmopi()[h]; n++) {
            occ[ count ] = n < ( wavefunction_->doccpi()[h] + wavefunction_->soccpi()[h] ) ? 1 : 0;
            occ[ count ] += n <  wavefunction_->doccpi()[h] ? extra : 0;
            eps[ count ] = Ea.get(h,n);
            sym[ count ] = h;
            for (int mu = 0; mu < nso; mu++) {
                Ca_pointer[mu*nmo + count] = Ca_old[mu][n];
            }
            count++;

        }
    }

    // dump to output file
    outfile->Printf("\n");
    if ( restricted_ )
        outfile->Printf("  ==> Molecular Orbitals <==\n");
    else
        outfile->Printf("  ==> Alpha-Spin Molecular Orbitals <==\n");
    outfile->Printf("\n");

    write_mos(mol);

    // now for beta spin
    if ( !restricted_ ) {


        // order orbitals in terms of energy
        for (int orb = 0; orb < nmo; orb++) skip[orb] = false;
        for (int orb = 0; orb < nmo; orb++) {

            int count = 0;
            double minen = 1.0e9;
            for (int h = 0; h < nirrep; h++) {
                for (int n = 0; n<wavefunction_->nmopi()[h]; n++) {

                    if ( skip[count] ) {
                        count++;
                        continue;
                    }

                    if ( Eb.get(h,n) <= minen ) {
                        minen = Eb.get(h,n);
                        minorb = count;
                    }
                    count++;

                }
            }
            map[ orb ] = minorb;
            skip[ minorb ] = true;
        }


        // reorder orbitals:
        for (int i = 0; i < nmo * nso; i++) Ca_pointer[i] = 0.0;
        count = 0;
        for (int h = 0; h < nirrep; h++) {
            double ** Ca_old = Cb_ao_mo->pointer(h);
            for (int n = 0; n<wavefunction_->nmopi()[h]; n++) {
                occ[ count ] = n < wavefunction_->doccpi()[h] ? 1 : 0;
                eps[ count ] = Eb.get(h,n);
                sym[ count ] = h;
                for (int mu = 0; mu < nso; mu++) {
                    Ca_pointer[mu*nmo + count] = Ca_old[mu][n];
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

void MOWriter::write_mos(Molecule & mol){

    CharacterTable ct = mol.point_group()->char_table();

    BasisSet& basisset = *wavefunction_->basisset().get();
    
    std::vector<std::string> ao_labels;
    for (int s = 0; s < basisset.nshell(); s++) {
        GaussianShell shell = basisset.shell(s);
        int center = shell.ncenter()+1;
        int am = shell.am();
        char amchar = shell.amchar();
        std::string basename = mol.symbol(shell.ncenter()) + std::to_string(center) + " ";
        basename += char(amchar);
       
        if (shell.is_pure()) {
            ao_labels.push_back(basename+"0");
            for (int i = 0; i < am; i++) {
                ao_labels.push_back(basename + "+" + std::to_string(i+1));
                ao_labels.push_back(basename + "-" + std::to_string(i+1));
            }
            continue;
        }
        
        for (int j = 0; j < am+1; j++) {
            int lx = am - j;
            for (int lz = 0; lz < j+1; lz++) {
                int ly = j - lz;
                std::string x = "";
                for (int count = 0; count < lx; count++) x+="x";
                std::string y = "";
                for (int count = 0; count < ly; count++) y+="y";
                std::string z = "";
                for (int count = 0; count < lz; count++) z+="z";
                std::string name = basename + x + y + z;
                ao_labels.push_back(basename+x+y+z);
            }
        }
    }

    // print mos (5 columns)
    int ncols = 5;
    int ncolsleft = nmo % ncols;
    int nrows = (nmo - ncolsleft ) / ncols;

    // print the full rows:
    int count = 0;
    for (int i = 0; i < nrows; i++) {

        // print blank space
        outfile->Printf("                ");
        // print mo number
        for (int j = 0; j < ncols; j++){
            outfile->Printf("%13d",count+j+1);
        }
        outfile->Printf("\n");
        outfile->Printf("\n");
        // print orbitals
        for (int mu = 0; mu < nso; mu++) {
            // print ao labels
            outfile->Printf(" %-4d %-10s",mu+1,ao_labels[mu].c_str());
            for (int j = 0; j < ncols; j++){

                outfile->Printf("%13.7lf",Ca_pointer[ mu*nmo + map[count + j] ]);
            }
            outfile->Printf("\n");
        }
        outfile->Printf("\n");
        // print energy
        outfile->Printf("            Ene ");
        for (int j = 0; j < ncols; j++){
            outfile->Printf("%13.7lf",eps[ map[count + j] ]);
        }
        outfile->Printf("\n");
        // print symmetry
        outfile->Printf("            Sym ");
        for (int j = 0; j < ncols; j++){
            outfile->Printf("%13s",ct.gamma(sym[map[count+j]]).symbol());
        }
        outfile->Printf("\n");
        // print occupancy
        outfile->Printf("            Occ ");
        for (int j = 0; j < ncols; j++){
            outfile->Printf("%13d",occ[map[count+j]]);
        }
        outfile->Printf("\n");
        outfile->Printf("\n");
        outfile->Printf("\n");
        count+=ncols;
    }

    // print the partial rows:
    if ( ncolsleft > 0 ) {

        // print blank space
        outfile->Printf("               ");
        // print mo number
        for (int j = 0; j < ncolsleft; j++){
            outfile->Printf("%13d",count+j+1);
        }
        outfile->Printf("\n");
        outfile->Printf("\n");
        // print orbitals
        for (int mu = 0; mu < nso; mu++) {
            // print ao labels
            outfile->Printf(" %-4d %-10s",mu+1,ao_labels[mu].c_str());
            for (int j = 0; j < ncolsleft; j++){
                outfile->Printf("%13.7lf",Ca_pointer[ mu*nmo + map[count + j] ]);
            }
            outfile->Printf("\n");
        }
        outfile->Printf("\n");
        // print energy
        outfile->Printf("            Ene ");
        for (int j = 0; j < ncolsleft; j++){
            outfile->Printf("%13.7lf",eps[ map[count + j] ]);
        }
        outfile->Printf("\n");
        outfile->Printf("            Sym ");
        for (int j = 0; j < ncolsleft; j++){
            outfile->Printf("%13s",ct.gamma(sym[map[count+j]]).symbol());
        }
        outfile->Printf("\n");
        // print occupancy
        outfile->Printf("            Occ ");
        for (int j = 0; j < ncolsleft; j++){
            outfile->Printf("%13d",occ[map[count+j]]);
        }
        outfile->Printf("\n");
        outfile->Printf("\n");

    }
}
