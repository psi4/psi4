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

/*! \file
    \ingroup OPTKING
    \brief fd_geoms_1_0(): returns geometries necessary for finite-difference
     computation of gradients from energies; puts undisplaced geometry last in list
*/

#include "findif.h"
#include "psi4/liboptions/liboptions_python.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/cdsalclist.h"

namespace psi {
namespace findif {

std::vector<SharedMatrix> fd_geoms_1_0(std::shared_ptr<Molecule> mol, Options &options)
{

    outfile->Printf("\n-------------------------------------------------------------\n\n");

    outfile->Printf("  Using finite-differences of energies to determine gradients (fd_geoms_1_0).\n");

    int pts = options.get_int("POINTS");
    outfile->Printf("\tGenerating geometries for use with %d-point formula.\n", pts);
    if (pts != 3 && pts != 5)
        throw PsiException("FINDIF: Invalid number of points!", __FILE__, __LINE__);

    double disp_size = options.get_double("DISP_SIZE");
    outfile->Printf("\tDisplacement size will be %6.2e.\n", disp_size);

    int Natom = mol->natom();
    outfile->Printf("\tNumber of atoms is %d.\n", Natom);

    // Get SALCS from libmints
    std::shared_ptr<MatrixFactory> fact;
    py::object pyExtern = dynamic_cast<PythonDataType *>(options["EXTERN"].get())->to_python();
    bool project = !pyExtern && !options.get_bool("PERTURB_H");
    CdSalcList cdsalc(mol, fact, 0x1, project, project);

    int Nsalc = cdsalc.ncd();
    outfile->Printf("\tNumber of symmetric SALC's is %d.\n", Nsalc);

    // Determine number of geometries (1 + # of displacements)
    int Ndisp = 1;
    if (pts == 3)
        Ndisp += 2 * Nsalc;
    else if (pts == 5)
        Ndisp += 4 * Nsalc;

    outfile->Printf("\tNumber of displacements (including reference) is %d.\n", Ndisp);

    if (options.get_int("PRINT") > 1)
        for (int i = 0; i < cdsalc.ncd(); ++i)
            cdsalc[i].print();

    // Get reference geometry
    Matrix ref_geom_temp = mol->geometry();
    SharedMatrix ref_geom(ref_geom_temp.clone());

    ref_geom->set_name("Reference geometry");

    // to be returned and converted into "matrix_vector" list in python
    std::vector<SharedMatrix> disp_geoms;

    if (pts == 3) {
        for (int i = 0; i < Nsalc; ++i) {

            // - displacement
            SharedMatrix geom_m(ref_geom->clone());
            geom_m->set_name("Displacement - SALC #" + to_string(i + 1));
            displace_cart(mol, geom_m, cdsalc, i, -1, disp_size);
            disp_geoms.push_back(geom_m);

            // + displacement
            SharedMatrix geom_p(ref_geom->clone());
            geom_p->set_name("Displacement + SALC #" + to_string(i + 1));
            displace_cart(mol, geom_p, cdsalc, i, +1, disp_size);
            disp_geoms.push_back(geom_p);

        }
    } // pts 3
    else if (pts == 5) {
        for (int i = 0; i < Nsalc; ++i) {

            SharedMatrix geom_m2(ref_geom->clone());
            geom_m2->set_name("Displacement - SALC #" + to_string(i + 1) + " * 2");
            displace_cart(mol, geom_m2, cdsalc, i, -2, disp_size);
            disp_geoms.push_back(geom_m2);

            SharedMatrix geom_m1(ref_geom->clone());
            geom_m1->set_name("Displacement - SALC #" + to_string(i + 1));
            displace_cart(mol, geom_m1, cdsalc, i, -1, disp_size);
            disp_geoms.push_back(geom_m1);

            SharedMatrix geom_p1(ref_geom->clone());
            geom_p1->set_name("Displacement + SALC #" + to_string(i + 1));
            displace_cart(mol, geom_p1, cdsalc, i, +1, disp_size);
            disp_geoms.push_back(geom_p1);

            SharedMatrix geom_p2(ref_geom->clone());
            geom_p2->set_name("Displacement + SALC #" + to_string(i + 1) + " * 2");
            displace_cart(mol, geom_p2, cdsalc, i, +2, disp_size);
            disp_geoms.push_back(geom_p2);

        }
    } // pts 3

    // put reference geometry list in list
    disp_geoms.push_back(ref_geom);

    outfile->Printf("\n-------------------------------------------------------------\n");

    return disp_geoms;
}

}
}
