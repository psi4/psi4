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
    \brief fd_1_0(): compute gradient using energies and finite-differences
*/

#include "findif.h"
#include "psi4/libmints/writer_file_prefix.h"
#include "psi4/liboptions/liboptions_python.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/writer.h"
#include "psi4/libmints/factory.h"
#include "psi4/libmints/cdsalclist.h"

#include "psi4/pybind11.h"

namespace psi {
namespace findif {

SharedMatrix fd_1_0(std::shared_ptr <Molecule> mol, Options &options, const py::list &python_energies)
{
    int pts = options.get_int("POINTS");
    double disp_size = options.get_double("DISP_SIZE");

    int Natom = mol->natom();
    std::shared_ptr <MatrixFactory> fact;

    py::object pyExtern = dynamic_cast<PythonDataType *>(options["EXTERN"].get())->to_python();
    bool project = !pyExtern && !options.get_bool("PERTURB_H");
    CdSalcList cdsalc(mol, fact, 0x1, project, project);
    int Nsalc = cdsalc.ncd();

    // Compute number of displacements - check with number of energies passed in
    // Determine number of geometries (1 + # of displacements)
    int Ndisp = 1;
    if (pts == 3)
        Ndisp += 2 * Nsalc;
    else if (pts == 5)
        Ndisp += 4 * Nsalc;
    else
        throw PSIEXCEPTION("fd_1_0: Unable to handle requested point formula. 3 or 5-point formula are supported.");

    if (len(python_energies) != Ndisp)
        throw PsiException("FINDIF: Incorrect number of energies passed in!", __FILE__, __LINE__);

    double *E = init_array(Ndisp);
    for (int i = 0; i < Ndisp; ++i)
        E[i] = python_energies[i].cast<double>();

    // Compute gradient in mass-weighted symmetry-adapted cartesians in ATOMIC units
    double *g_q = init_array(Nsalc);
    if (pts == 3) {
        for (int i = 0; i < Nsalc; ++i)
            g_q[i] = (E[2 * i + 1] - E[2 * i]) / (2.0 * disp_size);
    } else if (pts == 5) {
        for (int i = 0; i < Nsalc; ++i)
            g_q[i] = (E[4 * i] - 8.0 * E[4 * i + 1] + 8.0 * E[4 * i + 2] - E[4 * i + 3]) / (12.0 * disp_size);
    }

    outfile->Printf("\n-------------------------------------------------------------\n\n");
    outfile->Printf("  Computing gradient from energies (fd_1_0).\n");

    // Print out energies and gradients
    double energy_ref = E[Ndisp - 1];
    outfile->Printf("\tUsing %d-point formula.\n", pts);
    outfile->Printf("\tEnergy without displacement: %15.10lf\n", energy_ref);
    outfile->Printf("\tCheck energies below for precision!\n");
    outfile->Printf("\tForces are for mass-weighted, symmetry-adapted cartesians (in au).\n");

    int cnt;
    if (pts == 3) {
        cnt = -2;
        outfile->Printf("\n\t Coord      Energy(-)        Energy(+)        Force\n");
        for (int i = 0; i < Nsalc; ++i) {
            cnt += 2;
            outfile->Printf("\t%5d %17.10lf%17.10lf%17.10lf\n", i, E[cnt], E[cnt + 1], g_q[i]);
        }
        outfile->Printf("\n");
    } else if (pts == 5) {
        cnt = -4;
        outfile->Printf(
                "\n\t Coord      Energy(-2)        Energy(-1)        Energy(+1)        Energy(+2)            Force\n");
        for (int i = 0; i < Nsalc; ++i) {
            cnt += 4;
            outfile->Printf("\t%5d %17.10lf %17.10lf %17.10lf %17.10lf %17.10lf\n",
                            i, E[cnt], E[cnt + 1], E[cnt + 2], E[cnt + 3], g_q[i]);
        }
        outfile->Printf("\n");
    }

    // Build B matrix of salc coefficients
    SharedMatrix Bmat = cdsalc.matrix();
    double **B = Bmat->pointer();

    // compute gradient in mass-weighted (non-SALC) cartesians
    double *g_cart = init_array(3 * Natom);

    // B^t g_q^t = g_x^t -> g_q B = g_x
    C_DGEMM('n', 'n', 1, 3 * Natom, Nsalc, 1.0, g_q, Nsalc, B[0], 3 * Natom, 0, g_cart, 3 * Natom);

    free(g_q);

    // The undisplaced geometry should be in the global molecule, and the undisplaced
    // energy in globals["CURRENT ENERGY"], since we did that one last.  Clever, huh.

    // Un-massweight the gradient and save it
    Matrix gradient_matrix("F-D gradient", Natom, 3);

    for (int a = 0; a < Natom; ++a)
        for (int xyz = 0; xyz < 3; ++xyz)
            gradient_matrix.set(a, xyz, g_cart[3 * a + xyz] * sqrt(mol->mass(a)));

    free(g_cart);

    // Print a gradient file
    if (options.get_bool("GRADIENT_WRITE")) {
        GradientWriter grad(mol, gradient_matrix);
        std::string gradfile = get_writer_file_prefix(mol->name()) + ".grad";
        grad.write(gradfile);
        outfile->Printf("\tGradient written.\n");
    }

    SharedMatrix sgradient(gradient_matrix.clone());
    outfile->Printf("\n-------------------------------------------------------------\n");

    return sgradient;
}

}
}
