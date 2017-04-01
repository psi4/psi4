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
    \brief fd_freq_1(): compute frequencies from gradients
*/

#include "findif.h"
#include "psi4/libmints/writer_file_prefix.h"
#include "psi4/libparallel/ParallelPrinter.h"
#include "psi4/liboptions/liboptions_python.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/cdsalclist.h"
#include "psi4/libmints/pointgrp.h"
#include "psi4/libmints/petitelist.h"
#include "psi4/physconst.h"

#include "psi4/pybind11.h"

namespace psi {
namespace findif {

SharedMatrix fd_freq_1(std::shared_ptr <Molecule> mol, Options &options,
                       const py::list &grad_list, int freq_irrep_only)
{
    int pts = options.get_int("POINTS");
    double disp_size = options.get_double("DISP_SIZE");
    int print_lvl = options.get_int("PRINT");

    int Natom = mol->natom();
    std::shared_ptr <MatrixFactory> fact;
    py::object pyExtern = dynamic_cast<PythonDataType *>(options["EXTERN"].get())->to_python();
    bool project = !pyExtern && !options.get_bool("PERTURB_H");

    CdSalcList salc_list(mol, fact, 0xFF, project, project);
    int Nirrep = salc_list.nirrep();

    // *** Build vectors that list indices of salcs for each irrep
    std::vector <std::vector<int>> salcs_pi;
    for (int h = 0; h < Nirrep; ++h)
        salcs_pi.push_back(std::vector<int>());
    for (int i = 0; i < salc_list.ncd(); ++i)
        salcs_pi[salc_list[i].irrep()].push_back(i);

    // Now remove irreps that are not requested
    if (freq_irrep_only != -1) {
        for (int h = 0; h < Nirrep; ++h)
            if (h != freq_irrep_only)
                salcs_pi[h].clear();
    }

    // Determine total num of salcs and where each irrep starts
    int Nsalc_all = salcs_pi[0].size();
    int salc_irr_start[8];
    salc_irr_start[0] = 0;
    for (int h = 1; h < Nirrep; ++h) {
        Nsalc_all += salcs_pi[h].size();
        salc_irr_start[h] = salc_irr_start[h - 1] + salcs_pi[h - 1].size();
    }

    // ** Count displacements
    // symmetric displacements:
    std::vector<int> Ndisp_pi(Nirrep);
    if (pts == 3)
        Ndisp_pi[0] = 2 * salcs_pi[0].size();
    else if (pts == 5)
        Ndisp_pi[0] = 4 * salcs_pi[0].size();

    // asymmetric displacements:
    for (int h = 1; h < Nirrep; ++h) {
        if (pts == 3)
            Ndisp_pi[h] = salcs_pi[h].size();
        else if (pts == 5)
            Ndisp_pi[h] = 2 * salcs_pi[h].size();
    }
    int Ndisp_all = 0;
    for (int h = 0; h < Nirrep; ++h)
        Ndisp_all += Ndisp_pi[h];

    outfile->Printf("\n-------------------------------------------------------------\n\n");

    outfile->Printf("  Computing second-derivative from gradients using projected, \n");
    outfile->Printf("  symmetry-adapted, cartesian coordinates (fd_freq_1).\n\n");

    outfile->Printf("  %d gradients passed in, including the reference geometry.\n", (int) len(grad_list));

    // We are passing in the reference geometry at the moment, though we are not using
    // its gradient.  Could be removed later.

    if ((int) len(grad_list) != Ndisp_all + 1) { // last gradient is the reference, non-displaced one
        outfile->Printf("gradients.size() is %d\n", (int) len(grad_list));
        outfile->Printf("Ndisp_all is %d\n", Ndisp_all);
        throw PsiException("FINDIF: Incorrect number of gradients passed in!", __FILE__, __LINE__);
    }

    // *** Generate complete list of gradients from unique ones.
    outfile->Printf("  Generating complete list of displacements from unique ones.\n\n");

    std::shared_ptr <PointGroup> pg = mol->point_group();
    CharacterTable ct = mol->point_group()->char_table();
    int order = ct.order();

    // atom_map, how atoms are mapped to other atoms by operations
    int **atom_map = compute_atom_map(mol);
    if (print_lvl >= 3) {
        outfile->Printf("\tThe atom map:\n");
        for (int i = 0; i < Natom; ++i) {
            outfile->Printf("\t %d : ", i + 1);
            for (int j = 0; j < order; ++j)
                outfile->Printf("%4d", atom_map[i][j] + 1);
            outfile->Printf("\n");
        }
        outfile->Printf("\n");
    }

    // Extract the symmetric gradients.
    std::vector <SharedMatrix> gradients;
    for (int i = 0; i < Ndisp_pi[0]; ++i)
        gradients.push_back(grad_list[i].cast<SharedMatrix>());

    if (print_lvl >= 3) {
        outfile->Printf("\tSymmetric gradients\n");
        for (int i = 0; i < gradients.size(); ++i)
            gradients[i]->print();
    }

    // Extract the asymmetric gradients, one at a time and determine the gradient of the
    // non-computed displacements.

    int disp_cnt = Ndisp_pi[0]; // step through original list of gradients for non-symmetric ones

    for (int h = 1; h < Nirrep; ++h) { // loop over asymmetric irreps

        if (Ndisp_pi[h] == 0) continue;

        IrreducibleRepresentation gamma = ct.gamma(h);

        if (print_lvl >= 3) {
            outfile->Printf("Characters for irrep %d\n", h);
            for (int i = 0; i < order; ++i)
                outfile->Printf(" %5.1lf", gamma.character(i));
            outfile->Printf("\n");
        }

        // Find operation that takes + to - displacement.
        int op_disp;
        for (op_disp = 0; op_disp < order; ++op_disp)
            if (gamma.character(op_disp) == -1)
                break;
        outfile->Printf("\tOperation %d takes plus displacements of irrep %s to minus ones.\n",
                        op_disp + 1, gamma.symbol());

        // Get 3x3 matrix representation of operation.
        SymmetryOperation so = ct.symm_operation(op_disp);

        // Loop over coordinates of that irrep.
        for (int coord = 0; coord < salcs_pi[h].size(); ++coord) {

            // Read the - displacement and generate the +
            gradients.push_back(grad_list[disp_cnt].cast<SharedMatrix>());

            SharedMatrix new_grad(new Matrix(Natom, 3));

            for (int atom = 0; atom < Natom; ++atom) {
                int atom2 = atom_map[atom][op_disp]; // how this atom transforms under this op.

                for (int xyz2 = 0; xyz2 < 3; ++xyz2) { // target xyz
                    double tval = 0.0;
                    for (int xyz = 0; xyz < 3; ++xyz)   // original xyz
                        tval += so(xyz2, xyz) * gradients.back()->get(atom, xyz);
                    new_grad->set(atom2, xyz2, tval);
                }
            }
            ++disp_cnt;

            // if pts == 5, then read -1 displacement, generate +1 and insert it so order is -2,-1,+1,+2
            if (pts == 5) {
                gradients.push_back(grad_list[disp_cnt].cast<SharedMatrix>());

                SharedMatrix new_grad2(new Matrix(Natom, 3));

                for (int atom = 0; atom < Natom; ++atom) {
                    int atom2 = atom_map[atom][op_disp]; // how this atom transforms under this op.

                    for (int xyz2 = 0; xyz2 < 3; ++xyz2) { // target xyz
                        double tval = 0.0;
                        for (int xyz = 0; xyz < 3; ++xyz)   // original xyz
                            tval += so(xyz2, xyz) * gradients.back()->get(atom, xyz);
                        new_grad2->set(atom2, xyz2, tval);
                    }
                }
                ++disp_cnt;

                gradients.push_back(new_grad2); // put +1 gradient in list
            } // end extra gradient for 5-pt. formula

            gradients.push_back(new_grad); // put +1 (3pt.) or +2 (5pt.) gradient in list
        } // end coord
    }

    delete_atom_map(atom_map, mol);

    // Fix number of displacements for full list.
    for (int h = 0; h < Nirrep; ++h) {
        if (pts == 3)
            Ndisp_pi[h] = 2 * salcs_pi[h].size();
        else if (pts == 5)
            Ndisp_pi[h] = 4 * salcs_pi[h].size();
    }

    int disp_irr_start[8];
    disp_irr_start[0] = 0;
    Ndisp_all = Ndisp_pi[0];
    for (int h = 1; h < Nirrep; ++h) {
        Ndisp_all += Ndisp_pi[h];
        disp_irr_start[h] = disp_irr_start[h - 1] + Ndisp_pi[h - 1];
    }

    // Mass-weight all the gradients g_xm = 1/sqrt(m) g_x
    for (int i = 0; i < Ndisp_all; ++i) {
        double **disp = gradients[i]->pointer();
        for (int a = 0; a < Natom; ++a)
            for (int xyz = 0; xyz < 3; ++xyz)
                disp[a][xyz] /= sqrt(mol->mass(a));
    }

    if (print_lvl >= 3) {
        outfile->Printf("\tAll mass-weighted gradients\n");
        for (int i = 0; i < gradients.size(); ++i)
            gradients[i]->print();
    }

    char **irrep_lbls = mol->irrep_labels();
    double **H_irr[8];

    std::vector < VIBRATION * > modes;

    for (int h = 0; h < Nirrep; ++h) {

        if (salcs_pi[h].size() == 0) continue;

        // To store gradients in SALC displacement coordinates.
        double **grads_adapted = block_matrix(Ndisp_pi[h], salcs_pi[h].size());

        // Build B matrix / sqrt(masses).
        SharedMatrix B_irr_shared = salc_list.matrix_irrep(h);
        double **B_irr = B_irr_shared->pointer();

        // Compute forces in internal coordinates, g_q = G_inv B u g_x
        // In this case, B = c * masses^(1/2).  =>  G=I.
        // Thus, g_q = c * g_x / sqrt(masses) or B g_x = g_q.
        for (int disp = 0; disp < Ndisp_pi[h]; ++disp)
            for (int salc = 0; salc < salcs_pi[h].size(); ++salc)
                for (int a = 0; a < Natom; ++a)
                    for (int xyz = 0; xyz < 3; ++xyz)
                        grads_adapted[disp][salc] += B_irr[salc][3 * a + xyz] *
                                                     gradients[disp_irr_start[h] + disp]->get(a, xyz);

        if (print_lvl >= 3) {
            outfile->Printf("Gradients in B-matrix coordinates\n");
            for (int disp = 0; disp < Ndisp_pi[h]; ++disp) {
                outfile->Printf(" disp %d: ", disp);
                for (int salc = 0; salc < salcs_pi[h].size(); ++salc)
                    outfile->Printf("%15.10lf", grads_adapted[disp][salc]);
                outfile->Printf("\n");
            }
        }

        /* Test forces by recomputed cartesian, mass-weighted gradient: // B^t f_q = f_x
        outfile->Printf("Test gradients - recomputed\n");
        for (int disp=0; disp<Ndisp_pi[h]; ++disp) {
          outfile->Printf( "g_x %d : \n", disp);
          for (int a=0; a<Natom; ++a) {
            for (int xyz=0; xyz<3; ++xyz) {
              double tval = 0;
              for (int salc=0; salc<salcs_pi[h].size(); ++salc)
                tval += B_irr[salc][3*a+xyz] * grads_adapted[disp][salc];
              outfile->Printf("%15.10lf", tval);
            }
            outfile->Printf("\n");
          }
        }*/

        //** Construct force constant matrix from finite differences of forces
        H_irr[h] = init_matrix(salcs_pi[h].size(), salcs_pi[h].size());

        if (pts == 3) { // Hij = fj(i+1) - fj(i-1) / (2h)

            for (int i = 0; i < salcs_pi[h].size(); ++i)
                for (int j = 0; j < salcs_pi[h].size(); ++j)
                    H_irr[h][i][j] = (grads_adapted[2 * i + 1][j] - grads_adapted[2 * i][j]) / (2.0 * disp_size);

        } else if (pts == 5) { // fj(i-2) - 8fj(i-1) + 8fj(i+1) - fj(i+2) / (12h)

            for (int i = 0; i < salcs_pi[h].size(); ++i)
                for (int j = 0; j < salcs_pi[h].size(); ++j)
                    H_irr[h][i][j] = (1.0 * grads_adapted[4 * i][j] - 8.0 * grads_adapted[4 * i + 1][j]
                                      + 8.0 * grads_adapted[4 * i + 2][j] - 1.0 * grads_adapted[4 * i + 3][j])
                                     / (12.0 * disp_size);

        }

        if (print_lvl >= 3) {
            outfile->Printf("\n\tForce Constants for irrep %s in mass-weighted, ", irrep_lbls[h]);
            outfile->Printf("symmetry-adapted cartesian coordinates.\n");
            mat_print(H_irr[h], salcs_pi[h].size(), salcs_pi[h].size(), "outfile");
        }

        // diagonalize force constant matrix
        int dim = salcs_pi[h].size();
        double *evals = init_array(dim);
        double **evects = block_matrix(dim, dim);

        sq_rsp(dim, dim, H_irr[h], evals, 3, evects, 1e-14);

        // Bu^1/2 * evects -> normal mode
        for (int i = 0; i < dim; ++i)
            for (int a = 0; a < Natom; ++a)
                for (int xyz = 0; xyz < 3; ++xyz)
                    B_irr[i][3 * a + xyz] /= sqrt(mol->mass(a));

        double **normal_irr = block_matrix(3 * Natom, dim);
        C_DGEMM('t', 'n', 3 * Natom, dim, dim, 1.0, B_irr[0], 3 * Natom, evects[0],
                dim, 0, normal_irr[0], dim);

        if (print_lvl >= 2) {
            outfile->Printf("\n\tNormal coordinates (non-mass-weighted) for irrep %s:\n", irrep_lbls[h]);
            eivout(normal_irr, evals, 3 * Natom, dim, "outfile");
        }

        for (int i = 0; i < salcs_pi[h].size(); ++i) {
            double *v = init_array(3 * Natom);
            for (int x = 0; x < 3 * Natom; ++x)
                v[x] = normal_irr[x][i];
            VIBRATION *vib = new VIBRATION(h, evals[i], v);
            modes.push_back(vib);
        }

        free(evals);
        free_block(evects);
        free_block(normal_irr);
    }

    // This print function also saves frequencies in wavefunction.
    print_vibrations(mol, modes);

    // Optionally, save normal modes to file.
    if (options.get_bool("NORMAL_MODES_WRITE")) {
        save_normal_modes(mol, modes);
    }

    for (int i = 0; i < modes.size(); ++i)
        delete modes[i];
    modes.clear();

    // Build complete hessian for transformation to cartesians
    double **H = block_matrix(Nsalc_all, Nsalc_all);

    for (int h = 0; h < Nirrep; ++h)
        for (int i = 0; i < salcs_pi[h].size(); ++i) {
            int start = salc_irr_start[h];
            for (int j = 0; j <= i; ++j)
                H[start + i][start + j] = H[start + j][start + i] = H_irr[h][i][j];
        }

    for (int h = 0; h < Nirrep; ++h)
        if (salcs_pi[h].size()) free_block(H_irr[h]);

    // Transform Hessian into cartesian coordinates
    if (print_lvl >= 3) {
        outfile->Printf("\n\tFull force constant matrix in mass-weighted SALCS.\n");
        mat_print(H, Nsalc_all, Nsalc_all, "outfile");
    }

    // Build Bu^-1/2 matrix for the whole Hessian
    SharedMatrix B_shared = salc_list.matrix();
    double **B = B_shared->pointer();

    // double **Hx = block_matrix(3*Natom, 3*Natom);
    SharedMatrix mat_Hx = SharedMatrix(new Matrix("Hessian", 3 * Natom, 3 * Natom));
    double **Hx = mat_Hx->pointer();

    // Hx = Bt H B
    for (int i = 0; i < Nsalc_all; ++i)
        for (int j = 0; j < Nsalc_all; ++j)
            for (int x1 = 0; x1 < 3 * Natom; ++x1)
                for (int x2 = 0; x2 <= x1; ++x2)
                    Hx[x1][x2] += B[i][x1] * H[i][j] * B[j][x2];

    for (int x1 = 0; x1 < 3 * Natom; ++x1)
        for (int x2 = 0; x2 < x1; ++x2)
            Hx[x2][x1] = Hx[x1][x2];

    free_block(H);

    if (print_lvl >= 3) {
        outfile->Printf("\n\tForce Constants in mass-weighted cartesian coordinates.\n");
        mat_print(Hx, 3 * Natom, 3 * Natom, "outfile");
    }

    // Un-mass-weight Hessian
    for (int x1 = 0; x1 < 3 * Natom; ++x1)
        for (int x2 = 0; x2 < 3 * Natom; ++x2)
            Hx[x1][x2] *= sqrt(mol->mass(x1 / 3)) * sqrt(mol->mass(x2 / 3));

    if (print_lvl >= 3) {
        outfile->Printf("\n\tForce Constants in cartesian coordinates.\n");
        mat_print(Hx, 3 * Natom, 3 * Natom, "outfile");
    }

    // Print a hessian file
    if (options.get_bool("HESSIAN_WRITE")) {
        std::string hess_fname = get_writer_file_prefix(mol->name()) + ".hess";
        std::shared_ptr <OutFile> printer(new OutFile(hess_fname, TRUNCATE));
        //FILE *of_Hx = fopen(hess_fname.c_str(),"w");
        printer->Printf("%5d", Natom);
        printer->Printf("%5d\n", 6 * Natom);

        int cnt = -1;
        for (int i = 0; i < 3 * Natom; ++i) {
            for (int j = 0; j < 3 * Natom; ++j) {
                printer->Printf("%20.10lf", Hx[i][j]);
                if (++cnt == 2) {
                    printer->Printf("\n");
                    cnt = -1;
                }
            }
        }
    }
//  free_block(Hx);

    outfile->Printf("\n-------------------------------------------------------------\n");

    return mat_Hx;
}

}
}
