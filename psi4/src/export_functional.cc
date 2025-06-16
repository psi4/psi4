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

#include "psi4/pybind11.h"
#include "psi4/libfunctional/superfunctional.h"
#include "psi4/libfunctional/functional.h"
#include "psi4/libfunctional/LibXCfunctional.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/typedefs.h"
#include "psi4/libmints/numinthelper.h"
#include "psi4/libdisp/dispersion.h"
#include "psi4/libfock/v.h"
#include "psi4/libfock/points.h"
#include "psi4/libfock/cubature.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libsapt_solver/fdds_disp.h"
#include "psi4/libqt/qt.h"
#include "psi4/libpsi4util/process.h"

using namespace psi;
namespace py = pybind11;
using namespace pybind11::literals;

void export_functional(py::module &m) {
    py::class_<Functional, std::shared_ptr<Functional>>(m, "Functional", "docstring")
        .def_static("build_base", &Functional::build_base, "alias"_a, "docstring")
        .def("compute_functional", &Functional::compute_functional, "docstring")
        .def("name", &Functional::name, "docstring")
        .def("description", &Functional::description, "docstring")
        .def("citation", &Functional::citation, "docstring")
        .def("alpha", &Functional::alpha, "docstring")
        .def("omega", &Functional::omega, "docstring")
        .def("lsda_cutoff", &Functional::lsda_cutoff, "docstring")
        .def("meta_cutoff", &Functional::meta_cutoff, "docstring")
        .def("density_cutoff", &Functional::density_cutoff, "docstring")
        .def("is_gga", &Functional::is_gga, "docstring")
        .def("is_meta", &Functional::is_meta, "docstring")
        .def("is_lrc", &Functional::is_lrc, "docstring")
        .def("set_name", &Functional::set_name, "docstring")
        .def("set_description", &Functional::set_description, "docstring")
        .def("set_citation", &Functional::set_citation, "docstring")
        .def("set_gga", &Functional::set_gga, "docstring")
        .def("set_meta", &Functional::set_meta, "docstring")
        .def("set_alpha", &Functional::set_alpha, "docstring")
        .def("set_omega", &Functional::set_omega, "docstring")
        .def("set_lsda_cutoff", &Functional::set_lsda_cutoff, "docstring")
        .def("set_meta_cutoff", &Functional::set_meta_cutoff, "docstring")
        .def("set_density_cutoff", &Functional::set_density_cutoff, "docstring")
        .def("set_parameter", &Functional::set_parameter, "docstring")
        .def("print_out", &Functional::py_print, "docstring")
        .def("print_detail", &Functional::py_print_detail, "docstring");

    py::class_<BasisExtents, std::shared_ptr<BasisExtents>>(m, "BasisExtents", "docstring")
        .def(py::init<std::shared_ptr<BasisSet>, double>())
        .def("set_delta", &BasisExtents::set_delta, "docstring")
        .def("delta", &BasisExtents::delta, "docstring")
        .def("basis", &BasisExtents::basis, "docstring")
        .def("shell_extents", &BasisExtents::shell_extents, "docstring")
        .def("maxR", &BasisExtents::maxR, "docstring");

    py::class_<BlockOPoints, std::shared_ptr<BlockOPoints>>(m, "BlockOPoints", "docstring")
        .def(py::init<SharedVector, SharedVector, SharedVector, SharedVector, std::shared_ptr<BasisExtents>>())
        .def("x",
             [](BlockOPoints &grid) {
                 auto ret = std::make_shared<Vector>("X Grid points", grid.npoints());
                 C_DCOPY(grid.npoints(), grid.x(), 1, ret->pointer(), 1);
                 return ret;
             })
        .def("y",
             [](BlockOPoints &grid) {
                 auto ret = std::make_shared<Vector>("Y Grid points", grid.npoints());
                 C_DCOPY(grid.npoints(), grid.y(), 1, ret->pointer(), 1);
                 return ret;
             })
        .def("z",
             [](BlockOPoints &grid) {
                 auto ret = std::make_shared<Vector>("Z Grid points", grid.npoints());
                 C_DCOPY(grid.npoints(), grid.z(), 1, ret->pointer(), 1);
                 return ret;
             })
        .def("w",
             [](BlockOPoints &grid) {
                 auto ret = std::make_shared<Vector>("Grid Weights", grid.npoints());
                 C_DCOPY(grid.npoints(), grid.w(), 1, ret->pointer(), 1);
                 return ret;
             })
        .def("refresh", &BlockOPoints::refresh, "docstring")
        .def("npoints", &BlockOPoints::npoints, "docstring")
        .def("parent_atom", &BlockOPoints::parent_atom, "Returns the atom number this BlockOfPoints belongs to.")
        .def("print_out", &BlockOPoints::print, "out_fname"_a = "outfile", "print"_a = 2, "docstring")
        .def("shells_local_to_global", &BlockOPoints::shells_local_to_global, "docstring")
        .def("functions_local_to_global", &BlockOPoints::functions_local_to_global, "docstring");

    py::class_<SuperFunctional, std::shared_ptr<SuperFunctional>>(m, "SuperFunctional", "docstring")

        .def(py::init<>())
        .def_static("blank", &SuperFunctional::blank, "Initialize a blank SuperFunctional.")
        .def_static("XC_build", &SuperFunctional::XC_build, "name"_a, "unpolarized"_a, "tweak"_a = py::dict{}, "Builds a SuperFunctional from a XC string.")
        .def("allocate", &SuperFunctional::allocate,
             "Allocates the vectors, should be called after ansatz or npoint changes.")
        .def("compute_functional", &SuperFunctional::compute_functional,
             "vals"_a, "npoints"_a = -1, "singlet"_a = true,
             "Computes the SuperFunctional.")
        .def("x_functional", &SuperFunctional::x_functional, "Returns the desired X Functional.")
        .def("c_functional", &SuperFunctional::c_functional, "Returns the desired C Functional.")
        .def("x_functionals", &SuperFunctional::x_functionals, "Returns all X Functionals.")
        .def("c_functionals", &SuperFunctional::c_functionals, "Returns all C Functionals.")
        .def("add_x_functional", &SuperFunctional::add_x_functional, "Add a exchange Functional.")
        .def("add_c_functional", &SuperFunctional::add_c_functional, "Add a correlation Functional.")
        .def("test_functional", &SuperFunctional::test_functional, "Quick testing capabilities.")
        .def("values", &SuperFunctional::values, "Return all internal values.")
        .def("value", &SuperFunctional::value, "Returns a given internal value.")
        .def("name", &SuperFunctional::name, "The name of the SuperFunctional.")
        .def("description", &SuperFunctional::description, "The description of the SuperFunctional")
        .def("citation", &SuperFunctional::citation, "SuperFunctional citation.")
        .def("ansatz", &SuperFunctional::ansatz, "SuperFunctional rung.")
        .def("max_points", &SuperFunctional::max_points, "Maximum number of grid points per block.")
        .def("deriv", &SuperFunctional::deriv, "Maximum derivative to compute.")
        .def("x_omega", &SuperFunctional::x_omega, "Range-seperated exchange parameter.")
        .def("c_omega", &SuperFunctional::c_omega, "Range-seperated correlation parameter.")
        .def("x_alpha", &SuperFunctional::x_alpha, "Amount of exact HF exchange.")
        .def("x_beta", &SuperFunctional::x_beta, "Amount of exact HF exchange.")
        .def("c_alpha", &SuperFunctional::c_alpha, "Amount of MP2 correlation.")
        .def("c_os_alpha", &SuperFunctional::c_os_alpha, "Amount of SS MP2 correlation.")
        .def("c_ss_alpha", &SuperFunctional::c_ss_alpha, "Amount of OS MP2 correlation.")
        .def("vv10_b", &SuperFunctional::vv10_b, "The VV10 b parameter.")
        .def("vv10_c", &SuperFunctional::vv10_c, "The VV10 c parameter.")
        .def("grac_shift", &SuperFunctional::grac_shift, "Shift of the bulk potenital.")
        .def("grac_alpha", &SuperFunctional::grac_alpha, "GRAC Alpha.")
        .def("grac_beta", &SuperFunctional::grac_beta, "GRAC Beta.")
        .def("density_tolerance", &SuperFunctional::density_tolerance, "Density threshold for LibXC.")
        .def("is_gga", &SuperFunctional::is_gga, "Is this a GGA?")
        .def("is_meta", &SuperFunctional::is_meta, "Is this a MGGA?")
        .def("is_x_lrc", &SuperFunctional::is_x_lrc, "Contains range-seperated exchange?")
        .def("is_c_lrc", &SuperFunctional::is_c_lrc, "Contains range-seperated correlation?")
        .def("is_x_hybrid", &SuperFunctional::is_x_hybrid, "Requires exact exchange?")
        .def("is_c_hybrid", &SuperFunctional::is_c_hybrid, "Requires MP2 correlation?")
        .def("is_c_scs_hybrid", &SuperFunctional::is_c_scs_hybrid, "Requires SCS-MP2 correlation?")
        .def("is_libxc_func", &SuperFunctional::is_libxc_func, "A full SuperFunctional definition from LibXC.")
        .def("set_name", &SuperFunctional::set_name, "Sets the SuperFunctional name.")
        .def("set_description", &SuperFunctional::set_description, "Sets the SuperFunctional description.")
        .def("set_citation", &SuperFunctional::set_citation, "Sets the SuperFunctional citation.")
        .def("set_max_points", &SuperFunctional::set_max_points, "Sets the maximum number of points.")
        .def("set_deriv", &SuperFunctional::set_deriv, "Sets the derivative level.")
        .def("set_lock", &SuperFunctional::set_lock, "Locks the functional to prevent changes.")
        .def("set_x_omega", &SuperFunctional::set_x_omega, "Sets the range-seperation exchange parameter.")
        .def("set_c_omega", &SuperFunctional::set_c_omega, "Sets the range-seperation correlation parameter.")
        .def("set_x_alpha", &SuperFunctional::set_x_alpha, "Sets the amount of exact global HF exchange.")
        .def("set_x_beta", &SuperFunctional::set_x_beta,
             "Sets how much more long-range exchange than short-range exchange.")
        .def("set_c_alpha", &SuperFunctional::set_c_alpha, "Sets the amount of MP2 correlation.")
        .def("set_c_ss_alpha", &SuperFunctional::set_c_ss_alpha, "Sets the amount of SS MP2 correlation.")
        .def("set_c_os_alpha", &SuperFunctional::set_c_os_alpha, "Sets the amount of OS MP2 correlation.")
        .def("set_vv10_b", &SuperFunctional::set_vv10_b, "Sets the VV10 b parameter.")
        .def("set_vv10_c", &SuperFunctional::set_vv10_c, "Sets the VV10 c parameter.")
        .def("set_do_vv10", &SuperFunctional::set_do_vv10, "Sets whether to do VV10 correction.")
        .def("set_grac_shift", &SuperFunctional::set_grac_shift, "Sets the GRAC bulk shift value.")
        .def("set_grac_alpha", &SuperFunctional::set_grac_alpha, "Sets the GRAC alpha parameter.")
        .def("set_grac_beta", &SuperFunctional::set_grac_beta, "Sets the GRAC beta parameter.")
        .def("set_density_tolerance", &SuperFunctional::set_density_tolerance, "Sets the density threshold for LibXC.")
        .def("print_density_threshold", &SuperFunctional::py_print_density_threshold, "Queries the LibXCFunctionals for their density threshold values")
        .def("needs_xc", &SuperFunctional::needs_xc, "Does this functional need XC quantities.")
        .def("needs_vv10", &SuperFunctional::needs_vv10, "Does this functional need VV10 dispersion.")
        .def("needs_grac", &SuperFunctional::needs_grac, "Does this functional need GRAC.")
        .def("print_out", &SuperFunctional::py_print, "Prints out functional details.")
        .def("print_detail", &SuperFunctional::py_print_detail, "Prints all SuperFunctional information.")
        .def("xclib_description", &SuperFunctional::xclib_description, "LibXC version and citation string.")
        .def("set_xclib_description", &SuperFunctional::set_xclib_description, "Sets the LibXC version and citation string");

    typedef void (LibXCFunctional::*tweak_set1)(std::vector<double>, bool);
    typedef void (LibXCFunctional::*tweak_set2)(std::map<std::string, double>, bool);

    py::class_<LibXCFunctional, std::shared_ptr<LibXCFunctional>, Functional>(m, "LibXCFunctional", "docstring")
        .def(py::init<std::string, bool>())
        .def("get_mix_data", &LibXCFunctional::get_mix_data, "docstring")
        .def("set_tweak", tweak_set1(&LibXCFunctional::set_tweak), "tweaks"_a, "quiet"_a = false,
            "Set all tweaks on a LibXC functional through a list. Deprecated in v1.4")
        .def("set_tweak", tweak_set2(&LibXCFunctional::set_tweak), "tweaks"_a, "quiet"_a = false,
            "Set all tweaks on a LibXC functional through a dictionary of names (usually underscore prepended) and values. New in v1.4")
        .def("set_omega", &LibXCFunctional::set_omega, "docstring")
        .def("set_density_cutoff", &LibXCFunctional::set_density_cutoff, "docstring")
        .def("density_cutoff", &LibXCFunctional::density_cutoff, "docstring")
        .def("xclib_description", &LibXCFunctional::xclib_description, "query libxc for version and citation")
        .def("query_libxc", &LibXCFunctional::query_libxc, "query libxc regarding functional parameters.");

    py::class_<BasisFunctions, std::shared_ptr<BasisFunctions>>(m, "BasisFunctions", "docstring")
        .def(py::init<std::shared_ptr<BasisSet>, int, int>())
        .def("max_functions", &BasisFunctions::max_functions, "docstring")
        .def("max_points", &BasisFunctions::max_points, "docstring")
        .def("deriv", &BasisFunctions::deriv, "docstring")
        .def("set_deriv", &BasisFunctions::set_deriv, "docstring")
        .def("compute_functions", &BasisFunctions::compute_functions, "docstring")
        .def("basis_values", &BasisFunctions::basis_values, "docstring");

    typedef void (PointFunctions::*matrix_set1)(SharedMatrix);
    typedef void (PointFunctions::*matrix_set2)(SharedMatrix, SharedMatrix);

    py::class_<PointFunctions, std::shared_ptr<PointFunctions>, BasisFunctions>(m, "PointFunctions", "docstring")
        .def("print_out", &PointFunctions::print, "out_fname"_a = "outfile", "print"_a = 2, "docstring")
        .def("ansatz", &PointFunctions::ansatz, "docstring")
        .def("set_ansatz", &PointFunctions::set_ansatz, "docstring")
        .def("set_pointers", matrix_set1(&PointFunctions::set_pointers), "docstring")
        .def("set_pointers", matrix_set2(&PointFunctions::set_pointers), "docstring")
        .def("compute_points", &PointFunctions::compute_points, "block"_a, "force_compute"_a = true, "docstring")
        .def("point_values", &PointFunctions::point_values, "docstring")
        .def("orbital_values", &PointFunctions::orbital_values, "docstring");

    py::class_<MolecularGrid, std::shared_ptr<MolecularGrid>>(m, "MolecularGrid", "docstring")
        .def("print", &MolecularGrid::print, "Prints grid information.")
        .def("orientation", &MolecularGrid::orientation, "Returns the orientation of the grid.")
        .def("npoints", &MolecularGrid::npoints, "Returns the number of grid points.")
        .def("max_points", &MolecularGrid::max_points, "Returns the maximum number of points in a block.")
        .def("max_functions", &MolecularGrid::max_functions, "Returns the maximum number of functions in a block.")
        .def("collocation_size", &MolecularGrid::collocation_size, "Returns the total collocation size of all blocks.")
        .def("blocks", &MolecularGrid::blocks, "Returns a list of blocks.")
        .def("atomic_blocks", &MolecularGrid::atomic_blocks, "Returns a list of blocks.");

    py::class_<DFTGrid, std::shared_ptr<DFTGrid>, MolecularGrid>(m, "DFTGrid", "docstring")
        .def_static("build",
                    [](std::shared_ptr<Molecule> &mol, std::shared_ptr<BasisSet> &basis) {
                        return std::make_shared<DFTGrid>(mol, basis, Process::environment.options);
                    })
        .def_static("build", [](std::shared_ptr<Molecule> &mol, std::shared_ptr<BasisSet> &basis,
                                std::map<std::string, int> int_opts, std::map<std::string, std::string> string_opts) {
            return std::make_shared<DFTGrid>(mol, basis, int_opts, string_opts, Process::environment.options);
        });

    py::class_<VBase, std::shared_ptr<VBase>>(m, "VBase", "docstring")
        .def_static("build",
                    [](std::shared_ptr<BasisSet> &basis, std::shared_ptr<SuperFunctional> &func, std::string type) {
                        return VBase::build_V(basis, func, Process::environment.options, type);
                    })
        .def("initialize", &VBase::initialize, "doctsring")
        .def("finalize", &VBase::finalize, "doctsring")
        .def("basis", &VBase::basis, "Returns the internal basis set.")
        .def("functional", &VBase::functional, "Returns the interal superfunctional.")
        .def("properties", &VBase::properties, "Returns the properties computer.")
        .def("grid", &VBase::grid, "Returns the grid object.")
        .def("get_block", &VBase::get_block, "Returns the requested BlockOPoints.")
        .def("nblocks", &VBase::nblocks, "Total number of blocks.")
        .def("quadrature_values", &VBase::quadrature_values, "Returns the quadrature values.")
        .def("build_collocation_cache", &VBase::build_collocation_cache,
             "Constructs a collocation cache to prevent recomputation.")
        .def("clear_collocation_cache", &VBase::clear_collocation_cache, "Clears the collocation cache.")
        .def("set_D", &VBase::set_D, "Sets the internal density.")
        .def("Dao", &VBase::set_D, "Returns internal AO density.")
        .def("compute_V", &VBase::compute_V, "doctsring")
        .def("compute_Vx", &VBase::compute_Vx, "doctsring")
        .def("compute_gradient", &VBase::compute_gradient, "Compute the DFT nuclear gradient contribution.")
        .def("compute_hessain", &VBase::compute_hessian, "Compute the DFT nuclear Hessian contribution.")

        .def("set_print", &VBase::set_print, "Sets the print level of the object.")
        .def("set_debug", &VBase::set_debug, "Sets the debug level of the object.")

        .def("initialize", &VBase::initialize, "Initializes the V object.")
        .def("finalize", &VBase::finalize, "Finalizes the V object.")
        .def("print_header", &VBase::print_header, "Prints the objects header.");

    py::class_<RKSFunctions, std::shared_ptr<RKSFunctions>, PointFunctions>(m, "RKSFunctions", "docstring")
        .def(py::init<std::shared_ptr<BasisSet>, int, int>());

    py::class_<UKSFunctions, std::shared_ptr<UKSFunctions>, PointFunctions>(m, "UKSFunctions", "docstring")
        .def(py::init<std::shared_ptr<BasisSet>, int, int>());

    py::class_<Dispersion, std::shared_ptr<Dispersion>>(m, "Dispersion", "docstring")
        .def_static("build", &Dispersion::build, "type"_a, "s6"_a = 0.0, "alpha6"_a = 0.0, "sr6"_a = 0.0,
                    "Initialize instance capable of computing a dispersion correction of *type*")
        .def("name", &Dispersion::name, "docstring")
        .def("description", &Dispersion::description, "docstring")
        .def("citation", &Dispersion::citation, "docstring")
        .def("bibtex", &Dispersion::bibtex, "Get the BibTeX key for the literature reference.")
        .def("set_name", &Dispersion::set_name, "docstring")
        .def("set_description", &Dispersion::set_description, "docstring")
        .def("set_citation", &Dispersion::set_citation, "docstring")
        .def("set_bibtex", &Dispersion::set_bibtex, "Set the BibTeX key for the literature reference.")
        .def("print_energy", &Dispersion::print_energy, "docstring")
        .def("print_gradient", &Dispersion::print_gradient, "docstring")
        .def("print_hessian", &Dispersion::print_hessian, "docstring")
        .def("compute_energy", &Dispersion::compute_energy, "docstring")
        .def("compute_gradient", &Dispersion::compute_gradient, "docstring")
        .def("compute_hessian", &Dispersion::compute_hessian, "docstring")
        .def("d", &Dispersion::get_d, "docstring")
        .def("s6", &Dispersion::get_s6, "docstring")
        .def("sr6", &Dispersion::get_sr6, "docstring")
        .def("s8", &Dispersion::get_s8, "docstring")
        .def("a1", &Dispersion::get_a1, "docstring")
        .def("a2", &Dispersion::get_a2, "docstring")
        .def("print_out", &Dispersion::py_print, "docstring");

    py::class_<sapt::FDDS_Dispersion, std::shared_ptr<sapt::FDDS_Dispersion>>(m, "FDDS_Dispersion", "docstring")
        .def(py::init<std::shared_ptr<BasisSet>, std::shared_ptr<BasisSet>, std::map<std::string, SharedMatrix>,
                      std::map<std::string, SharedVector>, bool>())
        .def("metric", &sapt::FDDS_Dispersion::metric, "Obtains the FDDS metric.")
        .def("metric_inv", &sapt::FDDS_Dispersion::metric_inv, "Obtains the FDDS metric_inv.")
        .def("aux_overlap", &sapt::FDDS_Dispersion::aux_overlap, "Obtains the FDDS aux_overlap.")
        .def("project_densities", &sapt::FDDS_Dispersion::project_densities,
             "Projects a density from the primary AO to auxiliary AO space.")
        .def("form_unc_amplitude", &sapt::FDDS_Dispersion::form_unc_amplitude,
             "Forms the uncoupled amplitudes for either monomer.")
        .def("get_tensor_pqQ", &sapt::FDDS_Dispersion::get_tensor_pqQ,
             "Debug only: fetches 3-index intermediate from disk and return as matrix.")
        .def("print_tensor_pqQ", &sapt::FDDS_Dispersion::print_tensor_pqQ,
             "Debug only: prints formatted 3-index intermediate to file.")
        .def("form_aux_matrices", &sapt::FDDS_Dispersion::form_aux_matrices,
             "Forms the uncoupled amplitudes and other matrices for either monomer.")
        .def("R_A", &sapt::FDDS_Dispersion::R_A, "Obtains (R^t)^-1 for monomer A.")
        .def("R_B", &sapt::FDDS_Dispersion::R_B, "Obtains (R^t)^-1 for monomer B.");

     py::class_<NumIntHelper, std::shared_ptr<NumIntHelper>>(m, "NumIntHelper",
                                                             "Computes numerical integrals using a DFT grid.")
         .def(py::init<std::shared_ptr<DFTGrid>>())
         .def("numint_grid", &NumIntHelper::numint_grid)
         .def("density_integral", &NumIntHelper::density_integral,
              "Compute an integral \\int \\rho(r) f(r) where f is a vector-valued function. f is represented for each "
              "block of points of the integration grid as a matrix (n_data, n_points). Return has shape (n_data)",
              "grid_data"_a, "D"_a)
         .def("dd_density_integral", &NumIntHelper::dd_density_integral,
              "Compute an integral \\int \\rho(r) f(r) where f is a vector-valued function. f is represented for each "
              "block of points of the integration grid as a matrix (n_data, n_points). Return has shape (n_atoms, "
              "n_data)", "grid_data"_a, "D"_a)
         .def("potential_integral", &NumIntHelper::potential_integral,
              "Compute an integral \\int \\chi_\\mu(r) \\chi_\\nu(r) f(r) where f is a scalar function represented for "
              "each block of points of the integration grid as a vector of n_points.");
}
