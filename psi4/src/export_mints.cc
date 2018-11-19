/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2018 The Psi4 Developers.
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

#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl_bind.h>
#include <pybind11/operators.h>

#include "psi4/libmints/basisset.h"
#include "psi4/libmints/deriv.h"
#include "psi4/libmints/twobody.h"
#include "psi4/libmints/integralparameters.h"
#include "psi4/libmints/orbitalspace.h"
#include "psi4/libmints/local.h"
#include "psi4/libmints/vector3.h"
#include "psi4/libmints/pointgrp.h"
#include "psi4/libmints/extern.h"
#include "psi4/libmints/sobasis.h"
#include "psi4/libmints/petitelist.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/local.h"
#include "psi4/libmints/sointegral_onebody.h"
#include "psi4/libmints/sointegral_twobody.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/multipolesymmetry.h"
#include "psi4/libmints/eri.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/3coverlap.h"
#include "psi4/libmints/pseudospectral.h"
#include "psi4/libmints/oeprop.h"
#include "psi4/libmints/nabla.h"
#include "psi4/libmints/electrostatic.h"
#include "psi4/libmints/potential.h"
#include "psi4/libmints/kinetic.h"
#include "psi4/libmints/factory.h"
#include "psi4/libmints/writer.h"
#include "psi4/libmints/corrtab.h"
#include "psi4/libmints/electricfield.h"
#include "psi4/libmints/tracelessquadrupole.h"
#include "psi4/libmints/angularmomentum.h"
#include "psi4/libmints/multipoles.h"
#include "psi4/libmints/quadrupole.h"
#include "psi4/libmints/dipole.h"
#include "psi4/libmints/overlap.h"

#include <string>

using namespace psi;
namespace py = pybind11;

/** Returns a new basis set object
 * Constructs a basis set from the parsed information
 *
 * @param mol           Psi4 molecule
 * @param py::dict      Python dictionary containing the basis information
 * @param forced_puream Force puream or not
 **/
std::shared_ptr<BasisSet> construct_basisset_from_pydict(const std::shared_ptr<Molecule>& mol, py::dict& pybs,
                                                         const int forced_puream) {
    std::string key = pybs["key"].cast<std::string>();
    std::string name = pybs["name"].cast<std::string>();
    std::string label = pybs["blend"].cast<std::string>();

    // Handle mixed puream signals and seed parser with the resolution
    int native_puream = pybs["puream"].cast<int>();
    int user_puream = (Process::environment.options.get_global("PUREAM").has_changed())
                          ? ((Process::environment.options.get_global("PUREAM").to_integer()) ? Pure : Cartesian)
                          : -1;
    GaussianType shelltype;
    if (user_puream == -1)
        shelltype = static_cast<GaussianType>(forced_puream == -1 ? native_puream : forced_puream);
    else
        shelltype = static_cast<GaussianType>(user_puream);

    mol->set_basis_all_atoms(name, key);

    // Map of GaussianShells: basis_atom_shell[basisname][atomlabel] = gaussian_shells
    typedef std::map<std::string, std::map<std::string, std::vector<ShellInfo>>> map_ssv;
    map_ssv basis_atom_shell;
    // basisname is uniform; fill map with key/value (gbs entry) pairs of elements from pybs['shell_map']
    py::list basisinfo = pybs["shell_map"].cast<py::list>();
    if (len(basisinfo) == 0) throw PSIEXCEPTION("Empty information being used to construct BasisSet.");
    for (int atom = 0; atom < py::len(basisinfo); ++atom) {
        std::vector<ShellInfo> vec_shellinfo;
        py::list atominfo = basisinfo[atom].cast<py::list>();
        std::string atomlabel = atominfo[0].cast<std::string>();
        std::string hash = atominfo[1].cast<std::string>();
        for (int atomshells = 2; atomshells < py::len(atominfo); ++atomshells) {
            // Each shell entry has p primitives that look like
            // [ angmom, [ [ e1, c1 ], [ e2, c2 ], ...., [ ep, cp ] ] ]
            py::list shellinfo = atominfo[atomshells].cast<py::list>();
            int am = shellinfo[0].cast<int>();
            std::vector<double> coefficients;
            std::vector<double> exponents;
            int nprim = (pybind11::len(shellinfo)) - 1;  // The leading entry is the angular momentum
            for (int primitive = 1; primitive <= nprim; primitive++) {
                py::list primitiveinfo = shellinfo[primitive].cast<py::list>();
                exponents.push_back(primitiveinfo[0].cast<double>());
                coefficients.push_back(primitiveinfo[1].cast<double>());
            }
            vec_shellinfo.push_back(ShellInfo(am, coefficients, exponents, shelltype, Unnormalized));
        }
        mol->set_shell_by_label(atomlabel, hash, key);
        basis_atom_shell[name][atomlabel] = vec_shellinfo;
    }

    /*
     * Handle the ECP terms, if needed
     */
    map_ssv basis_atom_ecpshell;
    std::map<std::string, std::map<std::string, int>> basis_atom_ncore;
    // basisname is uniform; fill map with key/value (gbs entry) pairs of elements from pybs['shell_map']
    int totalncore = 0;
    if (pybs.contains("ecp_shell_map")) {
        py::list ecpbasisinfo = pybs["ecp_shell_map"].cast<py::list>();
        for (int atom = 0; atom < py::len(ecpbasisinfo); ++atom) {
            std::vector<ShellInfo> vec_shellinfo;
            py::list atominfo = ecpbasisinfo[atom].cast<py::list>();
            std::string atomlabel = atominfo[0].cast<std::string>();
            std::string hash = atominfo[1].cast<std::string>();
            int ncore = atominfo[2].cast<int>();
            for (int atomshells = 3; atomshells < py::len(atominfo); ++atomshells) {
                // Each shell entry has p primitives that look like
                // [ angmom, [ [ e1, c1, r1 ], [ e2, c2, r2 ], ...., [ ep, cp, rp ] ] ]
                py::list shellinfo = atominfo[atomshells].cast<py::list>();
                int am = shellinfo[0].cast<int>();
                std::vector<double> coefficients;
                std::vector<double> exponents;
                std::vector<int> ns;
                int nprim = (pybind11::len(shellinfo)) - 1;  // The leading entry is the angular momentum
                for (int primitive = 1; primitive <= nprim; primitive++) {
                    py::list primitiveinfo = shellinfo[primitive].cast<py::list>();
                    exponents.push_back(primitiveinfo[0].cast<double>());
                    coefficients.push_back(primitiveinfo[1].cast<double>());
                    ns.push_back(primitiveinfo[2].cast<int>());
                }
                vec_shellinfo.push_back(ShellInfo(am, coefficients, exponents, ns));
            }
            basis_atom_ncore[name][atomlabel] = ncore;
            basis_atom_ecpshell[name][atomlabel] = vec_shellinfo;
            totalncore += ncore;
        }
    }

    mol->update_geometry();  // update symmetry with basisset info

    auto basisset = std::make_shared<BasisSet>(key, mol, basis_atom_shell, basis_atom_ecpshell);

    // Modify the nuclear charges, to account for the ECP.
    if (totalncore) {
        for (int atom = 0; atom < mol->natom(); ++atom) {
            if (mol->Z(atom) > 0) {
                const std::string& basis = mol->basis_on_atom(atom);
                const std::string& label = mol->label(atom);
                int ncore = basis_atom_ncore[basis][label];
                int Z = mol->true_atomic_number(atom) - ncore;
                mol->set_nuclear_charge(atom, Z);
                basisset->set_n_ecp_core(label, ncore);
            }
        }
    }

    basisset->set_name(name);
    basisset->set_key(key);
    basisset->set_target(label);

    return basisset;
}

bool _has_key(const py::dict& data, const std::string& key) {
    for (auto item : data) {
        if (std::string(py::str(item.first)) == key) return true;
    }
    return false;
}

std::shared_ptr<Molecule> from_dict(py::dict molrec) {
    // Compromises for psi4.core.Molecule
    // * molecular_charge is int, not float
    // * fragment_charges are int, not float
    // * elez are float, not int, b/c Psi4 Z is approx. elez * real

    std::shared_ptr<Molecule> mol(new Molecule);
    mol->set_lock_frame(false);

    if (_has_key(molrec, "name")) mol->set_name(molrec["name"].cast<std::string>());

    if (molrec["units"].cast<std::string>() == "Angstrom")
        mol->set_units(Molecule::Angstrom);
    else if (molrec["units"].cast<std::string>() == "Bohr")
        mol->set_units(Molecule::Bohr);
    else
        throw PSIEXCEPTION("Invalid geometry units to construct Molecule.");
    if (_has_key(molrec, "input_units_to_au")) mol->set_input_units_to_au(molrec["input_units_to_au"].cast<double>());

    mol->set_com_fixed(molrec["fix_com"].cast<bool>());
    mol->set_orientation_fixed(molrec["fix_orientation"].cast<bool>());
    if (_has_key(molrec, "fix_symmetry")) mol->reset_point_group(molrec["fix_symmetry"].cast<std::string>());

    std::vector<int> elea = molrec["elea"].cast<std::vector<int>>();
    std::vector<double> elez = molrec["elez"].cast<std::vector<double>>();
    std::vector<std::string> elem = molrec["elem"].cast<std::vector<std::string>>();
    std::vector<double> mass = molrec["mass"].cast<std::vector<double>>();
    std::vector<int> real = molrec["real"].cast<std::vector<int>>();
    std::vector<std::string> elbl = molrec["elbl"].cast<std::vector<std::string>>();

    size_t nat;
    bool unsettled;
    if (_has_key(molrec, "geom_unsettled")) {
        std::vector<std::vector<std::string>> geom_unsettled =
            molrec["geom_unsettled"].cast<std::vector<std::vector<std::string>>>();
        nat = geom_unsettled.size();
        unsettled = true;

        for (size_t iat = 0; iat < nat; ++iat) {
            std::string symbol = elem.at(iat);
            std::string label = elbl.at(iat);
            std::transform(symbol.begin(), symbol.end(), symbol.begin(), ::toupper);
            std::transform(label.begin(), label.end(), label.begin(), ::toupper);
            mol->add_unsettled_atom(elez.at(iat) * real.at(iat), geom_unsettled.at(iat), symbol, mass.at(iat),
                                    elez.at(iat) * real.at(iat), symbol + label, elea.at(iat));
        }

        std::vector<std::pair<std::string, double>> variables =
            molrec["variables"].cast<std::vector<std::pair<std::string, double>>>();
        for (size_t iv = 0; iv < variables.size(); ++iv)
            mol->set_geometry_variable(variables[iv].first, variables[iv].second);

    } else {
        std::vector<double> geom = molrec["geom"].cast<std::vector<double>>();
        nat = geom.size() / 3;
        unsettled = false;

        for (size_t iat = 0; iat < nat; ++iat) {
            std::string symbol = elem.at(iat);
            std::string label = elbl.at(iat);
            std::transform(symbol.begin(), symbol.end(), symbol.begin(), ::toupper);
            std::transform(label.begin(), label.end(), label.begin(), ::toupper);
            mol->add_atom(elez.at(iat) * real.at(iat), geom.at(3 * iat), geom.at(3 * iat + 1), geom.at(3 * iat + 2),
                          symbol, mass.at(iat), elez.at(iat) * real.at(iat), symbol + label, elea.at(iat));
        }
    }

    std::vector<Molecule::FragmentType> fragment_types;
    std::vector<int> fragment_separators;
    fragment_separators = molrec["fragment_separators"].cast<std::vector<int>>();
    std::vector<std::pair<int, int>> fragments;
    fragment_separators.insert(fragment_separators.begin(), 0);
    fragment_separators.push_back(nat);
    for (size_t i = 1; i < fragment_separators.size(); ++i) {
        fragments.push_back(std::make_pair(fragment_separators[i - 1], fragment_separators[i]));
        fragment_types.push_back(Molecule::Real);
    }

    std::vector<int> fragment_charges;
    for (auto item : molrec["fragment_charges"]) fragment_charges.push_back(static_cast<int>(item.cast<double>()));

    mol->set_fragment_pattern(fragments, fragment_types, fragment_charges,
                              molrec["fragment_multiplicities"].cast<std::vector<int>>());

    mol->set_molecular_charge(static_cast<int>(molrec["molecular_charge"].cast<double>()));
    mol->set_multiplicity(molrec["molecular_multiplicity"].cast<int>());

    // hack to prevent update_geometry termination upon no atoms
    if (nat == 0) mol->set_lock_frame(true);

    if (!unsettled) mol->update_geometry();
    return mol;
}

void export_mints(py::module& m) {
    // This is needed to wrap an STL vector into Boost.Python. Since the vector
    // is going to contain std::shared_ptr's we MUST set the no_proxy flag to true
    // (as it is) to tell Boost.Python to not create a proxy class to handle
    // the vector's data type.
    py::bind_vector<std::vector<std::shared_ptr<Matrix>>>(m, "VectorMatrix");

    // Other vector types
    // py::class_<std::vector<double> >(m, "vector_of_doubles", "docstring").
    //        def(vector_indexing_suite<std::vector<double>, true >());
    //    py::bind_vector<double>(m, "VectorDouble");

    // Use typedefs to explicitly tell Boost.Python which function in the class
    // to use. In most cases, you should not be making Python specific versions
    // of functions.

    // For example in Vector there are 2 versions of set: a (double*) version and a
    // (int, int, double) version. We create a typedef function pointer to tell
    // Boost.Python we only want the (int, int, double) version.
    typedef void (Vector::*vector_setitem_1)(int, double);
    typedef void (Vector::*vector_setitem_2)(int, int, double);
    typedef double (Vector::*vector_getitem_1)(int);
    typedef double (Vector::*vector_getitem_2)(int, int);
    typedef void (Vector::*vector_setitem_n)(const py::tuple&, double);
    typedef double (Vector::*vector_getitem_n)(const py::tuple&);

    py::class_<Dimension>(m, "Dimension", "Initializes and defines Dimension Objects")
        .def(py::init<int>())
        .def(py::init<int, const std::string&>())
        .def(py::init<const std::vector<int>&>())
        .def("print_out", &Dimension::print, "Print out the dimension object to the output file")
        .def("init", &Dimension::init, "Re-initializes the dimension object")
        .def("sum", &Dimension::sum, "Gets the sum of the values in the dimension object")
        .def("max", &Dimension::max, "Gets the maximum value from the dimension object")
        .def("zero", &Dimension::zero, "Zeros all values in the dimension object")
        .def("n", &Dimension::n,
             // py::return_value_policy<copy_const_reference>(),
             py::return_value_policy::copy, "The order of the dimension")
        .def_property("name",
                      // make_function(&Dimension::name, return_value_policy<copy_const_reference>()),
                      py::cpp_function(&Dimension::name), py::cpp_function(&Dimension::set_name),
                      "The name of the dimension. Used in printing.")
        .def("sum", &Dimension::sum, "Return the sum of constituent dimensions")
        .def("max", &Dimension::max, "Return the maximum element")
        .def("zero", &Dimension::zero, "Zero all elements")
        .def("fill", &Dimension::fill, "Fill all elements with given value", py::arg("val"))
        .def(py::self += py::self)
        .def(py::self + py::self)
        .def(py::self -= py::self)
        .def(py::self - py::self)
        .def(py::self == py::self)
        .def(py::self != py::self)
        .def("__getitem__", &Dimension::get, py::return_value_policy::copy, "Get the i'th value", py::arg("i"))
        .def("__setitem__", &Dimension::set, "Set element i to value val", py::arg("i"), py::arg("val"));

    py::class_<Slice>(m, "Slice", "Slicing for Matrix and Vector objects")
        .def(py::init<Dimension&, Dimension&>())
        .def("begin", &Slice::begin, "Get the first element of this slice")
        .def("end", &Slice::end, "Get the past-the-end element of this slice");

    py::class_<Vector, std::shared_ptr<Vector>>(m, "Vector", "Class for creating and manipulating vectors",
                                                py::dynamic_attr())
        .def(py::init<int>())
        .def(py::init<const Dimension&>())
        .def(py::init<const std::string&, int>())
        .def(py::init<const std::string&, const Dimension&>())
        .def_property("name", py::cpp_function(&Vector::name), py::cpp_function(&Vector::set_name),
                      "The name of the Vector. Used in printing.")
        .def("get", vector_getitem_1(&Vector::get), "Returns a single element value located at m", py::arg("m"))
        .def("get", vector_getitem_2(&Vector::get), "Returns a single element value located at m in irrep h",
             py::arg("h"), py::arg("m"))
        .def("set", vector_setitem_1(&Vector::set), "Sets a single element value located at m", py::arg("m"),
             py::arg("val"))
        .def("set", vector_setitem_2(&Vector::set), "Sets a single element value located at m in irrep h", py::arg("h"),
             py::arg("m"), py::arg("val"))
        .def("print_out", &Vector::print_out, "Prints the vector to the output file")
        .def("scale", &Vector::scale, "Scales the elements of a vector by sc", py::arg("sc"))
        .def("dim", &Vector::dim, "Returns the dimensions of the vector per irrep h", py::arg("h"))
        .def("nirrep", &Vector::nirrep, "Returns the number of irreps")
        .def("get_block", &Vector::get_block, "Get a vector block", py::arg("slice"))
        .def("set_block", &Vector::set_block, "Set a vector block", py::arg("slice"), py::arg("block"))
        .def("array_interface",
             [](Vector& v) {

                 // Build a list of NumPy views, used for the .np and .nph accessors.Vy
                 py::list ret;

                 // If we set a NumPy shape
                 if (v.numpy_shape().size()) {
                     if (v.nirrep() > 1) {
                         throw PSIEXCEPTION(
                             "Vector::array_interface numpy shape with more than one irrep is not "
                             "valid.");
                     }

                     // Cast the NumPy shape vector
                     std::vector<size_t> shape;
                     for (int val : v.numpy_shape()) {
                         shape.push_back((size_t)val);
                     }

                     // Build the array
                     py::array arr(shape, v.pointer(0), py::cast(&v));
                     ret.append(arr);

                 } else {
                     for (size_t h = 0; h < v.nirrep(); h++) {
                         // Hmm, sometimes we need to handle empty ptr's correctly
                         double* ptr = nullptr;
                         if (v.dim(h) != 0) {
                             ptr = v.pointer(h);
                         }

                         // Build the array
                         std::vector<size_t> shape{(size_t)v.dim(h)};
                         py::array arr(shape, ptr, py::cast(&v));
                         ret.append(arr);
                     }
                 }

                 return ret;
             },
             py::return_value_policy::reference_internal);

    typedef void (IntVector::*int_vector_set)(int, int, int);
    py::class_<IntVector, std::shared_ptr<IntVector>>(m, "IntVector", "Class handling vectors with integer values")
        .def(py::init<int>())
        .def("get", &IntVector::get, "Returns a single element value located at m in irrep h", py::arg("h"),
             py::arg("m"))
        .def("set", int_vector_set(&IntVector::set), "Sets a single element value located at m in irrep h",
             py::arg("h"), py::arg("m"), py::arg("val"))
        .def("print_out", &IntVector::print_out, "Prints the vector to the output file")
        .def("dim", &IntVector::dim, "Returns the number of dimensions per irrep h", py::arg("h"))
        .def("nirrep", &IntVector::nirrep, "Returns the number of irreps");

    py::enum_<diagonalize_order>(m, "DiagonalizeOrder", "Defines ordering of eigenvalues after diagonalization")
        .value("Ascending", ascending)
        .value("Descending", descending)
        .export_values();

    py::enum_<Molecule::GeometryUnits>(m, "GeometryUnits", "The units used to define the geometry")
        .value("Angstrom", Molecule::Angstrom)
        .value("Bohr", Molecule::Bohr)
        .export_values();

    py::enum_<Molecule::FragmentType>(m, "FragmentType", "Fragment activation status")
        .value("Absent", Molecule::Absent)  // Neglect completely
        .value("Real", Molecule::Real)      // Include, as normal
        .value("Ghost", Molecule::Ghost)    // Include, but with ghost atoms
        .export_values();

    typedef void (Matrix::*matrix_multiply)(bool, bool, double, const SharedMatrix&, const SharedMatrix&, double);
    typedef void (Matrix::*matrix_diagonalize)(SharedMatrix&, std::shared_ptr<Vector>&, diagonalize_order);
    typedef void (Matrix::*matrix_one)(const SharedMatrix&);
    typedef double (Matrix::*double_matrix_one)(const SharedMatrix&);
    typedef void (Matrix::*matrix_two)(const SharedMatrix&, const SharedMatrix&);
    typedef void (Matrix::*matrix_save)(const std::string&, bool, bool, bool);
    typedef void (Matrix::*matrix_set1)(double);
    typedef void (Matrix::*matrix_set3)(int, int, double);
    typedef void (Matrix::*matrix_set4)(int, int, int, double);
    typedef double (Matrix::*matrix_get3)(const int&, const int&, const int&) const;
    typedef double (Matrix::*matrix_get2)(const int&, const int&) const;
    typedef void (Matrix::*matrix_load)(const std::string&);
    typedef bool (Matrix::*matrix_load_psio1)(std::shared_ptr<psi::PSIO>&, size_t, const std::string&, int);
    typedef void (Matrix::*matrix_load_psio2)(std::shared_ptr<psi::PSIO>&, size_t, Matrix::SaveType);
    typedef const Dimension& (Matrix::*matrix_ret_dimension)() const;

    py::enum_<Matrix::SaveType>(m, "SaveType", "The layout of the matrix for saving")
        .value("Full", Matrix::SaveType::Full)
        .value("SubBlocks", Matrix::SaveType::SubBlocks)
        .value("LowerTriangle", Matrix::SaveType::LowerTriangle)
        .export_values();

    py::class_<Matrix, std::shared_ptr<Matrix>>(m, "Matrix", "Class for creating and manipulating matrices",
                                                py::dynamic_attr(), py::buffer_protocol())
        .def(py::init<int, int>())
        .def(py::init<const std::string&, int, int>())
        .def(py::init<const std::string&, const Dimension&, const Dimension&>())
        .def(py::init<const std::string&>())
        .def("clone", &Matrix::clone, "Creates exact copy of the matrix and returns it")
        .def_property("name", py::cpp_function(&Matrix::name), py::cpp_function(&Matrix::set_name),
                      "The name of the Matrix. Used in printing.")

        // def("set_name", &Matrix::set_name, "docstring").
        // def("name", &Matrix::name, py::return_value_policy::copy, "docstring").
        .def("print_out", &Matrix::print_out, "Prints the matrix to the output file")
        .def("print_atom_vector", &Matrix::print_atom_vector, py::arg("RMRoutfile") = "outfile",
             "Print the matrix with atom labels, assuming it is an natom X 3 tensor")
        .def("rows", &Matrix::rowdim, "Returns the rows in irrep h", py::arg("h") = 0)
        .def("cols", &Matrix::coldim, "Returns the columns in irrep h", py::arg("h") = 0)
        .def("rowdim", matrix_ret_dimension(&Matrix::rowspi), py::return_value_policy::copy,
             "Returns the rows per irrep array")
        .def("coldim", matrix_ret_dimension(&Matrix::colspi), py::return_value_policy::copy,
             "Returns the columns per irrep array")
        .def("nirrep", &Matrix::nirrep, py::return_value_policy::copy, "Returns the number of irreps")
        .def("symmetry", &Matrix::symmetry, py::return_value_policy::copy, "Returns the overall symmetry of the matrix")
        .def("identity", &Matrix::identity, "Sets the matrix to the identity")
        .def("hermitize", &Matrix::hermitivitize,
             "Makes a real matrix symmetric by averaging the matrix and its transpose.")
        .def("copy_lower_to_upper", &Matrix::copy_lower_to_upper, "Copy the lower triangle to the upper triangle")
        .def("copy_upper_to_lower", &Matrix::copy_upper_to_lower, "Copy the upper triangle to the lower triangle")
        .def("zero_lower", &Matrix::zero_lower, "Zero the lower triangle")
        .def("zero_upper", &Matrix::zero_upper, "Zero the upper triangle")
        .def("zero", &Matrix::zero, "Zero all elements of the matrix")
        .def("zero_diagonal", &Matrix::zero_diagonal, "Zero the diagonal of the matrix")
        .def("trace", &Matrix::trace, "Returns the trace of the matrix")

        .def("transpose_this", &Matrix::transpose_this, "Transpose the matrix in-place")
        .def("transpose", &Matrix::transpose, "Creates a new matrix that is the transpose of this matrix")
        .def("hermitivitize", &Matrix::hermitivitize, "Average off-diagonal element in-place")
        .def("add", matrix_one(&Matrix::add), "Adds a matrix to this matrix")
        .def("add", matrix_set4(&Matrix::add), "Increments row m and column n of irrep h's block matrix by val.",
             py::arg("h"), py::arg("m"), py::arg("n"), py::arg("val"))
        .def("axpy", &Matrix::axpy, "Add to this matrix another matrix scaled by a", py::arg("a"), py::arg("X"))
        .def("subtract", matrix_one(&Matrix::subtract), "Substract a matrix from this matrix")
        .def("accumulate_product", matrix_two(&Matrix::accumulate_product),
             "Multiplies two arguments and adds the result to this matrix")
        .def("scale", &Matrix::scale, "Scales the matrix by the floating point value a", py::arg("a"))
        .def("sum_of_squares", &Matrix::sum_of_squares, "Returns the sum of the squares of this matrix")
        .def("add_and_orthogonalize_row", &Matrix::add_and_orthogonalize_row,
             "Expands the row dimension by one, \
              and then orthogonalizes vector v against the current rows \
              before setting the new row to the orthogonalized copy of v",
             py::arg("v"))
        .def("rms", &Matrix::rms, "Returns the rms of this matrix")
        .def("absmax", &Matrix::absmax, "Returns the absolute maximum value")
        .def("scale_row", &Matrix::scale_row, "Scales row m of irrep h by a", py::arg("h"), py::arg("m"), py::arg("a"))
        .def("scale_column", &Matrix::scale_column, "Scales column n of irrep h by a", py::arg("h"), py::arg("n"),
             py::arg("a"))
        .def("transform", matrix_one(&Matrix::transform), "Transform this matrix with transformer",
             py::arg("transformer"))
        .def("transform", matrix_two(&Matrix::transform), "Transform A with transformer", py::arg("a"),
             py::arg("transformer"))
        .def("back_transform", matrix_one(&Matrix::back_transform), "Backtransform this with transformer",
             py::arg("transformer"))
        .def("back_transform", matrix_two(&Matrix::back_transform), "Backtransform A with transformer", py::arg("a"),
             py::arg("transformer"))
        .def("vector_dot", double_matrix_one(&Matrix::vector_dot), "Returns the vector dot product of this with rhs",
             py::arg("rhs"))
        .def("gemm", matrix_multiply(&Matrix::gemm),
             "Generalized matrix multiplication \
              argument transa Transpose the left matrix? \
              argument transb Transpose the right matrix? \
              argument alpha Prefactor for the matrix multiplication \
              argument A Left matrix \
              argument B Right matrix \
              argument beta Prefactor for the resulting matrix",
             py::arg("transa"), py::arg("transb"), py::arg("alpha"), py::arg("a"), py::arg("b"), py::arg("beta"))
        .def("diagonalize", matrix_diagonalize(&Matrix::diagonalize),
             "Diagonalizes this matrix, space for the eigvectors and eigvalues must be created by caller. Only for "
             "symmetric matrices.",
             py::arg("eigvectors"), py::arg("eigvalues"), py::arg("order") = ascending)
        .def("cholesky_factorize", &Matrix::cholesky_factorize,
             "Computes the Cholesky factorization of a real symmetric positive definite matrix")
        .def(
            "partial_cholesky_factorize", &Matrix::partial_cholesky_factorize,
            "Computes the fully pivoted partial Cholesky factorization of a real symmetric positive semidefinite matrix, \
              to numerical precision delta",
            py::arg("delta") = 0.0, py::arg("throw_if_negative") = false)

        // def("canonical_orthogonalization", &Matrix::canonical_orthogonalization,
        // CanonicalOrthog()).
        // def("canonical_orthogonalization", &Matrix::canonical_orthogonalization, py::arg("delta")
        // = 0.0, py::arg("eigvec") = SharedMatrix()).
        .def("schmidt", &Matrix::schmidt, "Calls the libqt schmidt function")
        .def("invert", &Matrix::invert, "Computes the inverse of a real symmetric positive definite matrix")
        .def("general_invert", &Matrix::general_invert,
             "Computes the inverse of any nonsingular matrix using LU factorization")
        .def("pseudoinverse", &Matrix::pseudoinverse,
             "Computes the matrix which is the conditioned pseudoinverse of this matrix", py::arg("condition"),
             py::arg("nremoved"))
        .def("apply_denominator", matrix_one(&Matrix::apply_denominator), "Apply matrix of denominators to this matrix",
             py::arg("Matrix"))
        .def("copy", matrix_one(&Matrix::copy), "Returns a copy of the matrix")
        .def("power", &Matrix::power, "Takes the matrix to the alpha power with precision cutoff", py::arg("alpha"),
             py::arg("cutoff") = 1.0E-12)
        .def_static("doublet", &Matrix::doublet,
                    "Returns the multiplication of two matrices A and B, with options to transpose each beforehand",
                    py::arg("A"), py::arg("B"), py::arg("transA") = false, py::arg("transB") = false)
        .def_static(
            "triplet", &Matrix::triplet,
            "Returns the multiplication of three matrics A, B, and C, with options to transpose each beforehand",
            py::arg("A"), py::arg("B"), py::arg("C"), py::arg("transA") = false, py::arg("transB") = false,
            py::arg("transC") = false)
        .def("get", matrix_get3(&Matrix::get), "Returns a single element of a matrix in subblock h, row m, col n",
             py::arg("h"), py::arg("m"), py::arg("n"))
        .def("get", matrix_get2(&Matrix::get), "Returns a single element of a matrix, row m, col n", py::arg("m"),
             py::arg("n"))
        .def("get_block", &Matrix::get_block, "Get a matrix block", py::arg("rows"), py::arg("cols"))
        .def("set", matrix_set1(&Matrix::set), "Sets every element of a matrix to val", py::arg("val"))
        .def("set", matrix_set3(&Matrix::set), "Sets a single element of a matrix to val at row m, col n", py::arg("m"),
             py::arg("n"), py::arg("val"))
        .def("set", matrix_set4(&Matrix::set),
             "Sets a single element of a matrix, subblock h, row m, col n, with value val", py::arg("h"), py::arg("m"),
             py::arg("n"), py::arg("val"))
        .def("set_block", &Matrix::set_block, "Set a matrix block", py::arg("rows"), py::arg("cols"), py::arg("block"))
        // destroyed according to matrix.h file
        //.def("project_out", &Matrix::project_out, "docstring")
        .def("save", matrix_save(&Matrix::save),
             "Saves the matrix in ASCII format to filename, as symmetry blocks or full matrix", py::arg("filename"),
             py::arg("append") = true, py::arg("saveLowerTriangle") = true, py::arg("saveSubBlocks") = false)
        .def("load", matrix_load(&Matrix::load),
             "Loads a block matrix from an ASCII file (see tests/mints3 for format)", py::arg("filename"))
        .def("load_mpqc", &Matrix::load_mpqc, "Loads a matrix from an ASCII file in MPQC format", py::arg("filename"))
        .def("load", matrix_load_psio1(&Matrix::load),
             "Load a matrix from a PSIO object from fileno with tocentry of size nso", py::arg("psio"),
             py::arg("fileno"), py::arg("tocentry"), py::arg("nso"))
        .def("load", matrix_load_psio2(&Matrix::load),
             "Load a matrix from a PSIO object from fileno and with toc position of the name of the matrix",
             py::arg("psio"), py::arg("fileno"), py::arg("savetype") = Matrix::SaveType::LowerTriangle)
        // should this take Petite List's sotoao() function as a default transfomer argument? has to be set C++ side
        // first i think
        .def("remove_symmetry", &Matrix::remove_symmetry, "Remove symmetry from a matrix A with PetiteList::sotoao()",
             py::arg("a"), py::arg("transformer"))
        .def("symmetrize_gradient", &Matrix::symmetrize_gradient,
             "Symmetrizes a gradient-like matrix (N,3) using information from a given molecule", py::arg("mol"))
        .def("rotate_columns", &Matrix::rotate_columns, "Rotates columns i and j in irrep h by angle theta",
             py::arg("h"), py::arg("i"), py::arg("j"), py::arg("theta"))
        .def("array_interface",
             [](Matrix& m) {

                 // Build a list of NumPy views, used for the .np and .nph accessors.Vy
                 py::list ret;

                 // If we set a NumPy shape
                 if (m.numpy_shape().size()) {
                     if (m.nirrep() > 1) {
                         throw PSIEXCEPTION(
                             "Vector::array_interface numpy shape with more than one irrep is not "
                             "valid.");
                     }

                     // Cast the NumPy shape vector
                     std::vector<size_t> shape;
                     for (int val : m.numpy_shape()) {
                         shape.push_back((size_t)val);
                     }

                     // Build the array
                     py::array arr(shape, m.pointer(0)[0], py::cast(&m));
                     ret.append(arr);

                 } else {
                     for (size_t h = 0; h < m.nirrep(); h++) {
                         // Hmm, sometimes we need to overload to nullptr
                         double* ptr = nullptr;
                         if ((m.rowdim(h) * m.coldim(h)) != 0) {
                             ptr = m.pointer(h)[0];
                         }

                         // Build the array
                         py::array arr({(size_t)m.rowdim(h), (size_t)m.coldim(h)}, ptr, py::cast(&m));
                         ret.append(arr);
                     }
                 }

                 return ret;
             },
             py::return_value_policy::reference_internal);

    py::class_<Deriv, std::shared_ptr<Deriv>>(m, "Deriv", "Computes gradients of wavefunctions")
        .def(py::init<std::shared_ptr<Wavefunction>>())
        .def(py::init<std::shared_ptr<Wavefunction>, char, bool, bool>())
        .def("set_tpdm_presorted", &Deriv::set_tpdm_presorted, "Is the TPDM already presorted? Default is False",
             py::arg("val") = false)
        .def("set_ignore_reference", &Deriv::set_ignore_reference,
             "Ignore reference contributions to the gradient? Default is False", py::arg("val") = false)
        .def("set_deriv_density_backtransformed", &Deriv::set_deriv_density_backtransformed,
             "Is the deriv_density already backtransformed? Default is False", py::arg("val") = false)
        .def("compute", &Deriv::compute, "Compute the gradient");

    typedef SharedMatrix (MatrixFactory::*create_shared_matrix)() const;
    typedef SharedMatrix (MatrixFactory::*create_shared_matrix_name)(const std::string&) const;
    // Something here is wrong. These will not work with py::args defined, not sure why
    py::class_<MatrixFactory, std::shared_ptr<MatrixFactory>>(m, "MatrixFactory", "Creates Matrix objects")
        .def("create_matrix", create_shared_matrix(&MatrixFactory::create_shared_matrix),
             "Returns a new matrix object with default dimensions")
        //             "Returns a new matrix object with default dimensions", py::arg("symmetry"));
        .def("create_matrix", create_shared_matrix_name(&MatrixFactory::create_shared_matrix),
             "Returns a new Matrix object named name with default dimensions");
    //              py::arg("name"), py::arg("symmetry"));

    py::class_<CdSalc::Component, std::shared_ptr<CdSalc::Component>>(
        m, "SalcComponent", "Component of a Cartesian displacement SALC in the basis of atomic displacements.")
        .def_readwrite("coef", &CdSalc::Component::coef, "The coefficient of the displacement")
        .def_readwrite("atom", &CdSalc::Component::atom, "The index of the atom being displaced. 0-indexed.")
        .def_readwrite("xyz", &CdSalc::Component::xyz,
                       "The direction of the displacement, given by x as 0, y as 1, z as 2.");

    py::class_<CdSalc, std::shared_ptr<CdSalc>>(m, "CdSalc", "Cartesian displacement SALC")
        .def("irrep", &CdSalc::irrep, "Return the irrep bit representation")
        .def("irrep_index", [](const CdSalc& salc) { return static_cast<int>(salc.irrep()); }, "Return the irrep index")
        .def("print_out", &CdSalc::print,
             "Print the irrep index and the coordinates of the SALC of Cartesian displacements. \
                                           Irrep index is 0-indexed and Cotton ordered.")
        .def("__getitem__", [](const CdSalc& salc, size_t i) { return salc.component(i); })
        .def("__len__", [](const CdSalc& salc) { return salc.ncomponent(); })
        .def("__iter__", [](const CdSalc& salc) { return py::make_iterator(salc.get_components()); },
             py::keep_alive<0, 1>());

    py::class_<CdSalcList, std::shared_ptr<CdSalcList>>(
        m, "CdSalcList", "Class for generating symmetry adapted linear combinations of Cartesian displacements")
        .def(py::init<std::shared_ptr<Molecule>, int, bool, bool>())
        .def("ncd", &CdSalcList::ncd, "Return the number of cartesian displacements SALCs")
        .def("create_matrices", &CdSalcList::create_matrices,
             "Return a vector of matrices with the SALC symmetries. Dimensions determined by factory.",
             py::arg("basename"), py::arg("factory"))
        .def("salc_name", &CdSalcList::salc_name, "Return the name of SALC #i.", py::arg("i"))
        .def("nirrep", &CdSalcList::nirrep, "Return the number of irreps")
        .def("__getitem__", [](const CdSalcList& salclist, size_t i) { return salclist[i]; })
        .def("__len__", [](const CdSalcList& salclist) { return salclist.ncd(); })
        .def("__iter__", [](const CdSalcList& salclist) { return py::make_iterator(salclist.get_salcs()); },
             py::keep_alive<0, 1>())
        .def("print_out", &CdSalcList::print, "Print the SALCs to the output file")
        .def("matrix", &CdSalcList::matrix, "Return the matrix that transforms Cartesian displacements to SALCs")
        .def("matrix_irrep", &CdSalcList::matrix_irrep,
             "Return the matrix that transforms Cartesian displacements to SALCs of irrep h", py::arg("h"));

    py::class_<GaussianShell, std::shared_ptr<GaussianShell>>(m, "GaussianShell",
                                                              "Class containing information about basis functions")
        .def_property_readonly("nprimitive", py::cpp_function(&GaussianShell::nprimitive),
                               "The number of primitive gaussians")
        .def_property_readonly("nfunction", py::cpp_function(&GaussianShell::nfunction),
                               "Total number of basis functions")
        .def_property_readonly("ncartesian", py::cpp_function(&GaussianShell::ncartesian),
                               "Total number of basis functions if this shell was Cartesian")
        .def_property_readonly("am", py::cpp_function(&GaussianShell::am),
                               "The angular momentum of the given contraction")
        .def_property_readonly("amchar", py::cpp_function(&GaussianShell::amchar),
                               "The character symbol for the angular momentum of the given contraction")
        .def_property_readonly("AMCHAR", py::cpp_function(&GaussianShell::AMCHAR),
                               "The upper-case character symbol for the angular momentum of the given contraction")
        .def_property_readonly("ncenter", py::cpp_function(&GaussianShell::ncenter),
                               "Returns atom number this shell is on")
        .def_property("function_index", py::cpp_function(&GaussianShell::function_index),
                      py::cpp_function(&GaussianShell::set_function_index),
                      "Basis function index where this shell starts.")
        .
        //            add_property("center", &GaussianShell::center, "A double* representing the
        //            center of the GaussianShell.").
        //            add_property("exps", &GaussianShell::exps, "The exponents of all the
        //            primitives").
        //            add_property("coefs", &GaussianShell::coefs, "The coefficients of all the
        //            primitives").
        def("is_cartesian", &GaussianShell::is_cartesian, "Returns true if the contraction is Cartesian")
        .def("is_pure", &GaussianShell::is_pure,
             "Returns true if the contraction is pure, i.e. a spherical harmonic basis function")
        .
        //            def("normalize_shell", &GaussianShell::normalize_shell, "docstring").
        def("exp", &GaussianShell::exp, "Returns the exponent of the given primitive", py::arg("prim"))
        .def("original_coef", &GaussianShell::original_coef, "Return unnormalized coefficient of the pi'th primitive",
             py::arg("pi"))
        .def("erd_coef", &GaussianShell::erd_coef, "Return ERD normalized coefficient of pi'th primitive",
             py::arg("pi"))
        .def("coef", &GaussianShell::coef, "Return coefficient of the pi'th primitive", py::arg("pi"));

    py::enum_<PrimitiveType>(m, "PrimitiveType", "May be Normalized or Unnormalized")
        .value("Normalized", Normalized)
        .value("Unnormalized", Unnormalized)
        .export_values();

    py::enum_<GaussianType>(m, "GaussianType", "0 if Cartesian, 1 if Pure")
        .value("Cartesian", Cartesian)
        .value("Pure", Pure)
        .export_values();

    py::class_<ShellInfo, std::shared_ptr<ShellInfo>>(m, "ShellInfo")
        .def(py::init<int, const std::vector<double>&, const std::vector<double>&, GaussianType, PrimitiveType>());

    py::bind_vector<std::vector<ShellInfo>>(m, "BSVec");

    py::class_<OneBodyAOInt, std::shared_ptr<OneBodyAOInt>> pyOneBodyAOInt(
        m, "OneBodyAOInt", "Basis class for all one-electron integrals");
    pyOneBodyAOInt
        .def("compute_shell", &OneBodyAOInt::compute_shell,
             "Compute integrals between basis functions in the given shell pair")
        .def_property("origin", py::cpp_function(&OneBodyAOInt::origin), py::cpp_function(&OneBodyAOInt::set_origin),
                      "The origin about which the one body ints are being computed.")
        .def_property_readonly("basis", py::cpp_function(&OneBodyAOInt::basis), "The basis set on center one")
        .def_property_readonly("basis1", py::cpp_function(&OneBodyAOInt::basis1), "The basis set on center one")
        .def_property_readonly("basis2", py::cpp_function(&OneBodyAOInt::basis2),
                               "The basis set on center two");  // <-- Added semicolon

    // typedef void (OneBodySOInt::*matrix_version)(SharedMatrix) const;
    // typedef void (OneBodySOInt::*vector_version)(std::vector<SharedMatrix>) const;
    // class_<OneBodySOInt, std::shared_ptr<OneBodySOInt>, boost::noncopyable>("OneBodySOInt",
    // "docstring", no_init).
    //        def("compute", matrix_version(&OneBodySOInt::compute_shell), "docstring").
    //        def("compute_list", vector_version(&OneBodySOInt::compute), "docstring").
    //        add_property("basis", &OneBodySOInt::basis, "The basis set on center one").
    //        add_property("basis1", &OneBodySOInt::basis1, "The basis set on center one").
    //        add_property("basis2", &OneBodySOInt::basis2, "The basis set on center two");

    py::class_<OverlapInt, std::shared_ptr<OverlapInt>>(m, "OverlapInt", pyOneBodyAOInt, "Computes overlap integrals");
    py::class_<DipoleInt, std::shared_ptr<DipoleInt>>(m, "DipoleInt", pyOneBodyAOInt, "Computes dipole integrals");
    py::class_<QuadrupoleInt, std::shared_ptr<QuadrupoleInt>>(m, "QuadrupoleInt", pyOneBodyAOInt,
                                                              "Computes quadrupole integrals");
    py::class_<MultipoleInt, std::shared_ptr<MultipoleInt>>(m, "MultipoleInt", pyOneBodyAOInt,
                                                            "Computes arbitrary-order multipole integrals");
    py::class_<TracelessQuadrupoleInt, std::shared_ptr<TracelessQuadrupoleInt>>(
        m, "TracelessQuadrupoleInt", pyOneBodyAOInt, "Computes traceless quadrupole integrals");
    py::class_<ElectricFieldInt, std::shared_ptr<ElectricFieldInt>>(m, "ElectricFieldInt", pyOneBodyAOInt,
                                                                    "Computes electric field integrals");
    py::class_<KineticInt, std::shared_ptr<KineticInt>>(m, "KineticInt", pyOneBodyAOInt, "Computes kinetic integrals");
    py::class_<PotentialInt, std::shared_ptr<PotentialInt>>(m, "PotentialInt", pyOneBodyAOInt,
                                                            "Computes potential integrals");
    py::class_<PseudospectralInt, std::shared_ptr<PseudospectralInt>>(m, "PseudospectralInt", pyOneBodyAOInt,
                                                                      "Computes pseudospectral integrals");
    py::class_<ElectrostaticInt, std::shared_ptr<ElectrostaticInt>>(m, "ElectrostaticInt", pyOneBodyAOInt,
                                                                    "Computes electrostatic integrals");
    py::class_<NablaInt, std::shared_ptr<NablaInt>>(m, "NablaInt", pyOneBodyAOInt, "Computes nabla integrals");
    py::class_<AngularMomentumInt, std::shared_ptr<AngularMomentumInt>>(m, "AngularMomentumInt", pyOneBodyAOInt,
                                                                        "Computes angular momentum integrals");

    typedef size_t (TwoBodyAOInt::*compute_shell_ints)(int, int, int, int);
    py::class_<TwoBodyAOInt, std::shared_ptr<TwoBodyAOInt>> pyTwoBodyAOInt(m, "TwoBodyAOInt",
                                                                           "Two body integral base class");
    pyTwoBodyAOInt.def("compute_shell", compute_shell_ints(&TwoBodyAOInt::compute_shell),
                       "Compute ERIs between 4 shells");  // <-- Semicolon

    py::class_<TwoElectronInt, std::shared_ptr<TwoElectronInt>>(m, "TwoElectronInt", pyTwoBodyAOInt,
                                                                "Computes two-electron repulsion integrals")
        .def("compute_shell", compute_shell_ints(&TwoBodyAOInt::compute_shell), "Compute ERIs between 4 shells");

    py::class_<ERI, std::shared_ptr<ERI>>(m, "ERI", pyTwoBodyAOInt, "Computes normal two electron reuplsion integrals");
    py::class_<F12, std::shared_ptr<F12>>(m, "F12", pyTwoBodyAOInt, "Computes F12 electron repulsion integrals");
    py::class_<F12G12, std::shared_ptr<F12G12>>(m, "F12G12", pyTwoBodyAOInt,
                                                "Computes F12G12 electron repulsion integrals");
    py::class_<F12Squared, std::shared_ptr<F12Squared>>(m, "F12Squared", pyTwoBodyAOInt,
                                                        "Computes F12 Squared electron repulsion integrals");
    py::class_<F12DoubleCommutator, std::shared_ptr<F12DoubleCommutator>>(
        m, "F12DoubleCommutator", pyTwoBodyAOInt, "Computes F12 Double Commutator electron repulsion integrals");
    py::class_<ErfERI, std::shared_ptr<ErfERI>>(m, "ErfERI", pyTwoBodyAOInt,
                                                "Computes ERF electron repulsion integrals");
    py::class_<ErfComplementERI, std::shared_ptr<ErfComplementERI>>(
        m, "ErfComplementERI", pyTwoBodyAOInt, "Computes ERF complement electron repulsion integrals");

    py::class_<AOShellCombinationsIterator, std::shared_ptr<AOShellCombinationsIterator>>(m,
                                                                                          "AOShellCombinationsIterator")
        .def_property_readonly("p", py::cpp_function(&AOShellCombinationsIterator::p), "Returns current P index")
        .def_property_readonly("q", py::cpp_function(&AOShellCombinationsIterator::q), "Returns current Q index")
        .def_property_readonly("r", py::cpp_function(&AOShellCombinationsIterator::r), "Returns current R index")
        .def_property_readonly("s", py::cpp_function(&AOShellCombinationsIterator::s), "Returns current S index")
        .def("first", &AOShellCombinationsIterator::first, "docstring")
        .def("next", &AOShellCombinationsIterator::next, "docstring")  // these are not documented C++ side
        .def("is_done", &AOShellCombinationsIterator::is_done, "docstring");

    py::class_<ThreeCenterOverlapInt, std::shared_ptr<ThreeCenterOverlapInt>>(m, "ThreeCenterOverlapInt",
                                                                              "Three center overlap integrals")
        .def("compute_shell", &ThreeCenterOverlapInt::compute_shell, "Compute the integrals of the form (a|b|c)");

    py::class_<IntegralFactory, std::shared_ptr<IntegralFactory>>(m, "IntegralFactory", "Computes integrals")
        .def(py::init<std::shared_ptr<BasisSet>, std::shared_ptr<BasisSet>, std::shared_ptr<BasisSet>,
                      std::shared_ptr<BasisSet>>())
        .def(py::init<std::shared_ptr<BasisSet>>())
        // def("shells_iterator", &IntegralFactory::shells_iterator_ptr,
        // py::return_value_policy<manage_new_object>(), "docstring").
        .def("shells_iterator", &IntegralFactory::shells_iterator_ptr,
             "Returns an ERI iterator object, only coded for standard ERIs")
        .def("eri", &IntegralFactory::eri, "Returns an ERI integral object", py::arg("deriv") = 0,
             py::arg("use_shell_pairs") = true)
        .def("f12", &IntegralFactory::f12, "Returns an F12 integral object", py::arg("cf"), py::arg("deriv") = 0,
             py::arg("use_shell_pairs") = true)
        .def("f12g12", &IntegralFactory::f12g12, "Returns an F12G12 integral object", py::arg("cf"),
             py::arg("deriv") = 0, py::arg("use_shell_pairs") = true)
        .def("f12_double_commutator", &IntegralFactory::f12_double_commutator,
             "Returns an F12 double commutator integral object", py::arg("cf"), py::arg("deriv") = 0,
             py::arg("use_shell_pairs") = true)
        .def("f12_squared", &IntegralFactory::f12_squared, "Returns an F12 squared integral object", py::arg("cf"),
             py::arg("deriv") = 0, py::arg("use_shell_pairs") = true)
        .def("erf_eri", &IntegralFactory::erf_eri, "Returns and erf ERI integral object (omega integral)",
             py::arg("omega"), py::arg("deriv") = 0, py::arg("use_shell_pairs") = true)
        .def("erf_complement_eri", &IntegralFactory::erf_complement_eri,
             "Returns an erf complement ERI integral object (omega integral)", py::arg("omega"), py::arg("deriv") = 0,
             py::arg("use_shell_pairs") = true)
        .def("ao_overlap", &IntegralFactory::ao_overlap, "Returns a OneBodyInt that computes the AO overlap integrals",
             py::arg("deriv") = 0)
        .def("so_overlap", &IntegralFactory::so_overlap, "Returns a OneBodyInt that computes the SO overlap integrals",
             py::arg("deriv") = 0)
        .def("ao_dipole", &IntegralFactory::ao_dipole, "Returns a OneBodyInt that computes the AO dipole integrals",
             py::arg("deriv") = 0)
        .def("so_dipole", &IntegralFactory::so_dipole, "Returns a OneBodyInt that computes the SO dipole integrals",
             py::arg("deriv") = 0)
        .def("ao_kinetic", &IntegralFactory::ao_kinetic, "Returns a OneBodyInt that computes the AO kinetic integrals",
             py::arg("deriv") = 0)
        .def("so_kinetic", &IntegralFactory::so_kinetic, "Returns a OneBodyInt that computes the SO kinetic integrals",
             py::arg("deriv") = 0)
        .def("ao_potential", &IntegralFactory::ao_potential,
             "Returns a OneBodyInt that computes the AO nuclear attraction integral", py::arg("deriv") = 0)
        .def("so_potential", &IntegralFactory::so_potential,
             "Returns a OneBodyInt that computes the SO nuclear attraction integral", py::arg("deriv") = 0)
        .def("ao_pseudospectral", &IntegralFactory::ao_pseudospectral,
             "Returns a OneBodyInt that computes the AO pseudospectral grid integrals", py::arg("deriv") = 0)
        .def("so_pseudospectral", &IntegralFactory::so_pseudospectral,
             "Returns a OneBodyInt that computes the SO pseudospectral grid integrals", py::arg("deriv") = 0)
        .def("ao_nabla", &IntegralFactory::ao_nabla, "Returns a OneBodyInt that computes the AO nabla integral",
             py::arg("deriv") = 0)
        .def("so_nabla", &IntegralFactory::so_nabla, "Returns a OneBodyInt that computes the SO nabla integral",
             py::arg("deriv") = 0)
        .def("ao_angular_momentum", &IntegralFactory::ao_angular_momentum,
             "Returns a OneBodyInt that computes the AO angular momentum integral", py::arg("deriv") = 0)
        .def("so_angular_momentum", &IntegralFactory::so_angular_momentum,
             "Returns a OneBodyInt that computes the SO angular momentum integral", py::arg("deriv") = 0)
        .def("ao_quadrupole", &IntegralFactory::ao_quadrupole,
             "Returns a OneBodyInt that computes AO the quadrupole integral")
        .def("so_quadrupole", &IntegralFactory::so_quadrupole,
             "Returns a OneBodyInt that computes SO the quadrupole integral")
        .def("ao_multipoles", &IntegralFactory::ao_multipoles,
             "Returns a OneBodyInt that computes arbitrary-order AO multipole integrals", py::arg("order"))
        .def("so_multipoles", &IntegralFactory::so_multipoles,
             "Returns a OneBodyInt that computes arbitrary-order SO multipole integrals", py::arg("order"))
        .def("ao_traceless_quadrupole", &IntegralFactory::ao_traceless_quadrupole,
             "Returns a OneBodyInt that computes the traceless AO quadrupole integral")
        .def("so_traceless_quadrupole", &IntegralFactory::so_traceless_quadrupole,
             "Returns a OneBodyInt that computes the traceless SO quadrupole integral")
        .def("electric_field", &IntegralFactory::electric_field,
             "Returns a OneBodyInt that computes the electric field")
        .def("electrostatic", &IntegralFactory::electrostatic,
             "Returns a OneBodyInt that computes the point electrostatic potential")
        .def("overlap_3c", &IntegralFactory::overlap_3c,
             "Returns a OneBodyInt that computes the 3 center overlap integral");

    typedef std::shared_ptr<PetiteList> (MintsHelper::*petite_list_0)() const;
    typedef std::shared_ptr<PetiteList> (MintsHelper::*petite_list_1)(bool) const;

    typedef SharedMatrix (MintsHelper::*erf)(double, SharedMatrix, SharedMatrix, SharedMatrix, SharedMatrix);
    typedef SharedMatrix (MintsHelper::*eri)(SharedMatrix, SharedMatrix, SharedMatrix, SharedMatrix);
    typedef SharedMatrix (MintsHelper::*normal_eri)();
    typedef SharedMatrix (MintsHelper::*normal_eri_factory)(std::shared_ptr<IntegralFactory>);
    typedef SharedMatrix (MintsHelper::*normal_eri2)(std::shared_ptr<BasisSet>, std::shared_ptr<BasisSet>,
                                                     std::shared_ptr<BasisSet>, std::shared_ptr<BasisSet>);
    typedef SharedMatrix (MintsHelper::*normal_3c)(std::shared_ptr<BasisSet>, std::shared_ptr<BasisSet>,
                                                   std::shared_ptr<BasisSet>);

    typedef SharedMatrix (MintsHelper::*normal_f12)(std::shared_ptr<CorrelationFactor>);
    typedef SharedMatrix (MintsHelper::*normal_f122)(std::shared_ptr<CorrelationFactor>, std::shared_ptr<BasisSet>,
                                                     std::shared_ptr<BasisSet>, std::shared_ptr<BasisSet>,
                                                     std::shared_ptr<BasisSet>);

    typedef SharedMatrix (MintsHelper::*oneelectron)();
    typedef SharedMatrix (MintsHelper::*oneelectron_mixed_basis)(std::shared_ptr<BasisSet>, std::shared_ptr<BasisSet>);
    typedef SharedMatrix (MintsHelper::*perturb_grad_options)(SharedMatrix);
    typedef SharedMatrix (MintsHelper::*perturb_grad_xyz)(SharedMatrix, double, double, double);

    py::class_<MintsHelper, std::shared_ptr<MintsHelper>>(m, "MintsHelper", "Computes integrals")
        .def(py::init<std::shared_ptr<BasisSet>>())
        .def(py::init<std::shared_ptr<Wavefunction>>())

        // Options and attributes
        .def("nbf", &MintsHelper::nbf, "Returns the number of basis functions")
        .def("set_print", &MintsHelper::set_print, "Sets the print level")
        .def("basisset", &MintsHelper::basisset, "Returns the basis set being used")
        .def("sobasisset", &MintsHelper::sobasisset, "Returns the SO basis set being used")
        .def("factory", &MintsHelper::factory, "Returns the Matrix factory being used")
        .def("cdsalcs", &MintsHelper::cdsalcs, "Returns a CdSalcList object")
        .def("petite_list", petite_list_0(&MintsHelper::petite_list),
             "Returns petite list, which transforms AO basis functions to SO's")
        .def("petite_list1", petite_list_1(&MintsHelper::petite_list),
             "Returns petite list which transforms AO basis functions to SO's, \
              setting argument to true is for Cartesian basis, false is for Spherical Harmonic basis",
             py::arg("include_pure_transform"))

        // Integral builders
        .def("integral", &MintsHelper::integral, "Integral factory being used")
        .def("integrals", &MintsHelper::integrals, "Molecular integrals")
        .def("integrals_erf", &MintsHelper::integrals_erf, "ERF integrals", py::arg("w") = -1.0)
        .def("integrals_erfc", &MintsHelper::integrals_erfc, "ERFC integrals", py::arg("w") = -1.0)
        .def("one_electron_integrals", &MintsHelper::one_electron_integrals, "Standard one-electron integrals")

        // One-electron
        .def("ao_overlap", oneelectron(&MintsHelper::ao_overlap), "AO basis overlap integrals")
        .def("ao_overlap", oneelectron_mixed_basis(&MintsHelper::ao_overlap), "AO mixed basis overlap integrals")
        .def("so_overlap", &MintsHelper::so_overlap, "SO basis overlap integrals")
        .def("ao_kinetic", oneelectron(&MintsHelper::ao_kinetic), "AO basis kinetic integrals")
        .def("ao_kinetic", oneelectron_mixed_basis(&MintsHelper::ao_kinetic), "AO mixed basis kinetic integrals")
        .def("so_kinetic", &MintsHelper::so_kinetic, "SO basis kinetic integrals")
        .def("ao_potential", oneelectron(&MintsHelper::ao_potential), "AO potential integrals")
        .def("ao_potential", oneelectron_mixed_basis(&MintsHelper::ao_potential), "AO mixed basis potential integrals")
        .def("so_potential", &MintsHelper::so_potential, "SO basis potential integrals",
             py::arg("include_perturbations") = true)
        .def("ao_ecp", oneelectron(&MintsHelper::ao_ecp), "AO basis effective core potential integrals.")
        .def("ao_ecp", oneelectron_mixed_basis(&MintsHelper::ao_ecp), "AO basis effective core potential integrals.")
        .def("so_ecp", &MintsHelper::so_ecp, "SO basis effective core potential integrals.")

        // One-electron properties and
        .def("ao_pvp", &MintsHelper::ao_pvp, "AO pvp integrals")
        .def("ao_dkh", &MintsHelper::ao_dkh, "AO dkh integrals")
        .def("so_dkh", &MintsHelper::so_dkh, "SO dkh integrals")
        .def("ao_dipole", &MintsHelper::ao_dipole, "Vector AO dipole integrals")
        .def("so_dipole", &MintsHelper::so_dipole, "Vector SO dipole integrals")
        .def("ao_quadrupole", &MintsHelper::ao_quadrupole, "Vector AO quadrupole integrals")
        .def("so_quadrupole", &MintsHelper::so_quadrupole, "Vector SO quadrupole integrals")
        .def("ao_traceless_quadrupole", &MintsHelper::ao_traceless_quadrupole,
             "Vector AO traceless quadrupole integrals")
        .def("so_traceless_quadrupole", &MintsHelper::so_traceless_quadrupole,
             "Vector SO traceless quadrupole integrals")
        .def("ao_nabla", &MintsHelper::ao_nabla, "Vector AO nabla integrals")
        .def("so_nabla", &MintsHelper::so_nabla, "Vector SO nabla integrals")
        .def("ao_angular_momentum", &MintsHelper::ao_angular_momentum, "Vector AO angular momentum integrals")
        .def("so_angular_momentum", &MintsHelper::so_angular_momentum, "Vector SO angular momentum integrals")
        .def("ao_efp_multipole_potential", &MintsHelper::ao_efp_multipole_potential,
             "Vector AO EFP multipole integrals", py::arg("origin") = std::vector<double>{0, 0, 0},
             py::arg("deriv") = 0)
        .def("electric_field", &MintsHelper::electric_field, "Vector electric field integrals",
             py::arg("origin") = std::vector<double>{0, 0, 0}, py::arg("deriv") = 0)

        // Two-electron AO
        .def("ao_eri", normal_eri_factory(&MintsHelper::ao_eri), "AO ERI integrals", py::arg("factory") = nullptr)
        .def("ao_eri", normal_eri2(&MintsHelper::ao_eri), "AO ERI integrals", py::arg("bs1"), py::arg("bs2"),
             py::arg("bs3"), py::arg("bs4"))
        .def("ao_eri_shell", &MintsHelper::ao_eri_shell, "AO ERI Shell", py::arg("M"), py::arg("N"), py::arg("P"),
             py::arg("Q"))
        .def("ao_erf_eri", &MintsHelper::ao_erf_eri, "AO ERF integrals", py::arg("omega"), py::arg("factory") = nullptr)
        .def("ao_f12", normal_f12(&MintsHelper::ao_f12), "AO F12 integrals", py::arg("corr"))
        .def("ao_f12", normal_f122(&MintsHelper::ao_f12), "AO F12 integrals", py::arg("corr"), py::arg("bs1"),
             py::arg("bs2"), py::arg("bs3"), py::arg("bs4"))
        .def("ao_f12_scaled", normal_f12(&MintsHelper::ao_f12_scaled), "AO F12 intgerals", py::arg("corr"))
        .def("ao_f12_scaled", normal_f122(&MintsHelper::ao_f12_scaled), "AO F12 intgerals", py::arg("corr"),
             py::arg("bs1"), py::arg("bs2"), py::arg("bs3"), py::arg("bs4"))
        .def("ao_f12_squared", normal_f12(&MintsHelper::ao_f12_squared), "AO F12 squared integrals", py::arg("corr"))
        .def("ao_f12_squared", normal_f122(&MintsHelper::ao_f12_squared), "AO F12 squared integrals", py::arg("corr"),
             py::arg("bs1"), py::arg("bs2"), py::arg("bs3"), py::arg("bs4"))
        .def("ao_f12g12", &MintsHelper::ao_f12g12, "AO F12G12 integrals", py::arg("corr"))
        .def("ao_f12_double_commutator", &MintsHelper::ao_f12_double_commutator, "AO F12 double commutator integrals",
             py::arg("corr"))
        .def("ao_3coverlap", normal_eri(&MintsHelper::ao_3coverlap), "3 Center overlap integrals")
        .def("ao_3coverlap", normal_3c(&MintsHelper::ao_3coverlap), "3 Center overalp integrals", py::arg("bs1"),
             py::arg("bs2"), py::arg("bs3"))

        // Two-electron MO and transformers
        .def("mo_eri", eri(&MintsHelper::mo_eri), "MO ERI Integrals. Pass appropriate MO coefficients", py::arg("C1"),
             py::arg("C2"), py::arg("C3"), py::arg("C4"))
        .def("mo_erf_eri", erf(&MintsHelper::mo_erf_eri), "MO ERFC Omega Integrals", py::arg("omega"), py::arg("C1"),
             py::arg("C2"), py::arg("C3"), py::arg("C4"))
        .def("mo_f12", &MintsHelper::mo_f12, "MO F12 Integrals", py::arg("corr"), py::arg("C1"), py::arg("C2"),
             py::arg("C3"), py::arg("C4"))
        .def("mo_f12_squared", &MintsHelper::mo_f12_squared, "MO F12 squared integrals", py::arg("corr"), py::arg("C1"),
             py::arg("C2"), py::arg("C3"), py::arg("C4"))
        .def("mo_f12g12", &MintsHelper::mo_f12g12, "MO F12G12 integrals", py::arg("corr"), py::arg("C1"), py::arg("C2"),
             py::arg("C3"), py::arg("C4"))
        .def("mo_f12_double_commutator", &MintsHelper::mo_f12_double_commutator, "MO F12 double commutator integrals",
             py::arg("corr"), py::arg("C1"), py::arg("C2"), py::arg("C3"), py::arg("C4"))
        .def("mo_spin_eri", &MintsHelper::mo_spin_eri, "Symmetric MO Spin ERI Integrals", py::arg("C1"), py::arg("C2"))
        .def("mo_transform", &MintsHelper::mo_transform, "N^5 ao to mo transfrom, in memory", py::arg("Iso"),
             py::arg("C1"), py::arg("C2"), py::arg("C3"), py::arg("C4"))
        .def("set_rel_basisset", &MintsHelper::set_rel_basisset, "Sets the relativistic basis set",
             py::arg("rel_basis"))
        .def("play", &MintsHelper::play, "play function")

        // Contracted gradient terms
        .def("overlap_grad", &MintsHelper::overlap_grad, "First nuclear derivative overlap integrals")
        .def("kinetic_grad", &MintsHelper::kinetic_grad, "First nuclear derivative kinetic integrals")
        .def("potential_grad", &MintsHelper::potential_grad, "First nuclear derivative potential integrals")
        .def("perturb_grad", perturb_grad_options(&MintsHelper::perturb_grad),
             "First nuclear derivative perturb integrals")
        .def("perturb_grad", perturb_grad_xyz(&MintsHelper::perturb_grad), "First nuclear derivative perturb integrals")
        .def("core_hamiltonian_grad", &MintsHelper::core_hamiltonian_grad,
             "First nuclear derivative T + V + Perturb integrals")

        // First and second derivatives of one and two electron integrals in AO and MO basis.
        .def("ao_oei_deriv1", &MintsHelper::ao_oei_deriv1,
             "Gradient of AO basis OEI integrals: returns (3 * natoms) matrices")
        .def("ao_oei_deriv2", &MintsHelper::ao_oei_deriv2,
             "Hessian  of AO basis OEI integrals: returns (3 * natoms)^2 matrices")
        .def("ao_tei_deriv1", &MintsHelper::ao_tei_deriv1,
             "Gradient of AO basis TEI integrals: returns (3 * natoms) matrices", py::arg("atom"),
             py::arg("omega") = 0.0, py::arg("factory") = nullptr)
        .def("ao_tei_deriv2", &MintsHelper::ao_tei_deriv2,
             "Hessian  of AO basis TEI integrals: returns (3 * natoms)^2 matrices")
        .def("mo_oei_deriv1", &MintsHelper::mo_oei_deriv1,
             "Gradient of MO basis OEI integrals: returns (3 * natoms) matrices")
        .def("mo_oei_deriv2", &MintsHelper::mo_oei_deriv2,
             "Hessian  of MO basis OEI integrals: returns (3 * natoms)^2 matrices")
        .def("mo_tei_deriv1", &MintsHelper::mo_tei_deriv1,
             "Gradient of MO basis TEI integrals: returns (3 * natoms) matrices")
        .def("mo_tei_deriv2", &MintsHelper::mo_tei_deriv2,
             "Hessian  of MO basis TEI integrals: returns (3 * natoms)^2 matrices");

    py::class_<Vector3>(m, "Vector3",
                        "Class for vectors of length three, often Cartesian coordinate vectors, "
                        "and their common operations")
        .def(py::init<double>())
        .def(py::init<double, double, double>())
        .def(py::init<const Vector3&>())

        //      def(self = other<double>()).
        .def(py::self += py::self)
        .def(py::self -= py::self)
        .def(py::self *= double())

        // def(py::self + py::self).
        // def(py::self - py::self).
        // def(-py::self).
        .def("dot", &Vector3::dot, "Returns dot product of arg1 and arg2")
        .def("distance", &Vector3::distance, "Returns distance between two points represented by arg1 and arg2")
        .def("normalize", &Vector3::normalize, "Returns vector of unit length and arg1 direction")
        .def("norm", &Vector3::norm, "Returns Euclidean norm of arg1")
        .def("cross", &Vector3::cross, "Returns cross product of arg1 and arg2")
        .def("__str__", &Vector3::to_string, "Returns a string representation of arg1, suitable for printing.")
        .def(float() * py::self) /* The __mult__ operator */
        .def("__getitem__", &Vector3::get, "Returns the arg2-th element of arg1.")
        .def("__setitem__", &Vector3::set, "Set the ar2-th element of to val.");

    typedef void (SymmetryOperation::*intFunction)(int);
    typedef void (SymmetryOperation::*doubleFunction)(double);

    py::class_<SymmetryOperation>(m, "SymmetryOperation",
                                  "Class to provide a 3 by 3 matrix representation of a symmetry "
                                  "operation, such as a rotation or reflection.")
        .def(py::init<const SymmetryOperation&>())
        .def("trace", &SymmetryOperation::trace, "Returns trace of transformation matrix")
        .def("matrix", &SymmetryOperation::matrix, "Return the matrix for the operation on Cartesians")
        .def("zero", &SymmetryOperation::zero, "Zero out the symmetry operation")
        .def("operate", &SymmetryOperation::operate, "Performs the operation arg2 * arg1")
        .def("transform", &SymmetryOperation::transform, "Performs the transform arg2 * arg1 * arg2~")
        .def("unit", &SymmetryOperation::unit, "Set equal to a unit matrix")
        .def("E", &SymmetryOperation::E, "Set equal to E")
        .def("i", &SymmetryOperation::i, "Set equal to an inversion")
        .def("sigma_xy", &SymmetryOperation::sigma_xy, "Set equal to reflection in xy plane")
        .def("sigma_yz", &SymmetryOperation::sigma_yz, "Set equal to reflection in yz plane")
        .def("sigma_xz", &SymmetryOperation::sigma_xz, "Set equal to reflection in xz plane")
        .def("rotate_n", intFunction(&SymmetryOperation::rotation), "Set equal to a clockwise rotation by 2pi/n")
        .def("rotate_theta", doubleFunction(&SymmetryOperation::rotation), "Set equal to a clockwise rotation by theta")
        .def("c2_x", &SymmetryOperation::c2_x, "Set equal to C2 about the x axis")
        .def("c2_y", &SymmetryOperation::c2_y, "Set equal to C2 about the y axis")
        .def("c2_z", &SymmetryOperation::c2_z, "Set equal to C2 about the z axis")
        .def("transpose", &SymmetryOperation::transpose, "Performs transposition of matrix operation")
        // Warning! The __getitem__ function below copies the elements, it does not refer to them.
        .def("__getitem__",
             [](const SymmetryOperation& self, size_t i) { return std::vector<double>(self[i], self[i] + 3); });

    py::class_<OrbitalSpace>(m, "OrbitalSpace", "Contains information about the orbitals")
        .def(py::init<const std::string&, const std::string&, const SharedMatrix&, const SharedVector&,
                      const std::shared_ptr<BasisSet>&, const std::shared_ptr<IntegralFactory>&>())
        .def(py::init<const std::string&, const std::string&, const SharedMatrix&, const std::shared_ptr<BasisSet>&,
                      const std::shared_ptr<IntegralFactory>&>())
        .def(py::init<const std::string&, const std::string&, const std::shared_ptr<Wavefunction>&>())
        .def("nirrep", &OrbitalSpace::nirrep, "Returns number of irreps")
        .def("id", &OrbitalSpace::id, "Unique identifier")
        .def("name", &OrbitalSpace::name, "Name of the orbital space")
        .def("C", &OrbitalSpace::C, "MO coefficient matrix, AO->MO or SO->MO transformation matrix")
        .def("evals", &OrbitalSpace::evals, "Corresponding eigenvalues of the C matrix")
        .def("basisset", &OrbitalSpace::basisset, "The AO basis set used to create C")
        .def("integral", &OrbitalSpace::integral, "The integral factory used to create C")
        .def("dim", &OrbitalSpace::dim, "MO dimensions")
        .def("print_out", &OrbitalSpace::print, "Print information about the orbital space to the output file")
        .def_static("build_cabs_space", &OrbitalSpace::build_cabs_space,
                    "Given two spaces, it projects out one space from the other and returns the new spaces \
                    The first argument (orb_space) is the space to project out. The returned space will be orthogonal to this \
                    The second argument (ri_space) is the space that is being projected on. The returned space = this space - orb_space \
                    The third argument is the tolerance for linear dependencies",
                    py::arg("orb_space"), py::arg("ri_space"), py::arg("linear_tol"))
        .def_static("build_ri_space", &OrbitalSpace::build_ri_space,
                    "Given two basis sets, it merges the basis sets and then constructs an orthogonalized \
                    space with the same span. Linearly dependent orbitals are thrown out. \
                    The first argument, molecule, is the molecule to construct the basis for \
                    The second argument, obs_key, is the option keyword for orbital basis set 'BASIS' \
                    The third argument, aux_key, is the option keyword for auxiliery basis set 'DF_BASIS_MP2' \
                    The fourth argument, lindep_tol, is the tolerance for linear dependencies",
                    py::arg("molecule"), py::arg("obs_key"), py::arg("aux_key"), py::arg("lindep_tol"));

    py::class_<PointGroup, std::shared_ptr<PointGroup>>(m, "PointGroup", "Contains information about the point group")
        .def(py::init<const std::string&>())
        .def("symbol", &PointGroup::symbol, "Returns Schoenflies symbol for point group")
        .def("order", &PointGroup::order, "Return the order of the point group")
        .def("char_table", &PointGroup::char_table, "Return the CharacterTable of the point group");
    // def("origin", &PointGroup::origin).
    //            def("set_symbol", &PointGroup::set_symbol);

    py::class_<CharacterTable, std::shared_ptr<CharacterTable>>(m, "CharacterTable",
                                                                "Contains the character table of the point group")
        .def(py::init<const unsigned char>())
        .def(py::init<const std::string&>())
        .def("gamma", &CharacterTable::gamma, "Returns the irrep with the given index in the character table")
        .def("order", &CharacterTable::order, "Return the order of the point group")
        .def("symm_operation", &CharacterTable::symm_operation, "Return the i'th symmetry operation. 0-indexed.");

    py::class_<IrreducibleRepresentation, std::shared_ptr<IrreducibleRepresentation>>(
        m, "IrreducibleRepresentation", "An irreducible representation of the point group")
        .def("character", &IrreducibleRepresentation::character,
             "Return the character of the i'th symmetry operation for the irrep. 0-indexed.")
        .def("symbol", &IrreducibleRepresentation::symbol, "Return the symbol for the irrep");

    typedef void (Molecule::*matrix_set_geometry)(const Matrix&);
    typedef Vector3 (Molecule::*nuclear_dipole1)(const Vector3&) const;
    typedef Vector3 (Molecule::*nuclear_dipole2)() const;
    typedef Vector3 (Molecule::*vector_by_index)(int) const;

    py::class_<Molecule, std::shared_ptr<Molecule>>(m, "Molecule", py::dynamic_attr(),
                                                    "Class to store the elements, coordinates, "
                                                    "fragmentation pattern, basis sets, charge, "
                                                    "multiplicity, etc. of a molecule.")
        .def("set_geometry", matrix_set_geometry(&Molecule::set_geometry),
             "Sets the geometry, given a (Natom X 3) matrix arg0 of coordinates [a0] (excluding dummies)")
        .def("nuclear_dipole", nuclear_dipole1(&Molecule::nuclear_dipole),
             "Gets the nuclear contribution to the dipole, with respect to a specified origin atg0")
        .def("nuclear_dipole", nuclear_dipole2(&Molecule::nuclear_dipole),
             "Gets the nuclear contribution to the dipole, with respect to the origin")
        .def("set_name", &Molecule::set_name, "Sets molecule name")
        .def("name", &Molecule::name, "Gets molecule name")
        .def("reinterpret_coordentry", &Molecule::set_reinterpret_coordentry,
             "Do reinterpret coordinate entries during update_geometry().")
        .def("fix_orientation", &Molecule::set_orientation_fixed, "Fix the orientation at its current frame")
        .def("fix_com", &Molecule::set_com_fixed,
             "Sets whether to fix the Cartesian position, or to translate to the C.O.M.")
        .def("orientation_fixed", &Molecule::orientation_fixed, "Get whether or not orientation is fixed")
        .def("com_fixed", &Molecule::com_fixed, "Gets whether or not center of mass is fixed")
        .def("symmetry_from_input", &Molecule::symmetry_from_input, "Returns the symmetry specified in the input")
        .def("add_atom", &Molecule::add_atom,
             "Adds to self Molecule an atom with atomic number *Z*, Cartesian coordinates in Bohr "
             "(*x*, *y*, *z*), atomic *symbol*, *mass*, and *charge*, extended atomic *label*, and mass number *A*",
             py::arg("Z"), py::arg("x"), py::arg("y"), py::arg("z"), py::arg("symbol"), py::arg("mass"),
             py::arg("charge"), py::arg("label"), py::arg("A"))
        .def("natom", &Molecule::natom, "Number of real atoms")
        .def("nallatom", &Molecule::nallatom, "Number of real and dummy atoms")
        .def("multiplicity", &Molecule::multiplicity, "Gets the multiplicity (defined as 2Ms + 1)")
        .def("nfragments", &Molecule::nfragments, "Gets the number of fragments in the molecule")
        .def("set_nuclear_charge", &Molecule::set_nuclear_charge,
             "Set the nuclear charge of the given atom arg0 to the value arg1 (primarily for ECP).")
        .def("basis_on_atom", &Molecule::basis_on_atom, py::return_value_policy::copy,
             "Gets the label of the orbital basis set on a given atom arg0")
        .def("print_in_input_format", &Molecule::print_in_input_format,
             "Prints the molecule as Cartesian or ZMatrix entries, just as inputted.")
        .def("create_psi4_string_from_molecule", &Molecule::create_psi4_string_from_molecule,
             "Gets a string re-expressing in input format the current state of the molecule."
             "Contains Cartesian geometry info, fragmentation, charges and multiplicities, "
             "and any frame restriction.")
        .def("save_xyz_file", &Molecule::save_xyz_file, "Saves an XYZ file to arg0")
        .def("save_string_xyz_file", &Molecule::save_string_xyz_file, "Saves an XYZ file to arg2")
        .def("save_string_xyz", &Molecule::save_string_xyz, "Saves the string of an XYZ file to arg2")
        .def("Z", &Molecule::Z, py::return_value_policy::copy,
             "Nuclear charge of atom arg0 (0-indexed without dummies)")
        .def("mass_number", &Molecule::mass_number, py::return_value_policy::copy,
             "Mass number (A) of atom if known, else -1")
        .def("x", &Molecule::x, "x position [Bohr] of atom arg0 (0-indexed without dummies)")
        .def("y", &Molecule::y, "y position [Bohr] of atom arg0 (0-indexed without dummies)")
        .def("z", &Molecule::z, "z position [Bohr] of atom arg0 (0-indexed without dummies)")
        .def("xyz", vector_by_index(&Molecule::xyz), "Return the Vector3 for atom i (0-indexed without dummies)",
             py::arg("i"))
        .def("fZ", &Molecule::fZ, py::return_value_policy::copy,
             "Nuclear charge of atom arg1 (0-indexed including dummies)")
        .def("fx", &Molecule::fx, "x position of atom arg0 (0-indexed including dummies in Bohr)")
        .def("fy", &Molecule::fy, "y position of atom arg0 (0-indexed including dummies in Bohr)")
        .def("fz", &Molecule::fz, "z position of atom arg0 (0-indexed including dummies in Bohr)")
        .def("true_atomic_number", &Molecule::true_atomic_number,
             "Gets atomic number of "
             "*atom* from element (0-indexed without dummies)",
             py::arg("atom"))
        .def("ftrue_atomic_number", &Molecule::ftrue_atomic_number,
             "Gets atomic number of "
             "*atom* from element (0-indexed including dummies)",
             py::arg("atom"))
        .def("center_of_mass", &Molecule::center_of_mass,
             "Computes center of mass of molecule (does not translate molecule)")
        .def("translate", &Molecule::translate, "Translates molecule by arg0")
        .def("move_to_com", &Molecule::move_to_com, "Moves molecule to center of mass")
        .def("mass", &Molecule::mass, "Returns mass of *atom* (0-indexed)", py::arg("atom"))
        .def("set_mass", &Molecule::set_mass,
             "Sets mass of *atom* (0-indexed) to *mass* (good for isotopic substitutions)", py::arg("atom"),
             py::arg("mass"))
        .def("fmass", &Molecule::fmass, "Gets mass of *atom* (0-indexed including dummies)", py::arg("atom"))
        .def("symbol", &Molecule::symbol,
             "Gets the cleaned up label of *atom* (C2 => C, H4 = H) (0-indexed without dummies)", py::arg("atom"))
        .def("fsymbol", &Molecule::fsymbol,
             "Gets the cleaned up label of *atom* (C2 => C, H4 = H) (0-indexed including dummies)", py::arg("atom"))
        .def("label", &Molecule::label,
             "Gets the original label of the *atom* as given in the input file (C2, H4)"
             "(0-indexed without dummies)",
             py::arg("atom"))
        .def("flabel", &Molecule::flabel,
             "Gets the original label of the *atom* arg0 as given in the input file (C2, H4)"
             "(0-indexed including dummies)",
             py::arg("atom"))
        .def("charge", &Molecule::charge, "Gets charge of *atom* (0-indexed without dummies)", py::arg("atom"))
        .def("fcharge", &Molecule::fcharge, "Gets charge of *atom* (0-indexed including dummies)", py::arg("atom"))
        .def("molecular_charge", &Molecule::molecular_charge, "Gets the molecular charge")
        .def("extract_subsets", &Molecule::py_extract_subsets_1,
             "Returns copy of self with arg0 fragments Real and arg1 fragments Ghost")
        .def("extract_subsets", &Molecule::py_extract_subsets_2,
             "Returns copy of self with arg0 fragments Real and arg1 fragment Ghost")
        .def("extract_subsets", &Molecule::py_extract_subsets_3,
             "Returns copy of self with arg0 fragment Real and arg1 fragments Ghost")
        .def("extract_subsets", &Molecule::py_extract_subsets_4,
             "Returns copy of self with arg0 fragment Real and arg1 fragment Ghost")
        .def("extract_subsets", &Molecule::py_extract_subsets_5, "Returns copy of self with arg0 fragments Real")
        .def("extract_subsets", &Molecule::py_extract_subsets_6, "Returns copy of self with arg0 fragment Real")
        .def("activate_all_fragments", &Molecule::activate_all_fragments,
             "Sets all fragments in the molecule to be active")
        .def("deactivate_all_fragments", &Molecule::deactivate_all_fragments,
             "Sets all fragments in the molecule to be inactive")
        .def("set_active_fragments", &Molecule::set_active_fragments,
             "Sets the specified list arg0 of fragments to be Real")
        .def("set_active_fragment", &Molecule::set_active_fragment, "Sets the specified fragment arg0 to be Real")
        .def("set_ghost_fragments", &Molecule::set_ghost_fragments,
             "Sets the specified list arg0 of fragments to be Ghost")
        .def("set_ghost_fragment", &Molecule::set_ghost_fragment, "Sets the specified fragment arg0 to be Ghost")
        .def("get_fragments", &Molecule::get_fragments,
             "Returns list of pairs of atom ranges defining each fragment from parent molecule"
             "(fragments[frag_ind] = <Afirst,Alast+1>)")
        .def("get_fragment_types",
             [](Molecule& mol) {
                 const std::string FragmentTypeList[] = {"Absent", "Real", "Ghost"};
                 std::vector<std::string> srt;
                 for (auto item : mol.get_fragment_types()) srt.push_back(FragmentTypeList[item]);
                 return srt;
             },
             "Returns a list describing how to handle each fragment {Real, Ghost, Absent}")
        .def("get_fragment_charges", &Molecule::get_fragment_charges, "Gets the charge of each fragment")
        .def("get_fragment_multiplicities", &Molecule::get_fragment_multiplicities,
             "Gets the multiplicity of each fragment")
        .def("atom_at_position", &Molecule::atom_at_position1,
             "Tests to see if an atom is at the position *coord* with a given tolerance *tol*", py::arg("coord"),
             py::arg("tol"))
        .def("atom_at_position", &Molecule::atom_at_position3,
             "Tests to see if an atom is at the position *coord* with a given tolerance *tol*", py::arg("coord"),
             py::arg("tol"))
        .def("print_out", &Molecule::print, "Prints the molecule in Cartesians in input units to output file")
        .def("print_out_in_bohr", &Molecule::print_in_bohr, "Prints the molecule in Cartesians in Bohr to output file")
        .def("print_out_in_angstrom", &Molecule::print_in_angstrom,
             "Prints the molecule in Cartesians in Angstroms to output file")
        .def("print_cluster", &Molecule::print_cluster,
             "Prints the molecule in Cartesians in input units adding fragment separators")
        .def("rotational_constants", [](Molecule& mol) { return mol.rotational_constants(1.0e-8); },
             "Returns the rotational constants [cm^-1] of the molecule")
        .def("print_rotational_constants", &Molecule::print_rotational_constants,
             "Print the rotational constants to output file")
        .def("nuclear_repulsion_energy", &Molecule::nuclear_repulsion_energy,
             py::arg("dipole_field") = std::vector<double>(3, 0.0), "Computes nuclear repulsion energy")
        .def("nuclear_repulsion_energy_deriv1", &Molecule::nuclear_repulsion_energy_deriv1,
             py::arg("dipole_field") = std::vector<double>(3, 0.0),
             "Returns first derivative of nuclear repulsion energy as a matrix (natom, 3)")
        .def("nuclear_repulsion_energy_deriv2", &Molecule::nuclear_repulsion_energy_deriv2,
             "Returns second derivative of nuclear repulsion energy as a matrix (natom X 3, natom X 3)")
        .def("find_point_group", &Molecule::find_point_group,
             "Finds computational molecular point group, user can override this with the symmetry "
             "keyword")
        .def("reset_point_group", &Molecule::reset_point_group, "Overrides symmetry from outside the molecule string")
        .def("set_point_group", &Molecule::set_point_group,
             "Sets the molecular point group to the point group object arg0")
        .def("get_full_point_group", &Molecule::full_point_group, "Gets point group name such as C3v or S8")
        .def("get_full_point_group_with_n", &Molecule::full_point_group_with_n,
             "Gets point group name such as Cnv or Sn")
        .def("full_pg_n", &Molecule::full_pg_n,
             "Gets n in Cnv, etc.; If there is no n (e.g. Td) it's the highest-order rotation axis")
        .def("point_group", &Molecule::point_group, "Returns the current point group object")
        .def("schoenflies_symbol", &Molecule::schoenflies_symbol, "Returns the Schoenflies symbol")
        .def("form_symmetry_information", &Molecule::form_symmetry_information,
             "Uses the point group object obtain by calling point_group()")
        .def("symmetrize", &Molecule::symmetrize_to_abelian_group,
             "Finds the highest point Abelian point group within the specified tolerance, and "
             "forces the geometry to have that symmetry.")
        .def("inertia_tensor", &Molecule::inertia_tensor, "Returns intertial tensor")
        .def("is_variable", &Molecule::is_variable, "Checks if variable arg0 is in the structural variables list")
        .def("set_variable", &Molecule::set_variable,
             "Sets the value arg1 to the variable arg0 in the list of structure variables, then "
             "calls update_geometry()")
        .def("get_variable", &Molecule::get_variable,
             "Returns the value of variable arg0 in the structural variables list")
        .def("update_geometry", &Molecule::update_geometry,
             "Reevaluates the geometry with current variable values, orientation directives, etc. "
             "by clearing the atoms list and rebuilding it. Idempotent. Use liberally."
             "Must be called after initial Molecule definition by string.")
        .def("set_molecular_charge", &Molecule::set_molecular_charge,
             "Change the overall molecular charge. Setting in initial molecule string or constructor preferred.")
        .def(
            "set_multiplicity", &Molecule::set_multiplicity,
            "Change the multiplicity (defined as 2S + 1). Setting in initial molecule string or constructor preferred.")
        .def("set_basis_all_atoms", &Molecule::set_basis_all_atoms, "Sets basis set arg0 to all atoms")
        .def("set_basis_by_symbol", &Molecule::set_basis_by_symbol,
             "Sets basis set arg1 to all atoms with symbol (e.g., H) arg0")
        .def("set_basis_by_label", &Molecule::set_basis_by_label,
             "Sets basis set arg1 to all atoms with label (e.g., H4) arg0")
        .def("distance_matrix", &Molecule::distance_matrix, "Returns Matrix of interatom distances")
        .def("print_distances", &Molecule::print_distances, "Print the interatomic distance geometrical parameters")
        .def("print_bond_angles", &Molecule::print_bond_angles, "Print the bond angle geometrical parameters")
        .def("print_out_of_planes", &Molecule::print_out_of_planes,
             "Print the out-of-plane angle geometrical parameters to output file")
        .def("irrep_labels", &Molecule::irrep_labels, "Returns Irreducible Representation symmetry labels")
        .def("units",
             [](Molecule& mol) {
                 const std::string GeometryUnitsList[] = {"Angstrom", "Bohr"};
                 std::string srt = GeometryUnitsList[mol.units()];
                 return srt;
             },
             "Returns units used to define the geometry, i.e. 'Angstrom' or 'Bohr'")
        .def("set_units", &Molecule::set_units,
             "Sets units (Angstrom or Bohr) used to define the geometry. Imposes Psi4 physical constants conversion "
             "for input_units_to_au.")
        .def("input_units_to_au", &Molecule::input_units_to_au, "Returns unit conversion to [a0] for geometry")
        .def("set_input_units_to_au", &Molecule::set_input_units_to_au, "Sets unit conversion to [a0] for geometry")
        .def("clone", &Molecule::clone, "Returns a new Molecule identical to arg1")
        .def_static("from_dict", from_dict,
                    "Returns a new Molecule constructed from python dictionary. In progress: name and capabilities "
                    "should not be relied upon")
        .def("rotational_symmetry_number", &Molecule::rotational_symmetry_number,
             "Returns number of unique orientations of the rigid molecule that only interchange identical atoms")
        .def("rotor_type",
             [](Molecule& mol) {
                 const std::string RotorTypeList[] = {"RT_ASYMMETRIC_TOP", "RT_SYMMETRIC_TOP", "RT_SPHERICAL_TOP",
                                                      "RT_LINEAR", "RT_ATOM"};
                 std::string srt = RotorTypeList[mol.rotor_type()];
                 return srt;
             },
             "Returns rotor type, e.g. 'RT_ATOM' or 'RT_SYMMETRIC_TOP'")
        .def("geometry", &Molecule::geometry,
             "Gets the geometry [Bohr] as a (Natom X 3) matrix of coordinates (excluding dummies)")
        .def("full_geometry", &Molecule::full_geometry,
             "Gets the geometry [Bohr] as a (Natom X 3) matrix of coordinates (including dummies)")
        .def("set_full_geometry", matrix_set_geometry(&Molecule::set_full_geometry),
             "Sets the geometry, given a (Natom X 3) matrix arg0 of coordinates (in Bohr) (including dummies");

    py::class_<PetiteList, std::shared_ptr<PetiteList>>(m, "PetiteList", "Handles symmetry transformations")
        .def("aotoso", &PetiteList::aotoso, "Return the AO->SO coefficient matrix")
        .def("sotoao", &PetiteList::sotoao, "Return the SO->AO coefficient matrix")
        .def("print", &PetiteList::print, "Print to outfile");

    typedef void (BasisSet::*basis_print_out)() const;
    typedef const GaussianShell& (BasisSet::*no_center_version)(int)const;
    typedef const GaussianShell& (BasisSet::*center_version)(int, int)const;
    typedef std::shared_ptr<BasisSet> (BasisSet::*ptrversion)(const std::shared_ptr<BasisSet>&) const;
    typedef int (BasisSet::*ncore_no_args)() const;
    typedef int (BasisSet::*ncore_one_arg)(const std::string&) const;

    py::class_<BasisSet, std::shared_ptr<BasisSet>>(m, "BasisSet", "Contains basis set information")
        .def(py::init<const std::string&, std::shared_ptr<Molecule>,
                      std::map<std::string, std::map<std::string, std::vector<ShellInfo>>>&,
                      std::map<std::string, std::map<std::string, std::vector<ShellInfo>>>&>())
        .def("name", &BasisSet::name, "Callback handle, may represent string or function")
        .def("blend", &BasisSet::target, "Plus-separated string of [basisname] values")
        .def("molecule", &BasisSet::molecule, "Molecule object")
        .def("print_out", basis_print_out(&BasisSet::print), "Prints basis set info to outfile")
        .def("print_detail_out", basis_print_out(&BasisSet::print_detail), "Prints detailed basis set info to outfile")
        .def("genbas", &BasisSet::print_detail_cfour, "Returns basis set per atom in CFOUR format")
        .def_static("make_filename", &BasisSet::make_filename,
                    "Returns filename for basis name: pluses, stars, parentheses replaced and gbs "
                    "extension added")
        .def_static("zero_ao_basis_set", &BasisSet::zero_ao_basis_set,
                    "Returns a BasisSet object that actually has a single s-function at the origin "
                    "with an exponent of 0.0 and contraction of 1.0.")
        .def("nbf", &BasisSet::nbf,
             "Returns number of basis functions (Cartesian or spherical depending on has_puream)")
        .def("nao", &BasisSet::nao, "Returns number of atomic orbitals (Cartesian)")
        .def("nprimitive", &BasisSet::nprimitive, "Returns total number of primitives in all contractions")
        .def("nshell", &BasisSet::nshell, "Returns number of shells")
        .def("shell", no_center_version(&BasisSet::shell), py::return_value_policy::copy,
             "Return the si'th Gaussian shell", py::arg("si"))
        .def("shell", center_version(&BasisSet::shell), py::return_value_policy::copy,
             "Return the si'th Gaussian shell on center", py::arg("center"), py::arg("si"))
        .def("n_frozen_core", &BasisSet::n_frozen_core,
             "Returns the number of orbital (non-ECP) frozen core electrons. For a given molecule and "
             "|globals__freeze_core|, `(n_ecp_core()/2 + n_frozen_core()) = constant`.")
        .def("n_ecp_core", ncore_no_args(&BasisSet::n_ecp_core),
             "Returns the total number of core electrons associated with all ECPs in this basis.")
        .def("n_ecp_core", ncore_one_arg(&BasisSet::n_ecp_core),
             "Returns the number of core electrons associated with any ECP on the specified atom type for this basis "
             "set.")
        .def("has_ECP", &BasisSet::has_ECP, "Whether this basis set object has an ECP associated with it.")
        .def("max_am", &BasisSet::max_am, "Returns maximum angular momentum used")
        .def("has_puream", &BasisSet::has_puream, "Spherical harmonics?")
        .def("shell_to_basis_function", &BasisSet::shell_to_basis_function,
             "Given a shell return its first basis function", py::arg("i"))
        .def("shell_to_center", &BasisSet::shell_to_center, "Return the atomic center for the i'th shell",
             py::arg("i"))  // are shell_to_basis_function and shell_to_ao_function the same?
        .def("shell_to_ao_function", &BasisSet::shell_to_ao_function,
             "Return the function number for the first function for the i'th shell", py::arg("i"))
        .def("function_to_shell", &BasisSet::function_to_shell,
             "Given a function number what shell does it correspond to", py::arg("i"))
        .def("function_to_center", &BasisSet::function_to_center, "The atomic center for the i'th function",
             py::arg("i"))
        .def("nshell_on_center", &BasisSet::nshell_on_center, "Return the number of shells on a given center",
             py::arg("i"))
        //            def("decontract", &BasisSet::decontract, "docstring").
        .def("ao_to_shell", &BasisSet::ao_to_shell,
             "Given a cartesian function (AO) number what shell does it correspond to", py::arg("i"))
        .def("move_atom", &BasisSet::move_atom,
             "Translate a given atom by a given amount.  Does not affect the underlying molecule object.")
        .def("max_function_per_shell", &BasisSet::max_function_per_shell,
             "The max number of basis functions in a shell")
        .def("max_nprimitive", &BasisSet::max_nprimitive, "The max number of primitives in a shell")
        .def_static("construct_from_pydict", &construct_basisset_from_pydict, "docstring");

    py::class_<SOBasisSet, std::shared_ptr<SOBasisSet>>(
        m, "SOBasisSet",
        "An SOBasis object describes the transformation from an atomic orbital basis to a symmetry orbital basis.")
        .def("petite_list", &SOBasisSet::petite_list, "Return the PetiteList object used in creating this SO basis");

    py::class_<ExternalPotential, std::shared_ptr<ExternalPotential>>(
        m, "ExternalPotential", "Stores external potential field, computes external potential matrix")
        .def(py::init<>())
        .def("setName", &ExternalPotential::setName, "Sets the name")
        .def("addCharge", &ExternalPotential::addCharge, "Add a charge Z at (x,y,z)", py::arg("Z"), py::arg("x"),
             py::arg("y"), py::arg("z"))
        .def("addBasis", &ExternalPotential::addBasis, "Add a basis of S auxiliary functions iwth Df coefficients",
             py::arg("basis"), py::arg("coefs"))
        .def("clear", &ExternalPotential::clear, "Reset the field to zero (eliminates all entries)")
        .def("computePotentialMatrix", &ExternalPotential::computePotentialMatrix,
             "Compute the external potential matrix in the given basis set", py::arg("basis"))
        .def("print_out", &ExternalPotential::py_print, "Print python print helper to the outfile");

    typedef std::shared_ptr<Localizer> (*localizer_with_type)(const std::string&, std::shared_ptr<BasisSet>,
                                                              std::shared_ptr<Matrix>);

    py::class_<Localizer, std::shared_ptr<Localizer>>(m, "Localizer",
                                                      "Class containing orbital localization procedures")
        .def_static("build", localizer_with_type(&Localizer::build), "Build the localization scheme")
        .def("localize", &Localizer::localize, "Perform the localization procedure")
        .def_property_readonly("L", py::cpp_function(&Localizer::L), "Localized orbital coefficients")
        .def_property_readonly("U", py::cpp_function(&Localizer::U), "Orbital rotation matrix")
        .def_property_readonly("converged", py::cpp_function(&Localizer::converged),
                               "Did the localization procedure converge?");

    py::class_<BoysLocalizer, std::shared_ptr<BoysLocalizer>, Localizer>(m, "BoysLocalizer",
                                                                         "Performs Boys orbital localization");
    py::class_<PMLocalizer, std::shared_ptr<PMLocalizer>, Localizer>(m, "PMLocalizer",
                                                                     "Performs Pipek-Mezey orbital localization");

    py::class_<FCHKWriter, std::shared_ptr<FCHKWriter>>(m, "FCHKWriter",
                                                        "Extracts information from a wavefunction object, \
                                                                          and writes it to an FCHK file")
        .def(py::init<std::shared_ptr<Wavefunction>>())
        .def("write", &FCHKWriter::write, "Write wavefunction information to file", py::arg("filename"));

    py::class_<MoldenWriter, std::shared_ptr<MoldenWriter>>(m, "MoldenWriter",
                                                            "Writes wavefunction information in molden format")
        .def(py::init<std::shared_ptr<Wavefunction>>())
        .def("write", &MoldenWriter::write, "Writes wavefunction information in molden format", py::arg("filename"),
             py::arg("Ca"), py::arg("Cb"), py::arg("Ea"), py::arg("Eb"), py::arg("OccA"), py::arg("OccB"),
             py::arg("dovirtual"));

    py::class_<NBOWriter, std::shared_ptr<NBOWriter>>(m, "NBOWriter", "The Natural Bond Orbital Writer")
        .def(py::init<std::shared_ptr<Wavefunction>>())
        .def("write", &NBOWriter::write, "Write the natural bond orbitals to a file", py::arg("filename"));

    py::class_<MOWriter, std::shared_ptr<MOWriter>>(m, "MOWriter", "Writes the MOs")
        .def(py::init<std::shared_ptr<Wavefunction>>())
        .def("write", &MOWriter::write, "Write the MOs");  // should the writer.h file take a filename as an argument?

    py::class_<OperatorSymmetry, std::shared_ptr<OperatorSymmetry>>(m, "MultipoleSymmetry", "docstring")
        .def(py::init<int, const std::shared_ptr<Molecule>&, const std::shared_ptr<IntegralFactory>&,
                      const std::shared_ptr<MatrixFactory>&>())
        .def("create_matrices", &OperatorSymmetry::create_matrices, "docstring");

    py::class_<CorrelationFactor, std::shared_ptr<CorrelationFactor>>(m, "CorrelationFactor", "docstring")
        .def(py::init<size_t>())
        .def(py::init<std::shared_ptr<Vector>, std::shared_ptr<Vector>>())
        .def("set_params", &CorrelationFactor::set_params, "Set coefficient and exponent", py::arg("coeff"),
             py::arg("exponent"));

    py::class_<FittedSlaterCorrelationFactor, std::shared_ptr<FittedSlaterCorrelationFactor>, CorrelationFactor>(
        m, "FittedSlaterCorrelationFactor", "docstring")
        .def(py::init<double>())
        .def("exponent", &FittedSlaterCorrelationFactor::exponent);

    py::class_<CorrelationTable, std::shared_ptr<CorrelationTable>>(
        m, "CorrelationTable", "Provides a correlation table between two point groups")
        .def(py::init<std::shared_ptr<PointGroup>, std::shared_ptr<PointGroup>>())
        .def("group", &CorrelationTable::group, "Returns higher order point group")
        .def("subgroup", &CorrelationTable::subgroup, "Returns lower order pointgroup")
        .def("n", &CorrelationTable::n, "Returns the number of irreps in high order group")
        .def("subn", &CorrelationTable::subn, "Returns number of irreps in subgroup")
        .def("degen", &CorrelationTable::degen, "Returns the degenercy of the irrep")
        .def("subdegen", &CorrelationTable::subdegen, "Returns the degeneracy of the subgroup irrep")
        .def("ngamma", &CorrelationTable::ngamma,
             "Returns the number of irreps in the low order group that an irrep \
             from the high order group can be reduced to.")
        .def("group", &CorrelationTable::gamma, "Returns the higher order point group");
}
