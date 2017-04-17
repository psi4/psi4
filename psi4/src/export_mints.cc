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

#include "psi4/libmints/deriv.h"
#include "psi4/libmints/twobody.h"
#include "psi4/libmints/integralparameters.h"
#include "psi4/libmints/orbitalspace.h"
#include "psi4/libmints/view.h"
#include "psi4/libmints/local.h"
#include "psi4/libmints/vector3.h"
#include "psi4/libmints/pointgrp.h"
#include "psi4/libmints/extern.h"
#include "psi4/libmints/sobasis.h"
#include "psi4/libmints/basisset_parser.h"
#include "psi4/libmints/petitelist.h"
#include "psi4/libmints/integral.h"
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


void export_mints(py::module& m)
{
    // This is needed to wrap an STL vector into Boost.Python. Since the vector
    // is going to contain std::shared_ptr's we MUST set the no_proxy flag to true
    // (as it is) to tell Boost.Python to not create a proxy class to handle
    // the vector's data type.
    py::bind_vector<std::vector<std::shared_ptr<Matrix>>>(m, "VectorMatrix");

    // Other vector types
    //py::class_<std::vector<double> >(m, "vector_of_doubles", "docstring").
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

    py::class_<Dimension>(m, "Dimension", "docstring").
            def(py::init<int>()).
            def(py::init<int, const std::string&>()).
            def("print_out",
                &Dimension::print,
                "docstring").
            def("init",
                &Dimension::init,
                "Re-initializes the dimension object").
            def("n", &Dimension::n,
                //py::return_value_policy<copy_const_reference>(),
                py::return_value_policy::copy,
                "The order of the dimension").
            def_property("name",
                         //make_function(&Dimension::name, return_value_policy<copy_const_reference>()),
                         py::cpp_function(&Dimension::name),
                         py::cpp_function(&Dimension::set_name),
                         "The name of the dimension. Used in printing.").
            //def("__getitem__", &Dimension::get, return_value_policy<copy_const_reference>(), "docstring").
            def("__getitem__", &Dimension::get, py::return_value_policy::copy, "docstring").
            def("__setitem__", &Dimension::set, "docstring");

    py::class_<Vector, std::shared_ptr<Vector> >(m, "Vector", "docstring", py::dynamic_attr()).
            def(py::init<int>()).
            def(py::init<const Dimension&>()).
            def(py::init<const std::string&, int>()).
            def(py::init<const std::string&, const Dimension&>()).
            def_property("name",
                         py::cpp_function(&Vector::name),
                         py::cpp_function(&Vector::set_name),
                         "The name of the Vector. Used in printing.").
            def("get", vector_getitem_1(&Vector::get), "docstring").
            def("get", vector_getitem_2(&Vector::get), "docstring").
            def("set", vector_setitem_1(&Vector::set), "docstring").
            def("set", vector_setitem_2(&Vector::set), "docstring").
            def("print_out", &Vector::print_out, "docstring").
            def("scale", &Vector::scale, "docstring").
            def("dim", &Vector::dim, "docstring").
            def("__getitem__", vector_getitem_1(&Vector::pyget), "docstring").
            def("__setitem__", vector_setitem_1(&Vector::pyset), "docstring").
            def("__getitem__", vector_getitem_n(&Vector::pyget), "docstring").
            def("__setitem__", vector_setitem_n(&Vector::pyset), "docstring").
            def("nirrep", &Vector::nirrep, "docstring").
            def("array_interface", [](Vector &v){
                // Dont ask, hopefully pybind11 will work on this
                py::list ret;
                std::vector<py::buffer_info> buff_vec(v.array_interface());
                std::string typestr = "<";
                {
                   std::stringstream sstr;
                   sstr << (int)sizeof(double);
                   typestr += "f" + sstr.str();
                }

                for (auto const& buff: v.array_interface()){
                    py::dict interface;
                    interface["typestr"] = py::str(typestr);
                    interface["data"] = py::make_tuple((long)buff.ptr, false);
                    py::list shape;
                    for (auto const& s: buff.shape){
                        shape.append(py::int_(s));
                    }
                    interface["shape"] = shape;
                    ret.append(interface);

                }
                return ret;
            });

    typedef void  (IntVector::*int_vector_set)(int, int, int);
    py::class_<IntVector, std::shared_ptr<IntVector> >(m, "IntVector", "docstring").
            def(py::init<int>()).
            def("get", &IntVector::get, "docstring").
            def("set", int_vector_set(&IntVector::set), "docstring").
            def("print_out", &IntVector::print_out, "docstring").
            def("dim", &IntVector::dim, "docstring").
            def("nirrep", &IntVector::nirrep, "docstring");

    py::enum_<diagonalize_order>(m, "DiagonalizeOrder", "docstring")
            .value("Ascending", ascending)
            .value("Descending", descending)
            .export_values();

    py::enum_<Molecule::GeometryUnits>(m, "GeometryUnits", "docstring")
            .value("Angstrom", Molecule::Angstrom)
            .value("Bohr", Molecule::Bohr)
            .export_values();

    typedef void   (Matrix::*matrix_multiply)(bool, bool, double, const SharedMatrix&, const SharedMatrix&, double);
    typedef void   (Matrix::*matrix_diagonalize)(SharedMatrix&, std::shared_ptr<Vector>&, diagonalize_order);
    typedef void   (Matrix::*matrix_one)(const SharedMatrix&);
    typedef double (Matrix::*double_matrix_one)(const SharedMatrix&);
    typedef void   (Matrix::*matrix_two)(const SharedMatrix&, const SharedMatrix&);
    typedef void   (Matrix::*matrix_save)(const std::string&, bool, bool, bool);
    typedef void   (Matrix::*matrix_set1)(double);
    typedef void   (Matrix::*matrix_set3)(int, int, double);
    typedef void   (Matrix::*matrix_set4)(int, int, int, double);
    typedef double (Matrix::*matrix_get3)(const int&, const int&, const int&) const;
    typedef double (Matrix::*matrix_get2)(const int&, const int&) const;
    typedef void   (Matrix::*matrix_load)(const std::string&);
    typedef const Dimension& (Matrix::*matrix_ret_dimension)() const;

    py::class_<Matrix, std::shared_ptr<Matrix>>(m, "Matrix", "docstring", py::dynamic_attr()).
            def(py::init<int, int>()).
            def(py::init<const std::string&, int, int>()).
            def(py::init<const std::string&, const Dimension&, const Dimension&>()).
            def(py::init<const std::string&>()).
            def("clone", &Matrix::clone, "docstring").
            def_property("name",
                         py::cpp_function(&Matrix::name),
                         py::cpp_function(&Matrix::set_name),
                         "The name of the Matrix. Used in printing.").
            // def("set_name", &Matrix::set_name, "docstring").
            // def("name", &Matrix::name, py::return_value_policy::copy, "docstring").
            def("print_out", &Matrix::print_out, "docstring").
            def("rows", &Matrix::rowdim, "docstring").
            def("cols", &Matrix::coldim, "docstring").
            def("rowdim", matrix_ret_dimension(&Matrix::rowspi), py::return_value_policy::copy, "docstring").
            def("coldim", matrix_ret_dimension(&Matrix::colspi), py::return_value_policy::copy, "docstring").
            def("nirrep", &Matrix::nirrep, py::return_value_policy::copy, "docstring").
            def("symmetry", &Matrix::symmetry, py::return_value_policy::copy, "docstring").
            def("identity", &Matrix::identity, "docstring").
            def("copy_lower_to_upper", &Matrix::copy_lower_to_upper, "docstring").
            def("copy_upper_to_lower", &Matrix::copy_upper_to_lower, "docstring").
            def("zero_lower", &Matrix::zero_lower, "docstring").
            def("zero_upper", &Matrix::zero_upper, "docstring").
            def("zero", &Matrix::zero, "docstring").
            def("zero_diagonal", &Matrix::zero_diagonal, "docstring").
            def("trace", &Matrix::trace, "docstring").
            //            def("transpose", &Matrix::transpose).
            def("add", matrix_one(&Matrix::add), "docstring").
            def("axpy", &Matrix::axpy, "docstring").
            def("subtract", matrix_one(&Matrix::subtract), "docstring").
            def("accumulate_product", matrix_two(&Matrix::accumulate_product), "docstring").
            def("scale", &Matrix::scale, "docstring").
            def("sum_of_squares", &Matrix::sum_of_squares, "docstring").
            def("add_and_orthogonalize_row", &Matrix::add_and_orthogonalize_row, "docstring").
            def("rms", &Matrix::rms, "docstring").
            def("absmax", &Matrix::absmax, "docstring").
            def("scale_row", &Matrix::scale_row, "docstring").
            def("scale_column", &Matrix::scale_column, "docstring").
            def("transform", matrix_one(&Matrix::transform), "docstring").
            def("transform", matrix_two(&Matrix::transform), "docstring").
            def("transform", matrix_one(&Matrix::back_transform), "docstring").
            def("back_transform", matrix_two(&Matrix::back_transform), "docstring").
            def("vector_dot", double_matrix_one(&Matrix::vector_dot), "docstring").
            def("gemm", matrix_multiply(&Matrix::gemm), "docstring").
            def("diagonalize", matrix_diagonalize(&Matrix::diagonalize), "docstring").
            def("cholesky_factorize", &Matrix::cholesky_factorize, "docstring").
            def("partial_cholesky_factorize", &Matrix::partial_cholesky_factorize, "docstring").
            //def("canonical_orthogonalization", &Matrix::canonical_orthogonalization, CanonicalOrthog()).
            // def("canonical_orthogonalization", &Matrix::canonical_orthogonalization, py::arg("delta") = 0.0, py::arg("eigvec") = SharedMatrix()).
            def("schmidt", &Matrix::schmidt).
            def("invert", &Matrix::invert, "docstring").
            def("apply_denominator", matrix_one(&Matrix::apply_denominator), "docstring").
            def("copy", matrix_one(&Matrix::copy), "docstring").
            def("power", &Matrix::power, "docstring").
            // def("doublet", &Matrix::doublet, py::arg("transA") = false, py::arg("transB") = false).
            // def("triplet", &Matrix::triplet, py::arg("transA") = false, py::arg("transB") = false,
            //                                  py::arg("transC") = false, "docstring").
            def("doublet", &Matrix::doublet, "docstring").
            def("triplet", &Matrix::triplet, "docstring").
            def("get", matrix_get3(&Matrix::get), "docstring").
            def("get", matrix_get2(&Matrix::get), "docstring").
            def("set", matrix_set1(&Matrix::set), "docstring").
            def("set", matrix_set3(&Matrix::set), "docstring").
            def("set", matrix_set4(&Matrix::set), "docstring").
            def("set", &Matrix::set_by_python_list, "docstring").
            def("project_out", &Matrix::project_out, "docstring").
            def("__getitem__", &Matrix::pyget, "docstring").
            def("__setitem__", &Matrix::pyset, "docstring").
            def("save", matrix_save(&Matrix::save), "docstring").
            def("load", matrix_load(&Matrix::load), "docstring").
            def("load_mpqc", matrix_load(&Matrix::load_mpqc), "docstring").
            def("remove_symmetry", &Matrix::remove_symmetry, "docstring").
            def("symmetrize_gradient", &Matrix::symmetrize_gradient, "docstring").
            def("rotate_columns", &Matrix::rotate_columns, "docstring").
            def("array_interface", [](Matrix &m){
                // Dont ask, hopefully pybind11 will work on this
                py::list ret;
                std::vector<py::buffer_info> buff_vec(m.array_interface());
                std::string typestr = "<";
                {
                   std::stringstream sstr;
                   sstr << (int)sizeof(double);
                   typestr += "f" + sstr.str();
                }

                for (auto const& buff: m.array_interface()){
                    py::dict interface;
                    interface["typestr"] = py::str(typestr);
                    interface["data"] = py::make_tuple((long)buff.ptr, false);
                    py::list shape;
                    for (auto const& s: buff.shape){
                        shape.append(py::int_(s));
                    }
                    interface["shape"] = shape;
                    ret.append(interface);

                }
                return ret;
            });

    py::class_<View>(m, "View").
            def(py::init<SharedMatrix, const Dimension&, const Dimension&>()).
            def(py::init<SharedMatrix, const Dimension&, const Dimension&, const Dimension&, const Dimension&>()).
            def("__call__", &View::operator(), "docstring");

    py::class_<Deriv, std::shared_ptr<Deriv> >(m, "Deriv", "docstring").
            def(py::init<std::shared_ptr<Wavefunction> >()).
            def(py::init<std::shared_ptr<Wavefunction>, char, bool, bool>()).
            def("set_tpdm_presorted", &Deriv::set_tpdm_presorted, "docstring").
            def("set_ignore_reference", &Deriv::set_ignore_reference, "docstring").
            def("set_deriv_density_backtransformed", &Deriv::set_deriv_density_backtransformed, "docstring").
            def("compute", &Deriv::compute, "docstring");

    typedef SharedMatrix (MatrixFactory::*create_shared_matrix)();
    typedef SharedMatrix (MatrixFactory::*create_shared_matrix_name)(const std::string&);

    py::class_<MatrixFactory, std::shared_ptr<MatrixFactory> >(m, "MatrixFactory", "docstring").
            def("create_matrix", create_shared_matrix(&MatrixFactory::create_shared_matrix), "docstring").
            def("create_matrix", create_shared_matrix_name(&MatrixFactory::create_shared_matrix), "docstring");

    py::class_<CdSalcList, std::shared_ptr<CdSalcList>>(m, "CdSalcList", "docstring").
            def("print_out", &CdSalcList::print, "docstring").
            def("matrix", &CdSalcList::matrix, "docstring");

    py::class_<GaussianShell, std::shared_ptr<GaussianShell> >(m, "GaussianShell", "docstring").
            def_property_readonly("nprimitive", py::cpp_function(&GaussianShell::nprimitive), "docstring").
            def_property_readonly("nfunction", py::cpp_function(&GaussianShell::nfunction), "docstring").
            def_property_readonly("ncartesian", py::cpp_function(&GaussianShell::ncartesian), "docstring").
            def_property_readonly("am", py::cpp_function(&GaussianShell::am), "docstring").
            def_property_readonly("amchar", py::cpp_function(&GaussianShell::amchar), "docstring").
            def_property_readonly("AMCHAR", py::cpp_function(&GaussianShell::AMCHAR), "docstring").
            def_property_readonly("ncenter", py::cpp_function(&GaussianShell::ncenter), "docstring").
            def_property("function_index", py::cpp_function(&GaussianShell::function_index), py::cpp_function(&GaussianShell::set_function_index), "Basis function index where this shell starts.").
//            add_property("center", &GaussianShell::center, "A double* representing the center of the GaussianShell.").
//            add_property("exps", &GaussianShell::exps, "The exponents of all the primitives").
//            add_property("coefs", &GaussianShell::coefs, "The coefficients of all the primitives").
            def("is_cartesian", &GaussianShell::is_cartesian, "docstring").
            def("is_pure", &GaussianShell::is_pure, "docstring").
//            def("normalize_shell", &GaussianShell::normalize_shell, "docstring").
            def("exp", &GaussianShell::exp, "Returns the exponent of the given primitive").
            def("original_coef", &GaussianShell::original_coef, "docstring").
            def("erd_coef", &GaussianShell::erd_coef, "docstring").
            def("coef", &GaussianShell::coef, "docstring");

    py::enum_<PrimitiveType>(m, "PrimitiveType","docstring")
       .value("Normalized",Normalized)
       .value("Unnormalized",Unnormalized)
       .export_values()
       ;

    py::enum_ <GaussianType>(m, "GaussianType","docstring")
       .value("Cartesian",Cartesian)
       .value("Pure",Pure)
       .export_values()
       ;

    py::class_<ShellInfo, std::shared_ptr<ShellInfo>>(m, "ShellInfo")
        .def(py::init<int,
                  const std::vector<double>&,
                  const std::vector<double>&,
                  GaussianType,
                  int,
                  const Vector3&,
                  int,
                  PrimitiveType>());

    py::bind_vector<std::vector<ShellInfo>>(m, "BSVec");

    py::class_<OneBodyAOInt, std::shared_ptr<OneBodyAOInt>> pyOneBodyAOInt(m, "OneBodyAOInt", "docstring");
            pyOneBodyAOInt.
            def("compute_shell", &OneBodyAOInt::compute_shell, "docstring").
            def_property("origin", py::cpp_function(&OneBodyAOInt::origin), py::cpp_function(&OneBodyAOInt::set_origin), "The origin about which the one body ints are being computed.").
            def_property_readonly("basis", py::cpp_function(&OneBodyAOInt::basis), "The basis set on center one").
            def_property_readonly("basis1", py::cpp_function(&OneBodyAOInt::basis1), "The basis set on center one").
            def_property_readonly("basis2", py::cpp_function(&OneBodyAOInt::basis2), "The basis set on center two"); // <-- Added semicolon

    //typedef void (OneBodySOInt::*matrix_version)(SharedMatrix) const;
    //typedef void (OneBodySOInt::*vector_version)(std::vector<SharedMatrix>) const;
    //class_<OneBodySOInt, std::shared_ptr<OneBodySOInt>, boost::noncopyable>("OneBodySOInt", "docstring", no_init).
    //        def("compute", matrix_version(&OneBodySOInt::compute_shell), "docstring").
    //        def("compute_list", vector_version(&OneBodySOInt::compute), "docstring").
    //        add_property("basis", &OneBodySOInt::basis, "The basis set on center one").
    //        add_property("basis1", &OneBodySOInt::basis1, "The basis set on center one").
    //        add_property("basis2", &OneBodySOInt::basis2, "The basis set on center two");

    py::class_<OverlapInt, std::shared_ptr<OverlapInt>>(m, "OverlapInt", pyOneBodyAOInt, "docstring");
    py::class_<DipoleInt, std::shared_ptr<DipoleInt>>(m, "DipoleInt", pyOneBodyAOInt, "docstring");
    py::class_<QuadrupoleInt, std::shared_ptr<QuadrupoleInt>>(m, "QuadrupoleInt", pyOneBodyAOInt, "docstring");
    py::class_<MultipoleInt, std::shared_ptr<MultipoleInt>>(m, "MultipoleInt", pyOneBodyAOInt, "docstring");
    py::class_<TracelessQuadrupoleInt, std::shared_ptr<TracelessQuadrupoleInt>>(m, "TracelessQuadrupoleInt", pyOneBodyAOInt, "docstring");
    py::class_<ElectricFieldInt, std::shared_ptr<ElectricFieldInt>>(m, "ElectricFieldInt", pyOneBodyAOInt, "docstring");
    py::class_<KineticInt, std::shared_ptr<KineticInt>>(m, "KineticInt", pyOneBodyAOInt, "docstring");
    py::class_<PotentialInt, std::shared_ptr<PotentialInt>>(m, "PotentialInt", pyOneBodyAOInt, "docstring");
    py::class_<PseudospectralInt, std::shared_ptr<PseudospectralInt>>(m, "PseudospectralInt", pyOneBodyAOInt, "docstring");
    py::class_<ElectrostaticInt, std::shared_ptr<ElectrostaticInt>>(m, "ElectrostaticInt", pyOneBodyAOInt, "docstring");
    py::class_<NablaInt, std::shared_ptr<NablaInt>>(m, "NablaInt", pyOneBodyAOInt, "docstring");
    py::class_<AngularMomentumInt, std::shared_ptr<AngularMomentumInt>>(m, "AngularMomentumInt", pyOneBodyAOInt, "docstring");

    typedef size_t (TwoBodyAOInt::*compute_shell_ints)(int, int, int, int);
    py::class_<TwoBodyAOInt, std::shared_ptr<TwoBodyAOInt>> pyTwoBodyAOInt(m, "TwoBodyAOInt", "docstring");
            pyTwoBodyAOInt.def("compute_shell", compute_shell_ints(&TwoBodyAOInt::compute_shell), "docstring"); // <-- Semicolon

    py::class_<TwoElectronInt, std::shared_ptr<TwoElectronInt>>(m, "TwoElectronInt", pyTwoBodyAOInt, "docstring").
            def("compute_shell", compute_shell_ints(&TwoBodyAOInt::compute_shell), "docstring");

    py::class_<ERI, std::shared_ptr<ERI>>(m, "ERI", pyTwoBodyAOInt, "docstring");
    py::class_<F12, std::shared_ptr<F12>>(m, "F12", pyTwoBodyAOInt, "docstring");
    py::class_<F12G12, std::shared_ptr<F12G12>>(m, "F12G12", pyTwoBodyAOInt, "docstring");
    py::class_<F12Squared, std::shared_ptr<F12Squared>>(m, "F12Squared", pyTwoBodyAOInt, "docstring");
    py::class_<F12DoubleCommutator, std::shared_ptr<F12DoubleCommutator>>(m, "F12DoubleCommutator", pyTwoBodyAOInt, "docstring");
    py::class_<ErfERI, std::shared_ptr<ErfERI>>(m, "ErfERI", pyTwoBodyAOInt, "docstring");
    py::class_<ErfComplementERI, std::shared_ptr<ErfComplementERI>>(m, "ErfComplementERI", pyTwoBodyAOInt, "docstring");

    py::class_<AOShellCombinationsIterator, std::shared_ptr<AOShellCombinationsIterator>>(m, "AOShellCombinationsIterator").
            def_property_readonly("p", py::cpp_function(&AOShellCombinationsIterator::p), "docstring").
            def_property_readonly("q", py::cpp_function(&AOShellCombinationsIterator::q), "docstring").
            def_property_readonly("r", py::cpp_function(&AOShellCombinationsIterator::r), "docstring").
            def_property_readonly("s", py::cpp_function(&AOShellCombinationsIterator::s), "docstring").
            def("first", &AOShellCombinationsIterator::first, "docstring").
            def("next", &AOShellCombinationsIterator::next, "docstring").
            def("is_done", &AOShellCombinationsIterator::is_done, "docstring");

    py::class_<ThreeCenterOverlapInt, std::shared_ptr<ThreeCenterOverlapInt>>(m, "ThreeCenterOverlapInt", "docstring").
            def("compute_shell", &ThreeCenterOverlapInt::compute_shell, "docstring");

    py::class_<IntegralFactory, std::shared_ptr<IntegralFactory>>(m, "IntegralFactory", "docstring").
            def(py::init<std::shared_ptr<BasisSet>, std::shared_ptr<BasisSet>, std::shared_ptr<BasisSet>, std::shared_ptr<BasisSet> >()).
            def(py::init<std::shared_ptr<BasisSet>, std::shared_ptr<BasisSet> >()).
            def(py::init<std::shared_ptr<BasisSet> >()).
            //def("shells_iterator", &IntegralFactory::shells_iterator_ptr, py::return_value_policy<manage_new_object>(), "docstring").
            def("shells_iterator", &IntegralFactory::shells_iterator_ptr, "docstring").
            def("eri", &IntegralFactory::eri, "docstring", py::arg("deriv") = 0, py::arg("use_shell_pairs") = true).
            def("f12", &IntegralFactory::f12, "docstring", py::arg("cf"), py::arg("deriv") = 0, py::arg("use_shell_pairs") = true).
            def("f12g12", &IntegralFactory::f12g12, "docstring", py::arg("cf"), py::arg("deriv") = 0, py::arg("use_shell_pairs") = true).
            def("f12_double_commutator", &IntegralFactory::f12_double_commutator, "docstring", py::arg("cf"), py::arg("deriv") = 0, py::arg("use_shell_pairs") = true).
            def("f12_squared", &IntegralFactory::f12_squared, "docstring", py::arg("cf"), py::arg("deriv") = 0, py::arg("use_shell_pairs") = true).
            def("erf_eri", &IntegralFactory::erf_eri, "docstring", py::arg("omega"), py::arg("deriv") = 0, py::arg("use_shell_pairs") = true).
            def("erf_complement_eri", &IntegralFactory::erf_complement_eri, "docstring", py::arg("omega"), py::arg("deriv") = 0, py::arg("use_shell_pairs") = true).
            def("ao_overlap", &IntegralFactory::ao_overlap, "docstring", py::arg("deriv") = 0).
            def("so_overlap", &IntegralFactory::so_overlap, "docstring", py::arg("deriv") = 0).
            def("ao_dipole", &IntegralFactory::ao_dipole, "docstring", py::arg("deriv") = 0).
            def("so_dipole", &IntegralFactory::so_dipole, "docstring", py::arg("deriv") = 0).
            def("ao_kinetic", &IntegralFactory::ao_kinetic, "docstring", py::arg("deriv") = 0).
            def("so_kinetic", &IntegralFactory::so_kinetic, "docstring", py::arg("deriv") = 0).
            def("ao_potential", &IntegralFactory::ao_potential, "docstring", py::arg("deriv") = 0).
            def("so_potential", &IntegralFactory::so_potential, "docstring", py::arg("deriv") = 0).
            def("ao_pseudospectral", &IntegralFactory::ao_pseudospectral, "docstring", py::arg("deriv") = 0).
            def("so_pseudospectral", &IntegralFactory::so_pseudospectral, "docstring", py::arg("deriv") = 0).
            def("ao_nabla", &IntegralFactory::ao_nabla, "docstring", py::arg("deriv") = 0).
            def("so_nabla", &IntegralFactory::so_nabla, "docstring", py::arg("deriv") = 0).
            def("ao_angular_momentum", &IntegralFactory::ao_angular_momentum, "docstring", py::arg("deriv") = 0).
            def("so_angular_momentum", &IntegralFactory::so_angular_momentum, "docstring", py::arg("deriv") = 0).
            def("ao_quadrupole", &IntegralFactory::ao_quadrupole, "docstring").
            def("so_quadrupole", &IntegralFactory::so_quadrupole, "docstring").
            def("ao_multipoles", &IntegralFactory::ao_multipoles, "docstring", py::arg("order")).
            def("so_multipoles", &IntegralFactory::so_multipoles, "docstring", py::arg("order")).
            def("ao_traceless_quadrupole", &IntegralFactory::ao_traceless_quadrupole, "docstring").
            def("so_traceless_quadrupole", &IntegralFactory::so_traceless_quadrupole, "docstring").
            def("electric_field", &IntegralFactory::electric_field, "docstring").
            def("electrostatic", &IntegralFactory::electrostatic, "docstring").
            def("overlap_3c", &IntegralFactory::overlap_3c, "docstring");

    typedef std::shared_ptr<PetiteList> (MintsHelper::*petite_list_0)() const;
    typedef std::shared_ptr<PetiteList> (MintsHelper::*petite_list_1)(bool) const;

    typedef SharedMatrix (MintsHelper::*erf)(double, SharedMatrix, SharedMatrix, SharedMatrix, SharedMatrix);
    typedef SharedMatrix (MintsHelper::*eri)(SharedMatrix, SharedMatrix, SharedMatrix, SharedMatrix);
    typedef SharedMatrix (MintsHelper::*normal_eri)();
    typedef SharedMatrix (MintsHelper::*normal_eri2)(std::shared_ptr<BasisSet>,std::shared_ptr<BasisSet>,std::shared_ptr<BasisSet>,std::shared_ptr<BasisSet>);
    typedef SharedMatrix (MintsHelper::*normal_3c)(std::shared_ptr<BasisSet>,std::shared_ptr<BasisSet>,std::shared_ptr<BasisSet>);

    typedef SharedMatrix (MintsHelper::*normal_f12)(std::shared_ptr<CorrelationFactor>);
    typedef SharedMatrix (MintsHelper::*normal_f122)(std::shared_ptr<CorrelationFactor>, std::shared_ptr<BasisSet>,std::shared_ptr<BasisSet>,std::shared_ptr<BasisSet>,std::shared_ptr<BasisSet>);

    typedef SharedMatrix (MintsHelper::*oneelectron)();
    typedef SharedMatrix (MintsHelper::*oneelectron_mixed_basis)(std::shared_ptr<BasisSet>, std::shared_ptr<BasisSet>);

    py::class_<MintsHelper, std::shared_ptr<MintsHelper> >(m, "MintsHelper", "docstring").
            def(py::init<std::shared_ptr<BasisSet> >()).

            // Options and attributes
            def("nbf", &MintsHelper::nbf, "docstring").
            def("set_print", &MintsHelper::set_print, "docstring").
            def("basisset", &MintsHelper::basisset, "docstring").
            def("sobasisset", &MintsHelper::sobasisset, "docstring").
            def("factory", &MintsHelper::factory, "docstring").
            def("cdsalcs", &MintsHelper::cdsalcs, "docstring").
            def("petite_list", petite_list_0(&MintsHelper::petite_list), "docstring").
            def("petite_list1", petite_list_1(&MintsHelper::petite_list), "docstring").

            // Integral builders
            def("integral", &MintsHelper::integral, "docstring").
            def("integrals", &MintsHelper::integrals, "docstring").
            def("integrals_erf", &MintsHelper::integrals_erf, "docstring").
            def("integrals_erfc", &MintsHelper::integrals_erfc, "docstring").
            def("one_electron_integrals", &MintsHelper::one_electron_integrals, "docstring").
            def("one_electron_integrals", &MintsHelper::one_electron_integrals, "docstring").

            // One-electron
            def("ao_overlap", oneelectron(&MintsHelper::ao_overlap), "docstring").
            def("ao_overlap", oneelectron_mixed_basis(&MintsHelper::ao_overlap), "docstring").
            def("so_overlap", &MintsHelper::so_overlap, "docstring").
            def("ao_kinetic", oneelectron(&MintsHelper::ao_kinetic), "docstring").
            def("ao_kinetic", oneelectron_mixed_basis(&MintsHelper::ao_kinetic), "docstring").
            def("so_kinetic", &MintsHelper::so_kinetic, "docstring").
            def("ao_potential", oneelectron(&MintsHelper::ao_potential), "docstring").
            def("ao_potential", oneelectron_mixed_basis(&MintsHelper::ao_potential), "docstring").
            def("so_potential", &MintsHelper::so_potential, "docstring").

            // One-electron properties and
            def("ao_pvp", &MintsHelper::ao_pvp, "docstring").
            def("ao_dkh", &MintsHelper::ao_dkh, "docstring").
            def("so_dkh", &MintsHelper::so_dkh, "docstring").
            def("ao_dipole", &MintsHelper::ao_dipole, "docstring").
            def("so_dipole", &MintsHelper::so_dipole, "docstring").
            def("so_quadrupole", &MintsHelper::so_quadrupole, "docstring").
            def("so_traceless_quadrupole", &MintsHelper::so_traceless_quadrupole, "docstring").
            def("ao_nabla", &MintsHelper::ao_nabla, "docstring").
            def("so_nabla", &MintsHelper::so_nabla, "docstring").
            def("ao_angular_momentum", &MintsHelper::ao_angular_momentum, "docstring").
            def("so_angular_momentum", &MintsHelper::so_angular_momentum, "docstring").

            // Two-electron AO
            def("ao_eri", normal_eri(&MintsHelper::ao_eri), "docstring").
            def("ao_eri", normal_eri2(&MintsHelper::ao_eri), "docstring").
            def("ao_eri_shell", &MintsHelper::ao_eri_shell, "docstring").
            def("ao_erf_eri", &MintsHelper::ao_erf_eri, "docstring").
            def("ao_f12", normal_f12(&MintsHelper::ao_f12), "docstring").
            def("ao_f12", normal_f122(&MintsHelper::ao_f12), "docstring").
            def("ao_f12_scaled", normal_f12(&MintsHelper::ao_f12_scaled), "docstring").
            def("ao_f12_scaled", normal_f122(&MintsHelper::ao_f12_scaled), "docstring").
            def("ao_f12_squared", normal_f12(&MintsHelper::ao_f12_squared), "docstring").
            def("ao_f12_squared", normal_f122(&MintsHelper::ao_f12_squared), "docstring").
            def("ao_f12g12", &MintsHelper::ao_f12g12, "docstring").
            def("ao_f12_double_commutator", &MintsHelper::ao_f12_double_commutator, "docstring").
            def("ao_3coverlap", normal_eri(&MintsHelper::ao_3coverlap), "docstring").
            def("ao_3coverlap", normal_3c(&MintsHelper::ao_3coverlap), "docstring").

            // Two-electron MO and transformers
            def("mo_eri", eri(&MintsHelper::mo_eri), "docstring").
            def("mo_erf_eri", erf(&MintsHelper::mo_erf_eri), "docstring").
            def("mo_f12", &MintsHelper::mo_f12, "docstring").
            def("mo_f12_squared", &MintsHelper::mo_f12_squared, "docstring").
            def("mo_f12g12", &MintsHelper::mo_f12g12, "docstring").
            def("mo_f12_double_commutator", &MintsHelper::mo_f12_double_commutator, "docstring").
            def("mo_spin_eri", &MintsHelper::mo_spin_eri, "docstring").
            def("mo_transform", &MintsHelper::mo_transform, "docstring").
            def("set_rel_basisset", &MintsHelper::set_rel_basisset, "docstring").

            def("play", &MintsHelper::play, "docstring");

    py::class_<Vector3>(m, "Vector3", "Class for vectors of length three, often Cartesian coordinate vectors, and their common operations").
            def(py::init<double>()).
            def(py::init<double, double, double>()).
            def(py::init<const Vector3&>()).
            //      def(self = other<double>()).
            def(py::self += py::self).
            def(py::self -= py::self).
            def(py::self *= double()).
            //def(py::self + py::self).
            //def(py::self - py::self).
            //def(-py::self).
            def("dot", &Vector3::dot, "Returns dot product of arg1 and arg2").
            def("distance", &Vector3::distance, "Returns distance between two points represented by arg1 and arg2").
            def("normalize", &Vector3::normalize, "Returns vector of unit length and arg1 direction").
            def("norm", &Vector3::norm, "Returns Euclidean norm of arg1").
            def("cross", &Vector3::cross, "Returns cross product of arg1 and arg2").
            def("__str__", &Vector3::to_string, "Returns a string representation of arg1, suitable for printing.").
            def("__getitem__", &Vector3::get, "Returns the arg2-th element of arg1.");

    typedef void (SymmetryOperation::*intFunction)(int);
    typedef void (SymmetryOperation::*doubleFunction)(double);

    py::class_<SymmetryOperation>(m, "SymmetryOperation", "Class to provide a 3 by 3 matrix representation of a symmetry operation, such as a rotation or reflection.").
            def(py::init<const SymmetryOperation& >()).
            def("trace", &SymmetryOperation::trace, "Returns trace of transformation matrix").
            def("zero", &SymmetryOperation::zero, "Zero out the symmetry operation").
            def("operate", &SymmetryOperation::operate, "Performs the operation arg2 * arg1").
            def("transform", &SymmetryOperation::transform, "Performs the transform arg2 * arg1 * arg2~").
            def("unit", &SymmetryOperation::unit, "Set equal to a unit matrix").
            def("E", &SymmetryOperation::E, "Set equal to E").
            def("i", &SymmetryOperation::i, "Set equal to an inversion").
            def("sigma_xy", &SymmetryOperation::sigma_xy, "Set equal to reflection in xy plane").
            def("sigma_yz", &SymmetryOperation::sigma_yz, "Set equal to reflection in yz plane").
            def("sigma_xz", &SymmetryOperation::sigma_xz, "Set equal to reflection in xz plane").
            //        def("sigma_yz", &SymmetryOperation::sigma_yz).
            def("rotate_n", intFunction(&SymmetryOperation::rotation), "Set equal to a clockwise rotation by 2pi/n").
            def("rotate_theta", doubleFunction(&SymmetryOperation::rotation), "Set equal to a clockwise rotation by theta").
            def("c2_x", &SymmetryOperation::c2_x, "Set equal to C2 about the x axis").
            def("c2_y", &SymmetryOperation::c2_y, "Set equal to C2 about the y axis").
            def("c2_z", &SymmetryOperation::c2_z, "Set equal to C2 about the z axis").
            def("transpose", &SymmetryOperation::transpose, "Performs transposition of matrix operation");

    py::class_<OrbitalSpace>(m, "OrbitalSpace", "docstring").
            def(py::init<const std::string&, const std::string&, const SharedMatrix&, const SharedVector&, const std::shared_ptr<BasisSet>&, const std::shared_ptr<IntegralFactory>& >()).
            def(py::init<const std::string&, const std::string&, const SharedMatrix&, const std::shared_ptr<BasisSet>&, const std::shared_ptr<IntegralFactory>& >()).
            def(py::init<const std::string&, const std::string&, const std::shared_ptr<Wavefunction>& >()).
            def("nirrep", &OrbitalSpace::nirrep, "docstring").
            def("id", &OrbitalSpace::id, "docstring").
            def("name", &OrbitalSpace::name, "docstring").
            def("C", &OrbitalSpace::C, "docstring").
            def("evals", &OrbitalSpace::evals, "docstring").
            def("basisset", &OrbitalSpace::basisset, "docstring").
            def("integral", &OrbitalSpace::integral, "docstring").
            def("dim", &OrbitalSpace::dim, "docstring").
            def("print_out", &OrbitalSpace::print, "docstring").
            def_static("build_cabs_space", &OrbitalSpace::build_cabs_space, "docstring").
            def_static("build_ri_space", &OrbitalSpace::build_ri_space, "docstring");

    py::class_<PointGroup, std::shared_ptr<PointGroup> >(m, "PointGroup", "docstring").
            def(py::init<const std::string&>()).
            def("symbol", &PointGroup::symbol, "Returns Schoenflies symbol for point group");
            //def("origin", &PointGroup::origin).
//            def("set_symbol", &PointGroup::set_symbol);

    typedef void (Molecule::*matrix_set_geometry)(const Matrix &);
    typedef Vector3 (Molecule::*nuclear_dipole1)(const Vector3&) const;
    typedef Vector3 (Molecule::*nuclear_dipole2)() const;

    py::class_<Molecule, std::shared_ptr<Molecule> >(m, "Molecule", "Class to store the elements, coordinates, fragmentation pattern, basis sets, charge, multiplicity, etc. of a molecule.").
            def("set_geometry", matrix_set_geometry(&Molecule::set_geometry), "Sets the geometry, given a (Natom X 3) matrix arg2 of coordinates (in Bohr)").
            def("set_name", &Molecule::set_name, "Sets molecule name").
            def("nuclear_dipole", nuclear_dipole1(&Molecule::nuclear_dipole), "Gets the nuclear contribution to the dipole, withe respect to a specified origin").
            def("nuclear_dipole", nuclear_dipole2(&Molecule::nuclear_dipole), "Gets the nuclear contribution to the dipole, withe respect to the origin").
            def("name", &Molecule::name, "Gets molecule name").
            def("reinterpret_coordentry", &Molecule::set_reinterpret_coordentry, "Do reinterpret coordinate entries during update_geometry().").
            def("fix_orientation", &Molecule::set_orientation_fixed, "Fix the orientation at its current frame").
            def("fix_com", &Molecule::set_com_fixed, "Whether to fix the Cartesian position, or to translate to the C.O.M.").
            def("add_atom", &Molecule::add_atom, "Adds to Molecule arg1 an atom with atomic number arg2, Cartesian coordinates in Bohr (arg3, arg4, arg5), atomic symbol arg6, mass arg7, charge arg8 (optional), and lineno arg9 (optional)").
            def("natom", &Molecule::natom, "Number of real atoms").
            def("multiplicity", &Molecule::multiplicity, "Gets the multiplicity (defined as 2Ms + 1)").
            def("nfragments", &Molecule::nfragments, "Gets the number of fragments in the molecule").
            def("set_nuclear_charge", &Molecule::set_nuclear_charge, "Set the nuclear charge of the given atom to the value provided.").
            def("basis_on_atom", &Molecule::basis_on_atom, py::return_value_policy::copy, "Gets the label of the orbital basis set on a given atom.").
            def("print_in_input_format", &Molecule::print_in_input_format, "Prints the molecule as Cartesian or ZMatrix entries, just as inputted.").
            def("create_psi4_string_from_molecule", &Molecule::create_psi4_string_from_molecule, "Gets a string reexpressing in input format the current states of the molecule").
            def("save_xyz_file", &Molecule::save_xyz_file, "Saves an XYZ file to arg2").
            def("save_string_xyz_file", &Molecule::save_string_xyz_file, "Saves an XYZ file to arg2").
            def("save_string_xyz", &Molecule::save_string_xyz, "Saves the string of an XYZ file to arg2").
            def("Z", &Molecule::Z, py::return_value_policy::copy, "Nuclear charge of atom").
            def("x", &Molecule::x, "x position of atom").
            def("y", &Molecule::y, "y position of atom").
            def("z", &Molecule::z, "z position of atom").
            //def("xyz", &Molecule::xyz).
            def("center_of_mass", &Molecule::center_of_mass, "Computes center of mass of molecule (does not translate molecule)").
            def("translate", &Molecule::translate, "Translates molecule by arg2").
            def("move_to_com", &Molecule::move_to_com, "Moves molecule to center of mass").
            def("mass", &Molecule::mass, "Gets mass of atom arg2").
            def("set_mass", &Molecule::set_mass, "Gets mass of atom arg2").
            def("symbol", &Molecule::symbol, "Gets the cleaned up label of atom arg2 (C2 => C, H4 = H)").
            def("label", &Molecule::label, "Gets the original label of the atom as given in the input file (C2, H4)").
            def("charge", &Molecule::charge, "Gets charge of atom").
            def("molecular_charge", &Molecule::molecular_charge, "Gets the molecular charge").
            def("extract_subsets", &Molecule::py_extract_subsets_1, "Returns copy of arg1 with arg2 fragments Real and arg3 fragments Ghost").
            def("extract_subsets", &Molecule::py_extract_subsets_2, "Returns copy of arg1 with arg2 fragments Real and arg3 fragment Ghost").
            def("extract_subsets", &Molecule::py_extract_subsets_3, "Returns copy of arg1 with arg2 fragment Real and arg3 fragments Ghost").
            def("extract_subsets", &Molecule::py_extract_subsets_4, "Returns copy of arg1 with arg2 fragment Real and arg3 fragment Ghost").
            def("extract_subsets", &Molecule::py_extract_subsets_5, "Returns copy of arg1 with arg2 fragments Real").
            def("extract_subsets", &Molecule::py_extract_subsets_6, "Returns copy of arg1 with arg2 fragment Real").
            def("activate_all_fragments", &Molecule::activate_all_fragments, "Sets all fragments in the molecule to be active").
            def("deactivate_all_fragments", &Molecule::deactivate_all_fragments, "Sets all fragments in the molecule to be inactive").
            def("set_active_fragments", &Molecule::set_active_fragments, "Sets the specified list arg2 of fragments to be Real").
            def("set_active_fragment", &Molecule::set_active_fragment, "Sets the specified fragment arg2 to be Real").
            def("set_ghost_fragments", &Molecule::set_ghost_fragments, "Sets the specified list arg2 of fragments to be Ghost").
            def("set_ghost_fragment", &Molecule::set_ghost_fragment, "Sets the specified fragment arg2 to be Ghost").
            def("atom_at_position", &Molecule::atom_at_position1, "Tests to see if an atom is at the position arg2 with a given tolerance arg3").
            def("print_out", &Molecule::print, "Prints the molecule in Cartesians in input units").
            def("print_out_in_bohr", &Molecule::print_in_bohr, "Prints the molecule in Cartesians in Bohr").
            def("print_out_in_angstrom", &Molecule::print_in_angstrom, "Prints the molecule in Cartesians in Angstroms").
            def("print_cluster", &Molecule::print_cluster, "Prints the molecule in Cartesians in input units adding fragment separators").
            def("rotational_constants", &Molecule::rotational_constants, "Prints the rotational constants of the molecule").
            def("nuclear_repulsion_energy", &Molecule::nuclear_repulsion_energy, "Computes nuclear repulsion energy").
            def("find_point_group", &Molecule::find_point_group, "Finds computational molecular point group, user can override this with the symmetry keyword").
            def("reset_point_group", &Molecule::reset_point_group, "Overrides symmetry from outside the molecule string").
            def("set_point_group", &Molecule::set_point_group, "Sets the molecular point group to the point group object arg2").
            def("get_full_point_group", &Molecule::full_point_group, "Gets point group name such as C3v or S8").
            def("point_group", &Molecule::point_group, "Returns the current point group object").
            def("schoenflies_symbol", &Molecule::schoenflies_symbol, "Returns the Schoenflies symbol").
            def("form_symmetry_information", &Molecule::form_symmetry_information, "Uses the point group object obtain by calling point_group()").
            def("symmetrize", &Molecule::symmetrize_to_abelian_group, "Finds the highest point Abelian point group within the specified tolerance, and forces the geometry to have that symmetry.").
            def_static("create_molecule_from_string", &Molecule::create_molecule_from_string, "Returns a new Molecule with member data from the geometry string arg1 in psi4 format").
            def("is_variable", &Molecule::is_variable, "Checks if variable arg2 is in the list, returns true if it is, and returns false if not").
            def("set_variable", &Molecule::set_variable, "Assigns the value arg3 to the variable arg2 in the list of geometry variables, then calls update_geometry()").
            def("get_variable", &Molecule::get_variable, "Checks if variable arg2 is in the list, sets it to val and returns true if it is, and returns false if not").
            def("update_geometry", &Molecule::update_geometry, "Reevaluates the geometry with current variable values, orientation directives, etc. Must be called after initial Molecule definition by string.").
            def("set_molecular_charge", &Molecule::set_molecular_charge, "Sets the molecular charge").
            def("set_multiplicity", &Molecule::set_multiplicity, "Sets the multiplicity (defined as 2Ms + 1)").
            def("set_basis_all_atoms", &Molecule::set_basis_all_atoms, "Sets basis set arg2 to all atoms").
            def("set_basis_by_symbol", &Molecule::set_basis_by_symbol, "Sets basis set arg3 to all atoms with symbol (e.g., H) arg2").
            def("set_basis_by_label", &Molecule::set_basis_by_label, "Sets basis set arg3 to all atoms with label (e.g., H4) arg2").
            //def("set_basis_by_number", &Molecule::set_basis_by_number, "Sets basis set arg3 to atom number (1-indexed, incl. dummies) arg2").  // dangerous for user use
            def("irrep_labels", [](Molecule &mol){
                std::vector<std::string> ret;
                char **labels = mol.irrep_labels();
                int nirrep = mol.point_group()->char_table().nirrep();
                for (size_t h = 0; h<nirrep; h++){
                    std::string lh(labels[h]);
                    ret.push_back(lh);
                }
                return ret;
            }).
            def_property("units", py::cpp_function(&Molecule::units), py::cpp_function(&Molecule::set_units), "Units (Angstrom or Bohr) used to define the geometry").
            def("clone", &Molecule::clone, "Returns a new Molecule identical to arg1").
            def("geometry", &Molecule::geometry, "Gets the geometry as a (Natom X 3) matrix of coordinates (in Bohr)");

    py::class_<PetiteList, std::shared_ptr<PetiteList>>(m, "PetiteList", "docstring").
            def("aotoso", &PetiteList::aotoso, "docstring").
            def("sotoao", &PetiteList::sotoao, "docstring").
            def("print", &PetiteList::print, "docstring");

    py::class_<BasisSetParser, std::shared_ptr<BasisSetParser>>(m, "BasisSetParser", "docstring");
    py::class_<Gaussian94BasisSetParser, std::shared_ptr<Gaussian94BasisSetParser>, BasisSetParser>(m, "Gaussian94BasisSetParser", "docstring");

    //py::bind_map<std::string, std::vector>(m, "ShellInfo");

    //class_<ShellInfoMap>("ShellInfoMap")
    //   .def(map_indexing_suite<ShellInfoMap>());

    //using ShellInfoMapMap = std::map<std::string,ShellInfoMap>;
    //class_<ShellInfoMapMap>("ShellInfoMapMap")
    //   .def(map_indexing_suite<ShellInfoMapMap>());

    typedef void (BasisSet::*basis_print_out)() const;
    typedef const GaussianShell& (BasisSet::*no_center_version)(int) const;
    typedef const GaussianShell& (BasisSet::*center_version)(int, int) const;
    typedef std::shared_ptr<BasisSet> (BasisSet::*ptrversion)(const std::shared_ptr<BasisSet>&) const;
    py::class_<BasisSet, std::shared_ptr<BasisSet>>(m, "BasisSet", "docstring").
            def(py::init<const std::string&, std::shared_ptr<Molecule>, std::map<std::string, std::map<std::string, std::vector<ShellInfo>>>&>()).
            def("name", &BasisSet::name, "Callback handle, may represent string or function").
            def("blend", &BasisSet::target, "Plus-separated string of [basisname] values").
            def("molecule", &BasisSet::molecule, "docstring").
            def("print_out", basis_print_out(&BasisSet::print), "docstring").
            def("print_detail_out", basis_print_out(&BasisSet::print_detail), "docstring").
            def("genbas", &BasisSet::print_detail_cfour, "Returns basis set per atom in CFOUR format").
            def_static("make_filename", &BasisSet::make_filename, "Returns filename for basis name: pluses, stars, parentheses replaced and gbs extension added").
//            def_static("construct", &BasisSet::construct, "docstring").
            def_static("zero_ao_basis_set", &BasisSet::zero_ao_basis_set, "Returns a BasisSet object that actually has a single s-function at the origin with an exponent of 0.0 and contraction of 1.0.").
            def("nbf", &BasisSet::nbf, "Returns number of basis functions (Cartesian or spherical depending on has_puream)").
            def("nao", &BasisSet::nao, "Returns number of atomic orbitals (Cartesian)").
            def("nprimitive", &BasisSet::nprimitive, "Returns total number of primitives in all contractions").
            def("nshell", &BasisSet::nshell, "Returns number of shells").
            def("shell", no_center_version(&BasisSet::shell), py::return_value_policy::copy, "docstring").
            def("shell", center_version(&BasisSet::shell), py::return_value_policy::copy, "docstring").
            def("max_am", &BasisSet::max_am, "Returns maximum angular momentum used").
            def("has_puream", &BasisSet::has_puream, "Spherical harmonics?").
            def("shell_to_basis_function", &BasisSet::shell_to_basis_function, "docstring").
            def("shell_to_center", &BasisSet::shell_to_center, "docstring").
            def("shell_to_ao_function", &BasisSet::shell_to_ao_function, "docstring").
            def("function_to_shell", &BasisSet::function_to_shell, "docstring").
            def("function_to_center", &BasisSet::function_to_center, "Given a function number, return the number of the center it is on.").
            def("nshell_on_center", &BasisSet::nshell_on_center, "docstring").
//            def("decontract", &BasisSet::decontract, "docstring").
            def("ao_to_shell", &BasisSet::ao_to_shell, "docstring").
            def("max_function_per_shell", &BasisSet::max_function_per_shell, "docstring").
            def("max_nprimitive", &BasisSet::max_nprimitive, "docstring").
            def_static("construct_from_pydict", &BasisSet::construct_from_pydict, "docstring").
            def_static("construct_ecp_from_pydict", &BasisSet::construct_ecp_from_pydict, "docstring");


    py::class_<SOBasisSet, std::shared_ptr<SOBasisSet>>(m, "SOBasisSet", "docstring").
            def("petite_list", &SOBasisSet::petite_list, "docstring");

    py::class_<ExternalPotential, std::shared_ptr<ExternalPotential>>(m, "ExternalPotential", "docstring").
            def(py::init<>()).
            def("setName", &ExternalPotential::setName, "docstring").
            def("addCharge", &ExternalPotential::addCharge, "docstring").
            def("addBasis", &ExternalPotential::addBasis, "docstring").
            def("clear", &ExternalPotential::clear, "docstring").
            def("computePotentialMatrix", &ExternalPotential::computePotentialMatrix, "docstring").
            def("print_out", &ExternalPotential::py_print, "docstring");


    typedef std::shared_ptr<Localizer> (*localizer_with_type)(const std::string&, std::shared_ptr<BasisSet>, std::shared_ptr<Matrix>);

    py::class_<Localizer, std::shared_ptr<Localizer>>(m, "Localizer", "docstring").
            def_static("build", localizer_with_type(&Localizer::build), "docstring").
            def("localize", &Localizer::localize, "Perform the localization procedure").
            def_property_readonly("L", py::cpp_function(&Localizer::L), "Localized orbital coefficients").
            def_property_readonly("U", py::cpp_function(&Localizer::U), "Orbital rotation matrix").
            def_property_readonly("converged", py::cpp_function(&Localizer::converged), "Did the localization procedure converge?");

    py::class_<BoysLocalizer, std::shared_ptr<BoysLocalizer>, Localizer>(m, "BoysLocalizer", "docstring");
    py::class_<PMLocalizer, std::shared_ptr<PMLocalizer>, Localizer>(m, "PMLocalizer", "docstring");

    py::class_<FCHKWriter, std::shared_ptr<FCHKWriter> >(m, "FCHKWriter", "docstring").
            def(py::init<std::shared_ptr<Wavefunction> >()).
            def("write", &FCHKWriter::write, "docstring");

    py::class_<MoldenWriter, std::shared_ptr<MoldenWriter> >(m, "MoldenWriter", "docstring").
            def(py::init<std::shared_ptr<Wavefunction> >()).
            def("write", &MoldenWriter::write, "docstring");

    py::class_<NBOWriter, std::shared_ptr<NBOWriter> >(m, "NBOWriter", "docstring").
            def(py::init<std::shared_ptr<Wavefunction> >()).
            def("write", &NBOWriter::write, "docstring");

    py::class_<MOWriter, std::shared_ptr<MOWriter> >(m, "MOWriter", "docstring").
            def(py::init<std::shared_ptr<Wavefunction> >()).
            def("write", &MOWriter::write, "docstring");

    py::class_<OperatorSymmetry, std::shared_ptr<OperatorSymmetry> >(m, "MultipoleSymmetry", "docstring").
            def(py::init<int, const std::shared_ptr<Molecule>&,
                const std::shared_ptr<IntegralFactory>&,
                const std::shared_ptr<MatrixFactory>&>()).
            def("create_matrices", &OperatorSymmetry::create_matrices, "docstring");

    py::class_<CorrelationFactor, std::shared_ptr<CorrelationFactor>>(m, "CorrelationFactor", "docstring").
            def(py::init<unsigned int>()).
            def(py::init<std::shared_ptr<Vector>, std::shared_ptr<Vector> >()).
            def("set_params", &CorrelationFactor::set_params, "docstring");

    py::class_<FittedSlaterCorrelationFactor, CorrelationFactor>(m, "FittedSlaterCorrelationFactor", "docstring").
            def(py::init<double>()).
            def("exponent", &FittedSlaterCorrelationFactor::exponent);

    py::class_<CorrelationTable, std::shared_ptr<CorrelationTable>>(m, "CorrelationTable", "docstring").
            def(py::init<std::shared_ptr<PointGroup>, std::shared_ptr<PointGroup>>()).
            def("group", &CorrelationTable::group, "docstring").
            def("subgroup", &CorrelationTable::subgroup, "docstring").
            def("n", &CorrelationTable::n, "docstring").
            def("subn", &CorrelationTable::subn, "docstring").
            def("degen", &CorrelationTable::degen, "docstring").
            def("subdegen", &CorrelationTable::subdegen, "docstring").
            def("ngamma", &CorrelationTable::ngamma, "docstring").
            def("group", &CorrelationTable::gamma, "docstring");
}
