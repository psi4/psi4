/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2019 The Psi4 Developers.
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

#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>

#include <complex>
#include <string>

#define FORCE_IMPORT_ARRAY              // numpy C api loading
#include "xtensor-python/pytensor.hpp"  // Numpy bindings

#include "psi4/libmints/dimension.h"
#include "psi4/libmints/linalg.h"
#include "psi4/libmints/tensor.h"
#include "psi4/libpsi4util/exception.h"

using namespace psi;
namespace py = pybind11;
using namespace pybind11::literals;

// NOTE py::overload_cast is _broken_ on Intel < 19. However, using
// static_cast here would result in overly verbose binding declarations.

// NOTE Using complex, single-precision data will lead to errors as the bindings
// for std::complex<float> are not in place. Indeed, the C++ standard does not
// allow mixing float and std::complex<double> so I don't know if this can even
// be "programmed around".

namespace {
template <typename T, size_t Rank>
struct Decorator final {
    using Class = Tensor<T, Rank>;
    using PyClass = py::class_<Class, std::shared_ptr<Class>>;

    static void decorate(py::module&, PyClass&) {}
};

template <typename T>
struct Decorator<T, 1> final {
    using Class = Tensor<T, 1>;
    using PyClass = py::class_<Class, std::shared_ptr<Class>>;

    static void decorate(py::module& /* mod */, PyClass& cls) {
        cls.def(py::init<const std::string&, const Dimension&, T>(), "Labeled, blocked vector", "label"_a, "dimpi"_a,
                "fill_value"_a = static_cast<T>(0));
        cls.def(py::init<const std::string&, int, T>(), "Labeled, 1-irrep vector", "label"_a, "dim"_a,
                "fill_value"_a = static_cast<T>(0));
        cls.def(py::init<const Dimension&, T>(), "Unlabeled, blocked vector", "dimpi"_a,
                "fill_value"_a = static_cast<T>(0));
        cls.def(py::init<int, T>(), "Unlabeled, 1-irrep vector", "dim"_a, "fill_value"_a = static_cast<T>(0));
        cls.def_property_readonly("dimpi", [](const Class& obj) { return obj.dimpi(); }, py::return_value_policy::copy,
                                  "Return the Dimension object");
    }
};

template <typename T>
struct Decorator<T, 2> final {
    using Class = Tensor<T, 2>;
    using PyClass = py::class_<Class, std::shared_ptr<Class>>;

    static void decorate(py::module& mod, PyClass& cls) {
        cls.def(py::init<const std::string&, const Dimension&, const Dimension&, unsigned int, T>(),
                "Labeled, blocked, symmetry-assigned matrix", "label"_a, "rowspi"_a, "colspi"_a, "symmetry"_a,
                "fill_value"_a = static_cast<T>(0));
        cls.def(py::init<const std::string&, const Dimension&, const Dimension&, T>(), "Labeled, blocked matrix",
                "label"_a, "rowspi"_a, "colspi"_a, "fill_value"_a = static_cast<T>(0));
        cls.def(py::init<const std::string&, int, int, T>(), "Labeled, 1-irrep matrix", "label"_a, "rows"_a, "cols"_a,
                "fill_value"_a = static_cast<T>(0));
        cls.def(py::init<const Dimension&, const Dimension&, unsigned int, T>(),
                "Unlabeled, blocked, symmetry-assigned matrix", "rowspi"_a, "colspi"_a, "symmetry"_a,
                "fill_value"_a = static_cast<T>(0));
        cls.def(py::init<const Dimension&, const Dimension&, T>(), "Unlabeled, blocked matrix", "rowspi"_a, "colspi"_a,
                "fill_value"_a = static_cast<T>(0));
        cls.def(py::init<int, int, T>(), "Unlabeled, 1-irrep matrix", "rows"_a, "cols"_a,
                "fill_value"_a = static_cast<T>(0));
        cls.def_property_readonly("rowspi", [](const Class& obj) { return obj.rowspi(); },
                                  py::return_value_policy::copy, "Returns the rows per irrep array");
        cls.def("rows", [](const Class& obj, size_t h) { return obj.rows(h); },
                "Returns the number of rows in given irrep", "h"_a = 0);
        cls.def_property_readonly("colspi", [](const Class& obj) { return obj.colspi(); },
                                  py::return_value_policy::copy, "Returns the columns per irrep array");
        cls.def("cols", [](const Class& obj, size_t h) { return obj.cols(h); },
                "Returns the number of columns in given irrep", "h"_a = 0);

        // Bind free functions to module
        declareRank2FreeFunctions(mod);
    }

    static void declareRank2FreeFunctions(py::module& mod) {
        // Type-homogeneous GEMM-s
        mod.def("gemm",
                [](const SharedTensor<T, 2>& A, const SharedTensor<T, 2>& B, Operation opA, Operation opB, double alpha,
                   double beta) { return gemm(opA, opB, alpha, A, B, beta); },
                "A"_a, "B"_a, "opA"_a = Operation::None, "opB"_a = Operation::None, "alpha"_a = 1.0, "beta"_a = 0.0);
        // Type-inhomogeneous GEMM-s
        // T x double
        mod.def("gemm",
                [](const SharedTensor<T, 2>& A, const SharedTensor<double, 2>& B, Operation opA, Operation opB,
                   double alpha, double beta) { return gemm(opA, opB, alpha, A, B, beta); },
                "A"_a, "B"_a, "opA"_a = Operation::None, "opB"_a = Operation::None, "alpha"_a = 1.0, "beta"_a = 0.0);
        // double x T
        mod.def("gemm",
                [](const SharedTensor<double, 2>& A, const SharedTensor<T, 2>& B, Operation opA, Operation opB,
                   double alpha, double beta) { return gemm(opA, opB, alpha, A, B, beta); },
                "A"_a, "B"_a, "opA"_a = Operation::None, "opB"_a = Operation::None, "alpha"_a = 1.0, "beta"_a = 0.0);
        // Type-homogeneous doublet-s
        mod.def(
            "doublet",
            py::overload_cast<const SharedTensor<T, 2>&, const SharedTensor<T, 2>&, Operation, Operation>(&doublet<T>),
            "Returns the multiplication of two matrices A and B, with options to transpose/transpose-conjugate "
            "each beforehand",
            "A"_a, "B"_a, "opA"_a = Operation::None, "opB"_a = Operation::None,
            py::return_value_policy::reference_internal);
        mod.def("doublet",
                py::overload_cast<const SharedTensor<T, 2>&, const SharedTensor<T, 2>&, bool, bool>(&doublet<T>),
                "Returns the multiplication of two matrices A and B, with options to transpose each beforehand", "A"_a,
                "B"_a, "transA"_a = false, "transB"_a = false, py::return_value_policy::reference_internal);
        // Type-inhomogeneous doublet-s
        // T x double
        mod.def("doublet",
                py::overload_cast<const SharedTensor<T, 2>&, const SharedTensor<double, 2>&, Operation, Operation>(
                    &doublet<T, double>),
                "Returns the multiplication of two matrices A and B, with options to transpose/transpose-conjugate "
                "each beforehand",
                "A"_a, "B"_a, "opA"_a = Operation::None, "opB"_a = Operation::None,
                py::return_value_policy::reference_internal);
        mod.def("doublet",
                py::overload_cast<const SharedTensor<T, 2>&, const SharedTensor<double, 2>&, bool, bool>(
                    &doublet<T, double>),
                "Returns the multiplication of two matrices A and B, with options to transpose each beforehand", "A"_a,
                "B"_a, "transA"_a = false, "transB"_a = false, py::return_value_policy::reference_internal);
        // double x T
        mod.def("doublet",
                py::overload_cast<const SharedTensor<double, 2>&, const SharedTensor<T, 2>&, Operation, Operation>(
                    &doublet<double, T>),
                "Returns the multiplication of two matrices A and B, with options to transpose/transpose-conjugate "
                "each beforehand",
                "A"_a, "B"_a, "opA"_a = Operation::None, "opB"_a = Operation::None,
                py::return_value_policy::reference_internal);
        mod.def("doublet",
                py::overload_cast<const SharedTensor<double, 2>&, const SharedTensor<T, 2>&, bool, bool>(
                    &doublet<double, T>),
                "Returns the multiplication of two matrices A and B, with options to transpose each beforehand", "A"_a,
                "B"_a, "transA"_a = false, "transB"_a = false, py::return_value_policy::reference_internal);
        // Type-homogeneous GEMV-s
        mod.def("gemv",
                [](const SharedTensor<T, 2>& A, const SharedTensor<T, 1>& x, Operation opA, double alpha, double beta) {
                    return gemv(opA, alpha, A, x, beta);
                },
                "A"_a, "x"_a, "opA"_a = Operation::None, "alpha"_a = 1.0, "beta"_a = 0.0);
        // Type-inhomogeneous GEMV-s
        // T x double
        mod.def("gemv",
                [](const SharedTensor<T, 2>& A, const SharedTensor<double, 1>& x, Operation opA, double alpha,
                   double beta) { return gemv(opA, alpha, A, x, beta); },
                "A"_a, "x"_a, "opA"_a = Operation::None, "alpha"_a = 1.0, "beta"_a = 0.0);
        // double x T
        mod.def("gemv",
                [](const SharedTensor<double, 2>& A, const SharedTensor<T, 1>& x, Operation opA, double alpha,
                   double beta) { return gemv(opA, alpha, A, x, beta); },
                "A"_a, "x"_a, "opA"_a = Operation::None, "alpha"_a = 1.0, "beta"_a = 0.0);
        // Decompositions
        mod.def("cholesky", py::overload_cast<const SharedTensor<T, 2>&>(&cholesky<T>),
                "Compute the Cholesky decomposition of A", "A"_a);
        mod.def("qr", py::overload_cast<const SharedTensor<T, 2>&, xt::linalg::qrmode>(&qr<T>),
                "Compute the QR decomposition of A", "A"_a, "mode"_a = xt::linalg::qrmode::reduced);
        mod.def("svd", py::overload_cast<const SharedTensor<T, 2>&, bool, bool>(&svd<T>),
                "Compute the singular value decomposition of A", "A"_a, "full_matrices"_a = true,
                "compute_uv"_a = true);
        // Matrix eigenvalues
        mod.def("eig", py::overload_cast<const SharedTensor<T, 2>&>(&eig<T>),
                "Compute the eigenvalues and right eigenvectors of a square matrix.", "A"_a);
        mod.def("eigvals", py::overload_cast<const SharedTensor<T, 2>&>(&eigvals<T>),
                "Compute the eigenvalues of a general matrix.", "A"_a);
        mod.def("eigh", py::overload_cast<const SharedTensor<T, 2>&, char>(&eigh<T>),
                "Compute eigenvalues and eigenvectors of a complex Hermitian (conjugate symmetric) or a real "
                "symmetric matrix.",
                "A"_a, "UPLO"_a = 'L');
        mod.def("eigvalsh", py::overload_cast<const SharedTensor<T, 2>&, char>(&eigvalsh<T>),
                "Compute the eigenvalues of a complex Hermitian or real symmetric matrix.", "A"_a, "UPLO"_a = 'L');
    }
};

template <typename T, typename U>
void declareRankN_builders(py::module& mod) {
    mod.def("full_like",
            py::overload_cast<const std::shared_ptr<T>&, U>(&full_like<typename T::value_type, T::rank, U>),
            "Returns a tensor with all blocks filled with given value of same shape and value type as input", "mold"_a,
            "fill_value"_a);
    mod.def("zeros_like", py::overload_cast<const std::shared_ptr<T>&>(&zeros_like<typename T::value_type, T::rank, U>),
            "Return a tensor with all blocks filled with 0 of same shape and value type as input", "mold"_a);
    mod.def("ones_like", py::overload_cast<const std::shared_ptr<T>&>(&ones_like<typename T::value_type, T::rank, U>),
            "Return a tensor with all blocks filled with 1 of same shape and value type as input", "mold"_a);
}

template <typename T, size_t Rank>
void bind_tensor(py::module& mod) {
    using Class = Tensor<T, Rank>;
    using SharedClass = SharedTensor<T, Rank>;
    using PyClass = py::class_<Class, SharedClass>;

    std::string name = Class::crtp_base::pyClassName();

    PyClass cls(mod, name.c_str());

    // CTORs shared by all ranks
    cls.def(py::init<const std::string&, size_t, const std::array<Dimension, Rank>&, unsigned int, T>(),
            ("Labeled, blocked, symmetry-assigned " + name).c_str(), "label"_a, "nirrep"_a, "axes_dimpi"_a,
            "symmetry"_a, "fill_value"_a = static_cast<T>(0));
    cls.def(py::init<const std::string&, const std::array<Dimension, Rank>&, T>(), ("Labeled, 1-irrep " + name).c_str(),
            "label"_a, "axes_dimpi"_a, "fill_value"_a = static_cast<T>(0));
    cls.def(py::init<size_t, const std::array<Dimension, Rank>&, unsigned int, T>(),
            ("Unlabeled, blocked, symmetry-assigned " + name).c_str(), "nirrep"_a, "axes_dimpi"_a, "symmetry"_a,
            "fill_value"_a = static_cast<T>(0));
    cls.def(py::init<size_t, const std::array<Dimension, Rank>&, T>(), ("Unlabeled, blocked " + name).c_str(),
            "nirrep"_a, "axes_dimpi"_a, "fill_value"_a = static_cast<T>(0));
    cls.def(py::init<const std::array<Dimension, Rank>&, T>(), ("Unlabeled, 1-irrep " + name).c_str(), "axes_dimpi"_a,
            "fill_value"_a = static_cast<T>(0));

    // Member functions shared by all ranks
    cls.def_property_readonly("dim", [](const Class& obj) { return obj.dim(); }, "Total number of elements");
    cls.def_property_readonly("nirrep", &Class::nirrep, "Number of irreps");
    cls.def_property("label", &Class::label, &Class::set_label, ("The label of the " + name).c_str());
    cls.def("axes_dimpi", py::overload_cast<>(&Class::axes_dimpi, py::const_), "Returns Dimension objects for all axes",
            py::return_value_policy::copy);
    cls.def("axes_dimpi", py::overload_cast<size_t>(&Class::axes_dimpi, py::const_),
            "Returns the Dimension object for given axis", "axis"_a);
    cls.def_property_readonly("shapes", [](const Class& obj) { return obj.shapes(); }, py ::return_value_policy::copy,
                              "Shapes of blocks");
    cls.def_property("symmetry", &Class::symmetry, &Class::set_symmetry, ("The symmetry of " + name).c_str());
    cls.def("__repr__", &Class::repr);
    cls.def("__str__", &Class::str);
    cls.def("__format__", &Class::format, "extra"_a = "");
    cls.def("__getitem__", py::overload_cast<size_t>(&Class::operator[]), "Return block at given irrep", "h"_a,
            py::is_operator(), py::return_value_policy::reference_internal);
    cls.def("__setitem__",
            [](Class& obj, size_t h, const xt::pytensor<T, Rank>& block) {
                if (h >= obj.nirrep()) throw py::index_error();
                obj[h] = block;
            },
            "h"_a, "block"_a, "Set block at given irrep", py::is_operator());
    cls.def("__len__", &Class::nirrep);
    cls.def("__iter__", [](const Class& obj) { return py::make_iterator(obj.begin(), obj.end()); },
            py::keep_alive<0, 1>() /* Essential: keep object alive while iterator exists */);

    // Free functions shared by all ranks
    // Builders
    declareRankN_builders<Class, int>(mod);
    declareRankN_builders<Class, float>(mod);
    declareRankN_builders<Class, double>(mod);
    declareRankN_builders<Class, std::complex<double>>(mod);
    mod.def("real", &real<T, Rank>, "Return real part of tensor", "in"_a);
    mod.def("imag", &imag<T, Rank>, "Return imaginary part of tensor. For real tensors, returns zeros.", "in"_a);
    mod.def("conj", &conj<T, Rank>, "Return complex conjugate of tensor", "in"_a);

    // Rank-dependent bindings, e.g. CTORs, member and free functions
    Decorator<T, Rank>::decorate(mod, cls);
}
}  // namespace

void export_linalg(py::module& mod) {
    xt::import_numpy();

    py::module sub_mod = mod.def_submodule("linalg");

    py::enum_<Operation>(sub_mod, "Operation")
        .value("none", Operation::None)
        .value("transpose", Operation::Transpose)
        .value("transpose_conj", Operation::TransposeConj);

    py::enum_<xt::linalg::qrmode>(sub_mod, "QRMode")
        .value("reduced", xt::linalg::qrmode::reduced)
        .value("complete", xt::linalg::qrmode::complete)
        .value("r", xt::linalg::qrmode::r)
        .value("raw", xt::linalg::qrmode::raw);

    // Rank-1 tensor, aka blocked vector
    bind_tensor<float, 1>(sub_mod);
    bind_tensor<double, 1>(sub_mod);
    bind_tensor<std::complex<double>, 1>(sub_mod);

    // Rank-2 tensor, aka blocked matrix
    bind_tensor<float, 2>(sub_mod);
    bind_tensor<double, 2>(sub_mod);
    bind_tensor<std::complex<double>, 2>(sub_mod);

    // Rank-3 tensor
    bind_tensor<float, 3>(sub_mod);
    bind_tensor<double, 3>(sub_mod);
    bind_tensor<std::complex<double>, 3>(sub_mod);
}
