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

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>

#include <functional>
#include <string>
#include <vector>

#define FORCE_IMPORT_ARRAY              // numpy C api loading
#include "xtensor-python/pytensor.hpp"  // Numpy bindings

#include "psi4/libmints/dimension.h"
#include "psi4/libmints/linalg.h"
#include "psi4/libmints/vector.h"
#include "psi4/libpsi4util/exception.h"

using namespace psi;
namespace py = pybind11;
using namespace pybind11::literals;

namespace {
struct Rank1Decorator final {
    template <typename PyClass>
    void operator()(PyClass& cls, const std::string& /* name */) const {
        cls.def(py::init<const std::string&, const Dimension&>(), "Labeled, blocked vector", "label"_a, "dimpi"_a);
        cls.def(py::init<const std::string&, int>(), "Labeled, 1-irrep vector", "label"_a, "dim"_a);
        cls.def(py::init<const Dimension&>(), "Unlabeled, blocked vector", "dimpi"_a);
        cls.def(py::init<int>(), "Unlabeled, 1-irrep vector", "dim"_a);
        cls.def_property_readonly("dimpi", [](const typename PyClass::type& obj) { return obj.dimpi(); },
                                  py::return_value_policy::copy, "Return the Dimension object");
    }
};

struct Rank2Decorator final {
    template <typename PyClass>
    void operator()(PyClass& cls, const std::string& /* name */) const {
        cls.def(py::init<const std::string&, const Dimension&, const Dimension&, unsigned int>(),
                "Labeled, blocked, symmetry-assigned matrix", "label"_a, "rowspi"_a, "colspi"_a, "symmetry"_a);
        cls.def(py::init<const std::string&, const Dimension&, const Dimension&>(), "Labeled, blocked matrix",
                "label"_a, "rowspi"_a, "colspi"_a);
        cls.def(py::init<const std::string&, int, int>(), "Labeled, 1-irrep matrix", "label"_a, "rows"_a, "cols"_a);
        cls.def(py::init<const Dimension&, const Dimension&, unsigned int>(),
                "Unlabeled, blocked, symmetry-assigned matrix", "rowspi"_a, "colspi"_a, "symmetry"_a);
        cls.def(py::init<const Dimension&, const Dimension&>(), "Unlabeled, blocked matrix", "rowspi"_a, "colspi"_a);
        cls.def(py::init<int, int>(), "Unlabeled, 1-irrep matrix", "rows"_a, "cols"_a);
        cls.def_property_readonly("rowspi", [](const typename PyClass::type& obj) { return obj.rowspi(); },
                                  py::return_value_policy::copy, "Returns the rows per irrep array");
        cls.def("rows", [](const typename PyClass::type& obj, size_t h) { return obj.rows(h); },
                "Returns the number of rows in given irrep", "h"_a = 0);
        cls.def_property_readonly("colspi", [](const typename PyClass::type& obj) { return obj.colspi(); },
                                  py::return_value_policy::copy, "Returns the columns per irrep array");
        cls.def("cols", [](const typename PyClass::type& obj, size_t h) { return obj.cols(h); },
                "Returns the number of columns in given irrep", "h"_a = 0);
    }
};

template <size_t Rank>
struct RankNDecorator final {
    template <typename PyClass>
    void operator()(PyClass& cls, const std::string& name) const {
        cls.def(py::init<const std::string&, size_t, const std::array<Dimension, Rank>&, unsigned int>(),
                ("Labeled, blocked, symmetry-assigned " + name).c_str(), "label"_a, "blocks"_a, "axes_dimpi"_a,
                "symmetry"_a);
        cls.def(py::init<const std::string&, const std::array<Dimension, Rank>&>(),
                ("Labeled, 1-irrep " + name).c_str(), "label"_a, "axes_dimpi"_a);
        cls.def(py::init<size_t, const std::array<Dimension, Rank>&, unsigned int>(),
                ("Unlabeled, blocked, symmetry-assigned " + name).c_str(), "blocks"_a, "axes_dimpi"_a, "symmetry"_a);
        cls.def(py::init<size_t, const std::array<Dimension, Rank>&>(), ("Unlabeled, blocked " + name).c_str(),
                "blocks"_a, "axes_dimpi"_a);
        cls.def(py::init<const std::array<Dimension, Rank>&>(), ("Unlabeled, 1-irrep " + name).c_str(), "axes_dimpi"_a);
    }
};

template <typename T, size_t Rank>
struct DeclareTensor final {
    using Class = Tensor<T, Rank>;
    using PyClass = py::class_<Class, std::shared_ptr<Class>>;
    using SpecialBinder = std::function<void(PyClass&, const std::string&)>;

    static void bind_tensor(py::module& mod, const SpecialBinder& decorate) {
        std::string name = Class::pyClassName();

        PyClass cls(mod, name.c_str());

        cls.def_property_readonly("dim", &Class::dim, "Total number of elements");
        cls.def_property_readonly("nirrep", &Class::nirrep, "Number of irreps");
        cls.def_property("label", &Class::label, &Class::set_label, ("The label of the " + name).c_str());
        cls.def("axes_dimpi", &Class::axes_dimpi, "Returns the Dimension object for given axis", "axis"_a);
        cls.def_property_readonly("shapes", [](const Class& obj) { return obj.shapes(); },
                                  py ::return_value_policy::copy, "Shapes of blocks");
        cls.def("nph", py::overload_cast<size_t>(&Class::block), "Return block at given irrep", "h"_a = 0,
                py::return_value_policy::reference_internal);
        cls.def_property("symmetry", &Class::symmetry, &Class::set_symmetry, ("The symmetry of " + name).c_str());

        cls.def("__repr__", &Class::repr);
        cls.def("__str__", &Class::str);
        cls.def("__format__", &Class::format, "extra"_a = "");

        // Rank-dependent bindings, e.g. CTORs
        decorate(cls, name);
    }
};
}  // namespace

void export_linalg(py::module& mod) {
    xt::import_numpy();
    // Rank-1 tensor, aka blocked vector
    auto decorate_v = Rank1Decorator();
    DeclareTensor<float, 1>::bind_tensor(mod, decorate_v);
    DeclareTensor<double, 1>::bind_tensor(mod, decorate_v);
    // Rank-2 tensor, aka blocked matrix
    auto decorate_m = Rank2Decorator();
    DeclareTensor<float, 2>::bind_tensor(mod, decorate_m);
    DeclareTensor<double, 2>::bind_tensor(mod, decorate_m);
    // Rank-3 tensor
    // FIXME not ideal to have rank info spread out...
    auto decorate_rankn = RankNDecorator<3>();
    DeclareTensor<float, 3>::bind_tensor(mod, decorate_rankn);
    DeclareTensor<double, 3>::bind_tensor(mod, decorate_rankn);

    using Class = Vector;
    using PyClass = py::class_<Vector, std::shared_ptr<Vector>>;

    PyClass cls(mod, "Vector", "Class for creating and manipulating vectors", py::dynamic_attr());
    cls.def(py::init<int>())
        .def(py::init<const Dimension&>())
        .def(py::init<const std::string&, int>())
        .def(py::init<const std::string&, const Dimension&>())
        .def_property("name", py::cpp_function(&Vector::name), py::cpp_function(&Vector::set_name),
                      "The name of the Vector. Used in printing.")
        .def("get", py::overload_cast<int>(&Vector::get, py::const_), "Returns a single element value located at m",
             "m"_a)

        .def("get", py::overload_cast<int, int>(&Vector::get, py::const_),
             "Returns a single element value located at m in irrep h", "h"_a, "m"_a)
        .def("set", py::overload_cast<int, double>(&Vector::set), "Sets a single element value located at m", "m"_a,
             "val"_a)
        .def("set", py::overload_cast<int, int, double>(&Vector::set),
             "Sets a single element value located at m in irrep h", "h"_a, "m"_a, "val"_a)
        .def("print_out", &Vector::print_out, "Prints the vector to the output file")
        .def("scale", &Vector::scale, "Scales the elements of a vector by sc", "sc"_a)
        .def("dim", &Vector::dim, "Returns the dimensions of the vector per irrep h", "h"_a = 0)
        .def("nirrep", &Vector::nirrep, "Returns the number of irreps")
        .def("get_block", &Vector::get_block, "Get a vector block", "slice"_a)
        .def("set_block", &Vector::set_block, "Set a vector block", "slice"_a, "block"_a)
        .def("dimpi", &Vector::dimpi, "Returns the Dimension object")
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
}
