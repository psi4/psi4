//
// Created by Justin Turney on 9/2/16.
//

#ifndef PSI4_CORE_PYBIND11_H_H
#define PSI4_CORE_PYBIND11_H_H

#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/operators.h>
#include <pybind11/eval.h>

namespace py =  pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

#endif //PSI4_CORE_PYBIND11_H_H
