//
// Created by Justin Turney on 1/5/16.
//

#ifndef AMBIT_CONVERTER_H
#define AMBIT_CONVERTER_H

#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"

namespace ambit
{

class Tensor;

namespace helpers
{

namespace psi4
{

void convert(const psi::Matrix &matrix, ambit::Tensor *target);

void convert(const psi::Vector &vector, ambit::Tensor *target);

} // namespace psi4

} // namespace helpers

} // namespace ambit

#endif // AMBIT_CONVERTER_H
