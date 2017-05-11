//
// Created by Justin Turney on 12/17/15.
//

#ifndef AMBIT_INTEGRALS_H
#define AMBIT_INTEGRALS_H

#include "psi4/libmints/onebody.h"
#include "psi4/libmints/twobody.h"

namespace ambit
{

class Tensor;

namespace helpers
{

namespace psi4
{

void integrals(psi::OneBodyAOInt &integral, ambit::Tensor *target);

void integrals(psi::TwoBodyAOInt &integral, ambit::Tensor *target);

} // namespace psi4

} // namespace helpers

} // namespace ambit

#endif // AMBIT_INTEGRALS_H
