#ifndef SAPTENERGY_H
#define SAPTENERGY_H

#include <libpsio/psio.hpp>
#include <libmints/wavefunction.h>
#include <psi4-dec.h>

namespace psi { namespace sapt {

PsiReturnType sapt(Options & options);

}}

#endif
