#ifndef DFCCWRAPPER_H
#define DFCCWRAPPER_H

#include <libpsio/psio.hpp>
#include <libmints/wavefunction.h>
#include <psi4-dec.h>

namespace psi { namespace dfcc {

PsiReturnType dfcc(Options & options);

}}

#endif
