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

#include "jk_grad.h"

#include "psi4/libmints/mintshelper.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsi4util/process.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace psi;

namespace psi {
namespace scfgrad {

JKGrad::JKGrad(int deriv, std::shared_ptr<BasisSet> primary) :
    deriv_(deriv), primary_(primary)
{
    common_init();
}
JKGrad::~JKGrad()
{
}
std::shared_ptr<JKGrad> JKGrad::build_JKGrad(int deriv, std::shared_ptr<MintsHelper> mints)
{
    Options& options = Process::environment.options;

    if (options.get_str("SCF_TYPE").find("DF") != std::string::npos) {
        DFJKGrad* jk = new DFJKGrad(deriv, mints);

        if (options["INTS_TOLERANCE"].has_changed())
            jk->set_cutoff(options.get_double("INTS_TOLERANCE"));
        if (options["PRINT"].has_changed())
            jk->set_print(options.get_int("PRINT"));
        if (options["DEBUG"].has_changed())
            jk->set_debug(options.get_int("DEBUG"));
        if (options["BENCH"].has_changed())
            jk->set_bench(options.get_int("BENCH"));
        jk->set_condition(options.get_double("DF_FITTING_CONDITION"));
        if (options["DF_INTS_NUM_THREADS"].has_changed())
            jk->set_df_ints_num_threads(options.get_int("DF_INTS_NUM_THREADS"));

        return std::shared_ptr<JKGrad>(jk);
    } else if (options.get_str("SCF_TYPE") == "DIRECT" || options.get_str("SCF_TYPE") == "PK" || options.get_str("SCF_TYPE") == "OUT_OF_CORE") {

        DirectJKGrad* jk = new DirectJKGrad(deriv, mints->get_basisset("ORBITAL"));

        if (options["INTS_TOLERANCE"].has_changed())
            jk->set_cutoff(options.get_double("INTS_TOLERANCE"));
        if (options["PRINT"].has_changed())
            jk->set_print(options.get_int("PRINT"));
        if (options["DEBUG"].has_changed())
            jk->set_debug(options.get_int("DEBUG"));
        if (options["BENCH"].has_changed())
            jk->set_bench(options.get_int("BENCH"));
        // TODO: rename every DF case
        if (options["DF_INTS_NUM_THREADS"].has_changed())
            jk->set_ints_num_threads(options.get_int("DF_INTS_NUM_THREADS"));

        return std::shared_ptr<JKGrad>(jk);

    } else {
        throw PSIEXCEPTION("JKGrad::build_JKGrad: Unknown SCF Type");
    }
}
void JKGrad::common_init() {
    print_ = 1;
    debug_ = 0;
    bench_ = 0;

    memory_ = 32000000L;
    omp_num_threads_ = 1;
#ifdef _OPENMP
    omp_num_threads_ = Process::environment.get_n_threads();
#endif

    cutoff_ = Process::environment.options.get_double("INTS_TOLERANCE");

    do_J_ = true;
    do_K_ = true;
    do_wK_ = false;
    omega_ = 0.0;
}
}  // namespace scfgrad
}  // namespace psi
