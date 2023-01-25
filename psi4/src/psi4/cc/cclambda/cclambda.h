/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2023 The Psi4 Developers.
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

#ifndef CCLAMBDA_H
#define CCLAMBDA_H

#include "psi4/libmints/wavefunction.h"
#include "psi4/cc/ccwave.h"

namespace psi {
class Wavefunction;
class Options;
}  // namespace psi

namespace psi {
namespace cclambda {

class CCLambdaWavefunction final : public psi::ccenergy::CCEnergyWavefunction {
   public:
    CCLambdaWavefunction(std::shared_ptr<Wavefunction> reference_wavefunction, Options &options);
    ~CCLambdaWavefunction() override;

    double compute_energy() override;

   private:
    void init();
    void init_io();
    void init_amps(const struct L_Params&);
    int **cacheprep_uhf(int level, int *cachefiles);
    int **cacheprep_rhf(int level, int *cachefiles);
    void cachedone_rhf(int **cachelist);
    void cachedone_uhf(int **cachelist);
    void cleanup();
    void denom(const struct L_Params&);
    void get_params(psi::Options &);
    void local_init();
    void local_done();
    void exit_io();
    void title();
    void get_moinfo(std::shared_ptr<psi::Wavefunction> wfn);

    int converged(int);
    void diis(int, int);
    void sort_amps(int);
    void status(const char *, std::string);
    void update();

    void cc2_L2_build(const struct L_Params&);
    void L2_build(const struct L_Params&);
};

}  // namespace cclambda
}  // namespace psi

#endif  // CCLAMBDA_H
