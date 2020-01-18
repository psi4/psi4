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

#include "psi4/libpsi4util/process.h"

#include "occwave.h"

namespace psi {
namespace occwave {

void OCCWave::mp2_postprocessing(bool include_singles) {
    outfile->Printf("\n");
    outfile->Printf("\tComputing MP2 energy using SCF MOs (%sMP2)... \n", include_singles ? "ROHF-" : "Canonical ");
    outfile->Printf("\t============================================================================== \n");
    outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
    outfile->Printf("\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
    outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
    outfile->Printf("\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp2AA);
    outfile->Printf("\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp2AB);
    outfile->Printf("\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp2BB);
    outfile->Printf("\tScaled_SS Correlation Energy (a.u.): %20.14f\n", Escsmp2AA + Escsmp2BB);
    outfile->Printf("\tScaled_OS Correlation Energy (a.u.): %20.14f\n", Escsmp2AB);
    outfile->Printf("\tSCS-MP2 Total Energy (a.u.)        : %20.14f\n", Escsmp2);
    outfile->Printf("\tSOS-MP2 Total Energy (a.u.)        : %20.14f\n", Esosmp2);
    outfile->Printf("\tSCSN-MP2 Total Energy (a.u.)       : %20.14f\n", Escsnmp2);
    outfile->Printf("\tSCS-MP2-VDW Total Energy (a.u.)    : %20.14f\n", Escsmp2vdw);
    outfile->Printf("\tSOS-PI-MP2 Total Energy (a.u.)     : %20.14f\n", Esospimp2);
    if (include_singles) {
        outfile->Printf("\tMP2 Singles Energy (a.u.)          : %20.14f\n", Emp2_t1);
        outfile->Printf("\tMP2 Doubles Energy (a.u.)          : %20.14f\n", Ecorr - Emp2_t1);
    }
    outfile->Printf("\tMP2 Correlation Energy (a.u.)      : %20.14f\n", Ecorr);
    outfile->Printf("\tMP2 Total Energy (a.u.)            : %20.14f\n", Emp2);
    outfile->Printf("\t============================================================================== \n");
    outfile->Printf("\n");

    // Wavefunctions saved to variables are set on the wavefunction Py-side.
    // Other variables will be set Py-side via QCDB formulas, per a future Lori project.
    variables_["MP2 CORRELATION ENERGY"] = Emp2 - Escf;
    variables_["SCS-MP2 CORRELATION ENERGY"] = 6.0/5.0 * Emp2AB + 1.0/3.0 * (Emp2AA + Emp2BB) + Emp2_t1;
    variables_["CUSTOM SCS-MP2 CORRELATION ENERGY"] = os_scale * Emp2AB + ss_scale * (Emp2AA + Emp2BB) + Emp2_t1;
    Process::environment.globals["SOS-MP2 CORRELATION ENERGY"] = Esosmp2 - Escf;
    Process::environment.globals["SCSN-MP2 CORRELATION ENERGY"] = Escsnmp2 - Escf;
    Process::environment.globals["SCS-MP2-VDW CORRELATION ENERGY"] = Escsmp2vdw - Escf;
    Process::environment.globals["SOS-PI-MP2 CORRELATION ENERGY"] = Esospimp2 - Escf;

    variables_["MP2 TOTAL ENERGY"] = Emp2;
    variables_["SCS-MP2 TOTAL ENERGY"] =  Escf + variables_["SCS-MP2 CORRELATION ENERGY"];
    variables_["CUSTOM SCS-MP2 TOTAL ENERGY"] =  Escf + variables_["CUSTOM SCS-MP2 CORRELATION ENERGY"];
    Process::environment.globals["SOS-MP2 TOTAL ENERGY"] = Esosmp2;
    Process::environment.globals["SCSN-MP2 TOTAL ENERGY"] = Escsnmp2;
    Process::environment.globals["SCS-MP2-VDW TOTAL ENERGY"] = Escsmp2vdw;
    Process::environment.globals["SOS-PI-MP2 TOTAL ENERGY"] = Esospimp2;


    variables_["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp2AB;
    variables_["MP2 SAME-SPIN CORRELATION ENERGY"] = Emp2AA + Emp2BB;
}
}
}
