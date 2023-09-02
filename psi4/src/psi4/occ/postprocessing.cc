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

#include "psi4/libpsi4util/process.h"

#include "occwave.h"

namespace psi {
namespace occwave {

void OCCWave::mp2_printing(bool scf, bool include_singles, bool incomplete_singles) {
    outfile->Printf("\n");
    std::string which_mos = scf ? "SCF" : "optimized";
    std::string parenthetical = scf ? (include_singles ? " (ROHF-MP2)" : " (Canonical MP2)") : "";
    outfile->Printf("\tComputing MP2 energy using %s MOs%s... \n", which_mos.c_str(), parenthetical.c_str());
    outfile->Printf("\t============================================================================== \n");
    outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
    outfile->Printf("\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
    outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
    outfile->Printf("\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp2AA);
    outfile->Printf("\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp2AB);
    outfile->Printf("\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp2BB);
    outfile->Printf("\tScaled_SS Correlation Energy (a.u.): %20.14f\n", Escsmp2AA + Escsmp2BB);
    outfile->Printf("\tScaled_OS Correlation Energy (a.u.): %20.14f\n", Escsmp2AB);
    if (incomplete_singles) {
        outfile->Printf("\tMP2 Total Energy (a.u.)            : %s\n", "singles and total energies not available");
    } else {
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
    }
    outfile->Printf("\t============================================================================== \n");
    outfile->Printf("\n");
}

void OCCWave::mp2p5_printing(bool scf, bool include_singles, bool incomplete_singles) {
    outfile->Printf("\n");
    std::string which_mos = scf ? "SCF" : "optimized";
    std::string parenthetical = scf ? " (Canonical MP2.5)" : "";
    outfile->Printf("\tComputing MP2.5 energy using %s MOs%s... \n", which_mos.c_str(), parenthetical.c_str());
    outfile->Printf("\t============================================================================== \n");
    outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
    outfile->Printf("\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
    outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
    outfile->Printf("\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp3AA);
    outfile->Printf("\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp3AB);
    outfile->Printf("\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp3BB);
    if (incomplete_singles) {
        outfile->Printf("\tMP2.5 Total Energy (a.u.)            : %s\n", "singles and total energies not available");
    } else {
        outfile->Printf("\t0.5 Energy Correction (a.u.)       : %20.14f\n", Emp3 - Emp2);
        outfile->Printf("\tMP2.5 Correlation Energy (a.u.)    : %20.14f\n", Ecorr);
        outfile->Printf("\tMP2.5 Total Energy (a.u.)          : %20.14f\n", Emp3);
    }
    outfile->Printf("\t============================================================================== \n");
    outfile->Printf("\n");
}

void OCCWave::mp3_printing(bool scf, bool include_singles, bool incomplete_singles) {
    outfile->Printf("\n");
    std::string which_mos = scf ? "SCF" : "optimized";
    std::string parenthetical = scf ? " (Canonical MP3)" : "";
    outfile->Printf("\tComputing MP3 energy using %s MOs%s... \n", which_mos.c_str(), parenthetical.c_str());
    outfile->Printf("\t============================================================================== \n");
    outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
    outfile->Printf("\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
    outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
    outfile->Printf("\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp3AA);
    outfile->Printf("\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp3AB);
    outfile->Printf("\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp3BB);
    if (incomplete_singles) {
        outfile->Printf("\tMP3 Total Energy (a.u.)            : %s\n", "singles and total energies not available");
    } else {
        outfile->Printf("\tMP2.5 Correlation Energy (a.u.)    : %20.14f\n", (Emp2 - Escf) + 0.5 * (Emp3 - Emp2));
        outfile->Printf("\tMP2.5 Total Energy (a.u.)          : %20.14f\n", 0.5 * (Emp3 + Emp2));
        outfile->Printf("\tSCS-MP3 Total Energy (a.u.)        : %20.14f\n", Escsmp3);
        outfile->Printf("\t3rd Order Energy (a.u.)            : %20.14f\n", Emp3 - Emp2);
        outfile->Printf("\tMP3 Correlation Energy (a.u.)      : %20.14f\n", Ecorr);
        outfile->Printf("\tMP3 Total Energy (a.u.)            : %20.14f\n", Emp3);
    }
    outfile->Printf("\t============================================================================== \n");
    outfile->Printf("\n");
}

void OCCWave::mp2_postprocessing(bool include_singles, bool incomplete_singles) {
    mp2_printing(true, include_singles, incomplete_singles);

    if (incomplete_singles) {
        variables_["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp2AB;
        variables_["MP2 SAME-SPIN CORRELATION ENERGY"] = Emp2AA + Emp2BB;
        variables_["MP2 DOUBLES ENERGY"] = Emp2AB + Emp2AA + Emp2BB;

        return;
    }

    // Wavefunctions saved to variables are set on the wavefunction Py-side.
    // Other variables will be set Py-side via QCDB formulas, per a future Lori project.
    variables_["MP2 CORRELATION ENERGY"] = Emp2 - Escf;
    variables_["SCS-MP2 CORRELATION ENERGY"] = Escsmp2 - Escf;
    variables_["CUSTOM SCS-MP2 CORRELATION ENERGY"] = os_scale * Emp2AB + ss_scale * (Emp2AA + Emp2BB) + Emp2_t1;
    variables_["SOS-MP2 CORRELATION ENERGY"] = Esosmp2 - Escf;
    variables_["SCS(N)-MP2 CORRELATION ENERGY"] = Escsnmp2 - Escf;
    variables_["SCS-MP2-VDW CORRELATION ENERGY"] = Escsmp2vdw - Escf;
    variables_["SOS-PI-MP2 CORRELATION ENERGY"] = Esospimp2 - Escf;

    variables_["MP2 TOTAL ENERGY"] = Emp2;
    variables_["SCS-MP2 TOTAL ENERGY"] = Escf + variables_["SCS-MP2 CORRELATION ENERGY"];
    variables_["CUSTOM SCS-MP2 TOTAL ENERGY"] = Escf + variables_["CUSTOM SCS-MP2 CORRELATION ENERGY"];
    variables_["SOS-MP2 TOTAL ENERGY"] = Esosmp2;
    variables_["SCS(N)-MP2 TOTAL ENERGY"] = Escsnmp2;
    variables_["SCS-MP2-VDW TOTAL ENERGY"] = Escsmp2vdw;
    variables_["SOS-PI-MP2 TOTAL ENERGY"] = Esospimp2;

    variables_["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp2AB;
    variables_["MP2 SAME-SPIN CORRELATION ENERGY"] = Emp2AA + Emp2BB;
    variables_["MP2 SINGLES ENERGY"] = Emp2_t1;
    variables_["MP2 DOUBLES ENERGY"] = Emp2AB + Emp2AA + Emp2BB;
}

void OCCWave::mp3_postprocessing(bool include_singles, bool incomplete_singles) {
    mp3_printing(true, include_singles, incomplete_singles);

    if (incomplete_singles) {
        // by analogy to MP2, the vars below might be correct for ROHF-MP3 in the course of OOMP3,
        //   but since no reference values to confirm, compromise by printing but not saving.
        // variables_["MP3 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp3AB;
        // variables_["MP3 SAME-SPIN CORRELATION ENERGY"] = Emp3AA + Emp3BB;
        // variables_["MP3 DOUBLES ENERGY"] = Emp3AB + Emp3AA + Emp3BB;

        return;
    }

    variables_["MP3 TOTAL ENERGY"] = Emp3;
    variables_["SCS-MP3 TOTAL ENERGY"] = Escsmp3;

    variables_["MP2.5 CORRELATION ENERGY"] = (Emp2 - Escf) + 0.5 * (Emp3 - Emp2);
    variables_["MP2.5 TOTAL ENERGY"] = 0.5 * (Emp3 + Emp2);
    variables_["MP3 CORRELATION ENERGY"] = Emp3 - Escf;
    variables_["SCS-MP3 CORRELATION ENERGY"] = Escsmp3 - Escf;

    variables_["CUSTOM SCS-MP3 CORRELATION ENERGY"] = os_scale * Emp3AB + ss_scale * (Emp3AA + Emp3BB);
    variables_["CUSTOM SCS-MP3 TOTAL ENERGY"] = Escf + variables_["CUSTOM SCS-MP3 CORRELATION ENERGY"];

    variables_["MP3 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp3AB;
    variables_["MP3 SAME-SPIN CORRELATION ENERGY"] = Emp3AA + Emp3BB;
    variables_["MP3 SINGLES ENERGY"] = 0.0;  // RHF & UHF only
    variables_["MP3 DOUBLES ENERGY"] = Emp3AB + Emp3AA + Emp3BB;

    variables_["MP2.5 OPPOSITE-SPIN CORRELATION ENERGY"] = 0.5 * (Emp2AB + Emp3AB);
    variables_["MP2.5 SAME-SPIN CORRELATION ENERGY"] = 0.5 * (Emp2AA + Emp2BB + Emp3AA + Emp3BB);
    variables_["MP2.5 SINGLES ENERGY"] = 0.0;  // RHF & UHF only
    variables_["MP2.5 DOUBLES ENERGY"] = 0.5 * (Emp2AB + Emp2AA + Emp2BB + Emp3AB + Emp3AA + Emp3BB);
}

void OCCWave::mp2p5_postprocessing(bool include_singles, bool incomplete_singles) {
    mp2p5_printing(true, include_singles, incomplete_singles);

    if (incomplete_singles) {
        // by analogy to MP2, the vars below might be correct for ROHF-MP3 in the course of OOMP3,
        //   but since no reference values to confirm, compromise by printing but not saving.
        // variables_["MP3 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp3AB;
        // variables_["MP3 SAME-SPIN CORRELATION ENERGY"] = Emp3AA + Emp3BB;
        // variables_["MP3 DOUBLES ENERGY"] = Emp3AB + Emp3AA + Emp3BB;

        return;
    }

    variables_["MP2.5 TOTAL ENERGY"] = Emp3;
    variables_["MP2.5 CORRELATION ENERGY"] = Emp3 - Escf;
    variables_["MP3 TOTAL ENERGY"] = Emp2 + 2.0 * (Emp3 - Emp2);
    variables_["MP3 CORRELATION ENERGY"] = Emp2 + 2.0 * (Emp3 - Emp2) - Escf;

    variables_["CUSTOM SCS-MP2.5 CORRELATION ENERGY"] = os_scale * Emp3AB + ss_scale * (Emp3AA + Emp3BB);
    variables_["CUSTOM SCS-MP2.5 TOTAL ENERGY"] = Escf + variables_["CUSTOM SCS-MP2.5 CORRELATION ENERGY"];

    variables_["MP2.5 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp3AB;
    variables_["MP2.5 SAME-SPIN CORRELATION ENERGY"] = Emp3AA + Emp3BB;
    variables_["MP2.5 SINGLES ENERGY"] = 0.0;  // RHF & UHF only
    variables_["MP2.5 DOUBLES ENERGY"] = Emp3 - Escf;  // RHF & UHF only

    variables_["MP3 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp2AB + 2.0 * (Emp3AB - Emp2AB);
    variables_["MP3 SAME-SPIN CORRELATION ENERGY"] = Emp2AA + Emp2BB + 2.0 * (Emp3AA + Emp3BB - Emp2AA - Emp2BB);
    variables_["MP3 SINGLES ENERGY"] = 0.0;  // RHF & UHF only
    variables_["MP3 DOUBLES ENERGY"] = variables_["MP3 CORRELATION ENERGY"];  // RHF & UHF only
}
}
}
