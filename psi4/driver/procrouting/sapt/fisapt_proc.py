#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2025 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This file is part of Psi4.
#
# Psi4 is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3.
#
# Psi4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with Psi4; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

import os

import numpy as np

from psi4 import core

from .. import empirical_dispersion
from ...p4util.exceptions import ValidationError


def fisapt_compute_energy(self, jk_obj, *, external_potentials=None):
    """Computes the FSAPT energy. FISAPT::compute_energy"""

    # => Header <=

    self.print_header()

    # => Zero-th Order Wavefunction <=

    core.timer_on("FISAPT: Setup")
    self.localize()
    self.partition()
    self.overlap()
    self.kinetic()
    self.nuclear()
    self.coulomb(jk_obj)
    core.timer_off("FISAPT: Setup")
    core.timer_on("FISAPT: Monomer SCF")
    self.scf()
    core.timer_off("FISAPT: Monomer SCF")
    self.freeze_core()
    self.unify()
    self.nuclear()
    self.unify_part2()
    self.freeze_core()
    self.do_cubes()
    core.timer_on("FISAPT: Subsys E")
    self.dHF()
    core.timer_off("FISAPT: Subsys E")

    # => SAPT0 <=

    core.timer_on("FISAPT:SAPT:elst")
    self.elst()
    core.timer_off("FISAPT:SAPT:elst")
    core.timer_on("FISAPT:SAPT:exch")
    self.exch()
    core.timer_off("FISAPT:SAPT:exch")
    core.timer_on("FISAPT:SAPT:ind")
    self.ind()
    core.timer_off("FISAPT:SAPT:ind")
    if not core.get_option("FISAPT", "FISAPT_DO_FSAPT"):
        core.timer_on("FISAPT:SAPT:disp")
        # Expensive, only do if needed  # unteseted translation of below
        self.disp(self.matrices(), self.vectors(), True)
        # self.disp(matrices_, vectors_, true)  # Expensive, only do if needed
        core.timer_off("FISAPT:SAPT:disp")

    # => F-SAPT0 <=

    if core.get_option("FISAPT", "FISAPT_DO_FSAPT"):
        core.timer_on("FISAPT:FSAPT:loc")
        self.flocalize()
        core.timer_off("FISAPT:FSAPT:loc")
        core.timer_on("FISAPT:FSAPT:elst")
        self.felst()
        core.timer_off("FISAPT:FSAPT:elst")
        core.timer_on("FISAPT:FSAPT:exch")
        self.fexch()
        core.timer_off("FISAPT:FSAPT:exch")
        core.timer_on("FISAPT:FSAPT:ind")
        self.find()
        core.timer_off("FISAPT:FSAPT:ind")
        if core.get_option("FISAPT", "FISAPT_DO_FSAPT_DISP"):
            core.timer_on("FISAPT:FSAPT:disp")
            self.fdisp()
            core.timer_off("FISAPT:FSAPT:disp")
        # else:
        #    # Build Empirical Dispersion
        #    dashD = empirical_dispersion.EmpiricalDispersion(name_hint='SAPT0-D3M')
        #    dashD.print_out()
        #    # Compute -D
        #    Edisp = dashD.compute_energy(core.get_active_molecule())
        #    core.set_variable(
        # '{} DISPERSION CORRECTION ENERGY'.format(dashD.fctldash), Edisp)
        # Printing
        #    text = []
        #    text.append("   => {}: Empirical Dispersion <=".format(dashD.fctldash.upper()))
        #    text.append(" ")
        #    text.append(dashD.description)
        #    text.append(dashD.dashlevel_citation.rstrip())
        #    text.append("\n    Empirical Dispersion Energy [Eh] =     {:24.16f}\n".format(Edisp))
        #    text.append('\n')
        #    core.print_out('\n'.join(text))
        self.fdrop(external_potentials)

    # => Scalar-Field Analysis <=

    if core.get_option("FISAPT", "FISAPT_DO_PLOT"):
        core.timer_on("FISAPT:FSAPT:cubeplot")
        self.plot()
        core.timer_off("FISAPT:FSAPT:cubeplot")

    # => Summary <=

    self.print_trailer()


def fisapt_fdrop(self, external_potentials=None):
    """Drop output files from FSAPT calculation. FISAPT::fdrop"""

    core.print_out("  ==> F-SAPT Output <==\n\n")
    write_output_files = core.get_option("FISAPT", "FISAPT_FSAPT_FILEPATH").lower() != "none"

    if write_output_files:
        filepath = core.get_option("FISAPT", "FISAPT_FSAPT_FILEPATH")
        os.makedirs(filepath, exist_ok=True)

        core.print_out("    F-SAPT Data Filepath = {}\n\n".format(filepath))

        geomfile = filepath + os.sep + "geom.xyz"
        xyz = self.molecule().to_string(dtype="xyz", units="Angstrom")
        with open(geomfile, "w") as fh:
            fh.write(xyz)

    # write external potential geometries
    if external_potentials is not None and isinstance(external_potentials, dict):
        for frag in "ABC":
            potential = external_potentials.get(frag, None)
            if potential is not None:
                xyz = str(len(potential)) + "\n\n"
                potential_lst = []
                for qxyz in potential:
                    if len(qxyz) == 2:
                        xyz += "Ch %f %f %f\n" % (qxyz[1][0], qxyz[1][1], qxyz[1][2])
                        potential_lst.append(qxyz[1])
                    elif len(qxyz) == 4:
                        xyz += "Ch %f %f %f\n" % (qxyz[1], qxyz[2], qxyz[3])
                        potential_lst.append(qxyz[1:])
                    else:
                        raise ValidationError(
                            f"Point charge '{qxyz}' not mapping into 'chg, [x, y, z]' or 'chg, x, y, z'"
                        )
                potential_lst = np.array(potential_lst)
                core.set_variable("FSAPT_EXTERN_POTENTIAL_{}".format(frag), potential_lst)
                if write_output_files:
                    with open(filepath + os.sep + "Extern_%s.xyz" % frag, "w") as fh:
                        fh.write(xyz)

    vectors = self.vectors()
    matrices = self.matrices()

    matrices["Qocc0A"].name = "QA"
    matrices["Qocc0B"].name = "QB"
    matrices["Elst_AB"].name = "Elst"
    matrices["Exch_AB"].name = "Exch"
    matrices["IndAB_AB"].name = "IndAB"
    matrices["IndBA_AB"].name = "IndBA"
    core.set_variable("FSAPT_QA", matrices["Qocc0A"])
    core.set_variable("FSAPT_QB", matrices["Qocc0B"])
    core.set_variable("FSAPT_ELST_AB", matrices["Elst_AB"])
    core.set_variable("FSAPT_AB_SIZE", np.array(matrices["Elst_AB"].np.shape).reshape(1, -1))
    core.set_variable("FSAPT_EXCH_AB", matrices["Exch_AB"])
    core.set_variable("FSAPT_INDAB_AB", matrices["IndAB_AB"])
    core.set_variable("FSAPT_INDBA_AB", matrices["IndBA_AB"])

    if write_output_files:
        _drop(vectors["ZA"], filepath)
        _drop(vectors["ZB"], filepath)
        _drop(matrices["Qocc0A"], filepath)
        _drop(matrices["Qocc0B"], filepath)
        _drop(matrices["Elst_AB"], filepath)
        _drop(matrices["Exch_AB"], filepath)
        _drop(matrices["IndAB_AB"], filepath)
        _drop(matrices["IndBA_AB"], filepath)

    if core.get_option("FISAPT", "FISAPT_DO_FSAPT_DISP"):
        # In SAPT(DFT) case, you might not do dispersion
        if "Disp_AB" not in matrices:
            matrices["Disp_AB"] = core.Matrix.from_array(np.zeros_like(matrices["Elst_AB"]))
        matrices["Disp_AB"].name = "Disp"
        core.set_variable("FSAPT_DISP_AB", matrices["Disp_AB"])
        if write_output_files:
            _drop(matrices["Disp_AB"], filepath)

    if core.get_option("FISAPT", "SSAPT0_SCALE"):
        # NOTE: do same as above for conditionally writing
        ssapt_filepath = core.get_option("FISAPT", "FISAPT_FSSAPT_FILEPATH")
        write_ssapt_files = ssapt_filepath.lower() != "none"

        if write_ssapt_files:
            os.makedirs(ssapt_filepath, exist_ok=True)
            core.print_out("    F-sSAPT Data Filepath = {}\n\n".format(ssapt_filepath))
            geomfile = ssapt_filepath + os.sep + "geom.xyz"
            with open(geomfile, "w") as fh:
                fh.write(xyz)

        matrices["sIndAB_AB"].name = "IndAB"
        matrices["sIndBA_AB"].name = "IndBA"
        core.set_variable("FSAPT_SINDAB_AB", matrices["sIndAB_AB"])
        core.set_variable("FSAPT_SINDBA_AB", matrices["sIndBA_AB"])

        if write_ssapt_files:
            _drop(vectors["ZA"], ssapt_filepath)
            _drop(vectors["ZB"], ssapt_filepath)
            _drop(matrices["Qocc0A"], ssapt_filepath)
            _drop(matrices["Qocc0B"], ssapt_filepath)
            _drop(matrices["Elst_AB"], ssapt_filepath)
            _drop(matrices["Exch_AB"], ssapt_filepath)
            _drop(matrices["sIndAB_AB"], ssapt_filepath)
            _drop(matrices["sIndBA_AB"], ssapt_filepath)

        if core.get_option("FISAPT", "FISAPT_DO_FSAPT_DISP"):
            matrices["sDisp_AB"].name = "Disp"
            core.set_variable("FSAPT_SDISP_AB", matrices["sDisp_AB"])
            if write_ssapt_files:
                _drop(matrices["sDisp_AB"], ssapt_filepath)


def fisapt_variables_to_wfn(self, ref_wfn, external_potentials=None):
    """
    Stores FISAPT variables to the wavefunction for AtomicResults to
    store results.
    """
    # First Scalars
    scalars = self.scalars()
    ref_wfn.set_variable("SAPT ELST ENERGY", scalars["Electrostatics"])
    ref_wfn.set_variable("SAPT ELST10,R ENERGY", scalars["Elst10,r"])
    if "Extern-Extern" in scalars:
        ref_wfn.set_variable("SAPT ELST EXTERN-EXTERN ENERGY", scalars["Extern-Extern"])
    if core.has_variable("FSAPT_EXTERN_POTENTIAL_A"):
        ref_wfn.set_variable("FSAPT_EXTERN_POTENTIAL_A", core.variable("FSAPT_EXTERN_POTENTIAL_A"))
    if core.has_variable("FSAPT_EXTERN_POTENTIAL_B"):
        ref_wfn.set_variable("FSAPT_EXTERN_POTENTIAL_B", core.variable("FSAPT_EXTERN_POTENTIAL_B"))
    if core.has_variable("FSAPT_EXTERN_POTENTIAL_C"):
        ref_wfn.set_variable("FSAPT_EXTERN_POTENTIAL_C", core.variable("FSAPT_EXTERN_POTENTIAL_C"))
    ref_wfn.set_variable("SAPT EXCH ENERGY", scalars["Exchange"])
    ref_wfn.set_variable("SAPT EXCH10 ENERGY", scalars["Exch10"])
    ref_wfn.set_variable("SAPT EXCH10(S^2) ENERGY", scalars["Exch10(S^2)"])
    ref_wfn.set_variable("SAPT IND ENERGY", scalars["Induction"])
    ref_wfn.set_variable("SAPT IND20,R ENERGY", scalars["Ind20,r"])
    ref_wfn.set_variable("SAPT EXCH-IND20,R ENERGY", scalars["Exch-Ind20,r"])
    ref_wfn.set_variable("SAPT IND20,U ENERGY", scalars["Ind20,u"])
    ref_wfn.set_variable("SAPT EXCH-IND20,U ENERGY", scalars["Exch-Ind20,u"])
    ref_wfn.set_variable("SAPT DISP ENERGY", scalars["Dispersion"])
    ref_wfn.set_variable("SAPT DISP20 ENERGY", scalars["Disp20"])
    ref_wfn.set_variable("SAPT EXCH-DISP20 ENERGY", scalars["Exch-Disp20"])
    ref_wfn.set_variable("SAPT0 TOTAL ENERGY", scalars["SAPT"])
    ref_wfn.set_variable("SAPT TOTAL ENERGY", scalars["SAPT"])
    ref_wfn.set_variable("CURRENT ENERGY", scalars["SAPT"])
    # dHF to ref_wfn
    ref_wfn.set_variable("SAPT HF(2) ENERGY ABC(HF)", scalars["E_ABC_HF"])
    ref_wfn.set_variable("SAPT HF(2) ENERGY AC(0)", scalars["E_AC"])
    ref_wfn.set_variable("SAPT HF(2) ENERGY BC(0)", scalars["E_BC"])
    ref_wfn.set_variable("SAPT HF(2) ENERGY A(0)", scalars["E_A"])
    ref_wfn.set_variable("SAPT HF(2) ENERGY B(0)", scalars["E_B"])
    ref_wfn.set_variable("SAPT HF(2) ENERGY AC(HF)", scalars["E_AC_HF"])
    ref_wfn.set_variable("SAPT HF(2) ENERGY BC(HF)", scalars["E_BC_HF"])
    ref_wfn.set_variable("SAPT HF(2) ENERGY AB(HF)", scalars["E_AB_HF"])
    ref_wfn.set_variable("SAPT HF(2) ENERGY A(HF)", scalars["E_A_HF"])
    ref_wfn.set_variable("SAPT HF(2) ENERGY B(HF)", scalars["E_B_HF"])
    ref_wfn.set_variable("SAPT HF(2) ENERGY C", scalars["E_C"])
    ref_wfn.set_variable("SAPT HF(2) ENERGY HF", scalars["HF"])

    # Then matrices
    matrices = self.matrices()
    ref_wfn.set_variable("FSAPT_QA", matrices["Qocc0A"])
    ref_wfn.set_variable("FSAPT_QB", matrices["Qocc0B"])
    ref_wfn.set_variable("FSAPT_ELST_AB", matrices["Elst_AB"])
    ref_wfn.set_variable("FSAPT_EXCH_AB", matrices["Exch_AB"])
    ref_wfn.set_variable("FSAPT_INDAB_AB", matrices["IndAB_AB"])
    ref_wfn.set_variable("FSAPT_INDBA_AB", matrices["IndBA_AB"])

    # Handle conditional cases
    if core.get_option("FISAPT", "FISAPT_DO_FSAPT_DISP"):
        ref_wfn.set_variable("FSAPT_DISP_AB", matrices["Disp_AB"])

    if core.get_option("FISAPT", "SSAPT0_SCALE"):
        ref_wfn.set_variable("FSAPT_SINDAB_AB", matrices["sIndAB_AB"])
        ref_wfn.set_variable("FSAPT_SINDBA_AB", matrices["sIndBA_AB"])

        if core.get_option("FISAPT", "FISAPT_DO_FSAPT_DISP"):
            ref_wfn.set_variable("FSAPT_SDISP_AB", matrices["sDisp_AB"])
    return


def fisapt_plot(self):
    """Filesystem wrapper for FISAPT::plot."""

    filepath = core.get_option("FISAPT", "FISAPT_PLOT_FILEPATH")
    os.makedirs(filepath, exist_ok=True)

    geomfile = filepath + os.sep + "geom.xyz"
    xyz = self.molecule().to_string(dtype="xyz", units="Angstrom")
    with open(geomfile, "w") as fh:
        fh.write(xyz)

    self.raw_plot(filepath)


def _drop(array, filepath):
    """Helper to drop array to disk. FISAPT::drop

    Parameters
    ----------
    array : psi4.core.Matrix or psi4.core.Vector
        Matrix or vector to be written disk in plain text.
    filepath : str
        Full or partial file path. `array` will be written
        to <filepath>/<array.name>.dat.

    Returns
    -------
    None

    Notes
    -----
    Equivalent to https://github.com/psi4/psi4archive/blob/master/psi4/src/psi4/fisapt/fisapt.cc#L4389-L4420

    """
    filename = filepath + os.sep + array.name + ".dat"
    print("    Writing F-SAPT output file: {}".format(filename))
    with open(filename, "wb") as handle:
        np.savetxt(handle, array.to_array(), fmt="%24.16E", delimiter=" ", newline="\n")


core.FISAPT.compute_energy = fisapt_compute_energy
core.FISAPT.fdrop = fisapt_fdrop
core.FISAPT.plot = fisapt_plot
core.FISAPT.save_variables_to_wfn = fisapt_variables_to_wfn
