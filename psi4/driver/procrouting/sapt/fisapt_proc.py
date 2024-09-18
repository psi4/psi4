#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2024 The Psi4 Developers.
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


def fisapt_compute_energy(self, external_potentials=None):
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
    self.coulomb()
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
        self.disp(self.matrices(), self.vectors(), True)  # Expensive, only do if needed  # unteseted translation of below
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
        #else:
        #    # Build Empirical Dispersion
        #    dashD = empirical_dispersion.EmpiricalDispersion(name_hint='SAPT0-D3M')
        #    dashD.print_out()
        #    # Compute -D
        #    Edisp = dashD.compute_energy(core.get_active_molecule())
        #    core.set_variable('{} DISPERSION CORRECTION ENERGY'.format(dashD.fctldash), Edisp)            # Printing
        #    text = []
        #    text.append("   => {}: Empirical Dispersion <=".format(dashD.fctldash.upper()))
        #    text.append(" ")
        #    text.append(dashD.description)
        #    text.append(dashD.dashlevel_citation.rstrip())
        #    text.append("\n    Empirical Dispersion Energy [Eh] =     {:24.16f}\n".format(Edisp))
        #    text.append('\n')
        #    core.print_out('\n'.join(text))
        if core.get_option("FISAPT", "FISAPT_FSAPT_PSI_VARIABLES"):
            self.save_fsapt_variables(external_potentials)
        else:
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

    filepath = core.get_option("FISAPT", "FISAPT_FSAPT_FILEPATH")
    os.makedirs(filepath, exist_ok=True)

    core.print_out("    F-SAPT Data Filepath = {}\n\n".format(filepath))

    geomfile = filepath + os.sep + 'geom.xyz'
    xyz = self.molecule().to_string(dtype='xyz', units='Angstrom')
    with open(geomfile, 'w') as fh:
        fh.write(xyz)

    # write external potential geometries
    if external_potentials is not None and isinstance(external_potentials, dict):
        for frag in "ABC":
            potential = external_potentials.get(frag, None)
            if potential is not None:
                xyz = str(len(potential)) + "\n\n"
                for qxyz in potential:
                    if len(qxyz) == 2:
                        xyz += "Ch %f %f %f\n" % (qxyz[1][0], qxyz[1][1], qxyz[1][2])
                    elif len(qxyz) == 4:
                        xyz += "Ch %f %f %f\n" % (qxyz[1], qxyz[2], qxyz[3])
                    else:
                        raise ValidationError(f"Point charge '{qxyz}' not mapping into 'chg, [x, y, z]' or 'chg, x, y, z'")

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

    _drop(vectors["ZA"], filepath)
    _drop(vectors["ZB"], filepath)
    _drop(matrices["Qocc0A"], filepath)
    _drop(matrices["Qocc0B"], filepath)
    _drop(matrices["Elst_AB"], filepath)
    _drop(matrices["Exch_AB"], filepath)
    _drop(matrices["IndAB_AB"], filepath)
    _drop(matrices["IndBA_AB"], filepath)

    if core.get_option("FISAPT", "FISAPT_DO_FSAPT_DISP"):
        matrices["Disp_AB"].name = "Disp"
        _drop(matrices["Disp_AB"], filepath)

    if core.get_option("FISAPT", "SSAPT0_SCALE"):
        ssapt_filepath = core.get_option("FISAPT", "FISAPT_FSSAPT_FILEPATH")
        os.makedirs(ssapt_filepath, exist_ok=True)

        core.print_out("    sF-SAPT Data Filepath = {}\n\n".format(ssapt_filepath))

        geomfile = ssapt_filepath + os.sep + 'geom.xyz'
        with open(geomfile, 'w') as fh:
            fh.write(xyz)

        matrices["sIndAB_AB"].name = "IndAB"
        matrices["sIndBA_AB"].name = "IndBA"

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
            _drop(matrices["sDisp_AB"], ssapt_filepath)

def fisapt_save_fsapt_variables(self, external_potentials=None):
    core.print_out("  ==> F-SAPT Output (to psi vars) <==\n\n")

    # write external potential geometries
    external_pot_str = ""
    if external_potentials is not None and isinstance(external_potentials, dict):
        for frag in "ABC":
            potential = external_potentials.get(frag, None)
            if potential is not None:
                xyz = str(len(potential)) + "\n\n"
                for qxyz in potential:
                    if len(qxyz) == 2:
                        xyz += "Ch %f %f %f\n" % (qxyz[1][0], qxyz[1][1], qxyz[1][2])
                    elif len(qxyz) == 4:
                        xyz += "Ch %f %f %f\n" % (qxyz[1], qxyz[2], qxyz[3])
                    else:
                        raise ValidationError(f"Point charge '{qxyz}' not mapping into 'chg, [x, y, z]' or 'chg, x, y, z'")
                external_pot_str += xyz

    vectors = self.vectors()
    matrices = self.matrices()

    core.set_variable("QA", matrices["Qocc0A"].to_array())
    core.set_variable("QB", matrices["Qocc0B"].to_array())
    core.set_variable("Elst_AB", matrices["Elst_AB"].to_array())
    core.set_variable("Exch_AB", matrices["Exch_AB"].to_array())
    core.set_variable("IndAB_AB", matrices["IndAB_AB"].to_array())
    core.set_variable("IndBA_AB", matrices["IndBA_AB"].to_array())

    if core.get_option("FISAPT", "FISAPT_DO_FSAPT_DISP"):
        matrices["Disp_AB"].name = "Disp"
        core.set_variable("Disp_AB", matrices["Disp_AB"].to_array())

    if core.get_option("FISAPT", "SSAPT0_SCALE"):
        matrices["sIndAB_AB"].name = "IndAB"
        matrices["sIndBA_AB"].name = "IndBA"
        core.set_variable("sIndAB_AB", matrices["sIndAB_AB"].to_array())
        core.set_variable("sIndBA_AB", matrices["sIndBA_AB"].to_array())

        if core.get_option("FISAPT", "FISAPT_DO_FSAPT_DISP"):
            matrices["sDisp_AB"].name = "Disp"
            core.set_variable("sDisp_AB", matrices["sDisp_AB"].to_array())

def fisapt_plot(self):
    """Filesystem wrapper for FISAPT::plot."""

    filepath = core.get_option("FISAPT", "FISAPT_PLOT_FILEPATH")
    os.makedirs(filepath, exist_ok=True)

    geomfile = filepath + os.sep + 'geom.xyz'
    xyz = self.molecule().to_string(dtype='xyz', units='Angstrom')
    with open(geomfile, 'w') as fh:
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
    filename = filepath + os.sep + array.name + '.dat'
    with open(filename, 'wb') as handle:
        np.savetxt(handle, array.to_array(), fmt="%24.16E", delimiter=' ', newline='\n')


core.FISAPT.compute_energy = fisapt_compute_energy
core.FISAPT.fdrop = fisapt_fdrop
core.FISAPT.plot = fisapt_plot
core.FISAPT.save_fsapt_variables = fisapt_save_fsapt_variables
