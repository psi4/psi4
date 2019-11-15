#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2019 The Psi4 Developers.
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

def fisapt_compute_energy(self):
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
        self.disp(matrices_, vectors_, true)  # Expensive, only do if needed
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
        if core.get_option("FISAPT", "FISAPT_DO_EMPIRICAL_DISP"):
            # Build Empirical Dispersion
            dashlvl = core.get_option("FISAPT", "FISAPT_EMPIRICAL_DISP")
            dashD = empirical_dispersion.EmpiricalDispersion(name_hint=f'SAPT0-{dashlvl}')
            dashD.print_out()
            # Compute -D
            NA = self.molecule().extract_subsets(1).natom()
            NB = self.molecule().extract_subsets(2).natom()
            dimer_D3 = dashD.compute_energy(self.molecule())
            d3pairs = core.variables()['PAIRWISE DISPERSION CORRECTION ANALYSIS'].to_array()

            # Rebuild dispersion energy
            Edisp = 0.0
            for i in range(NA):
                for j in range(NB):
                    Edisp += d3pairs[i, NA + j]

            core.set_variable('SAPT EMPIRICAL {} DISP ENERGY'.format(dashD.dashlevel.upper()), Edisp)
            # Printing
            core.print_out("    Empirical Dispersion Interaction Energy [Eh] =     {:24.16f}\n\n".format(Edisp))
        self.fdrop()

    # => Scalar-Field Analysis <=

    if core.get_option("FISAPT", "FISAPT_DO_PLOT"):
        core.timer_on("FISAPT:FSAPT:cubeplot")
        self.plot()
        core.timer_off("FISAPT:FSAPT:cubeplot")

    # => Summary <=

    self.print_trailer()


def fisapt_fdrop(self):
    """Drop output files from FSAPT calculation. FISAPT::fdrop"""

    core.print_out("  ==> F-SAPT Output <==\n\n")

    filepath = core.get_option("FISAPT", "FISAPT_FSAPT_FILEPATH")
    os.makedirs(filepath, exist_ok=True)

    core.print_out("    F-SAPT Data Filepath = {}\n\n".format(filepath))

    geomfile = filepath + os.sep + 'geom.xyz'
    xyz = self.molecule().to_string(dtype='xyz', units='Angstrom')
    with open(geomfile, 'w') as fh:
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
    else:
        # Populate (NA, NB) chunk of Disp_AB
        d3pairs = core.variables()['PAIRWISE DISPERSION CORRECTION ANALYSIS'].to_array()
        Ntot = self.molecule().natom()
        NA = self.molecule().extract_subsets(1).natom()
        NB = self.molecule().extract_subsets(2).natom()

        Disp_AB = np.zeros((NA, NB))
        for a in range(NA):
           for b in range(NB):
                B = b + NA - 1
                Disp_AB[a,b] = d3pairs[a,B]

        matrices["Disp_AB"] = core.Matrix.from_array(Disp_AB)
        matrices["Disp_AB"].name = 'D3Disp'
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
            _drop(matrices["Disp_AB"], ssapt_filepath)
        else:
            # Populate (NA, NB) chunk of Disp_AB
            d3pairs = core.variables()['PAIRWISE DISPERSION CORRECTION ANALYSIS'].to_array()
            Ntot = self.molecule().natom()
            NA = self.molecule().extract_subsets(1).natom()
            NB = self.molecule().extract_subsets(2).natom()

            Disp_AB = np.zeros((NA, NB))
            for a in range(NA):
               A = a + 1
               for b in range(NB):
                   B = b + NA + 1
                   Disp_AB[a,b] = d3pairs[A,B]

            matrices["sDisp_AB"] = core.Matrix.from_array(Disp_AB)
            matrices["sDisp_AB"].name = 'D3Disp'
            _drop(matrices["sDisp_AB"], ssapt_filepath)

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
