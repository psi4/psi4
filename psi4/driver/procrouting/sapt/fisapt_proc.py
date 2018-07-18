#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2018 The Psi4 Developers.
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

#from psi4.driver import p4util
#from psi4.driver import constants
#from psi4.driver.p4util.exceptions import ConvergenceError, ValidationError
from psi4 import core


def fisapt_compute_energy(self):
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
        core.timer_on("FISAPT:FSAPT:disp")
        self.fdisp()
        core.timer_off("FISAPT:FSAPT:disp")
        self.fdrop()

    # => Scalar-Field Analysis <=

    if core.get_option("FISAPT", "FISAPT_DO_PLOT"):
        core.timer_on("FISAPT:FSAPT:cubeplot")
        self.plot()
        core.timer_off("FISAPT:FSAPT:cubeplot")

    # => Summary <=

    self.print_trailer()


core.FISAPT.compute_energy = fisapt_compute_energy
