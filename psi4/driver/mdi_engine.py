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

#from psi4 import core
import psi4
from psi4.driver import p4util
from psi4.driver import driver
from psi4.driver.p4util.exceptions import *

from MDI_Library.mdi import MDI_Init, MDI_Get_Intra_Code_MPI_Comm, MDI_Accept_Communicator
from MDI_Library.mdi import MDI_Send, MDI_Recv, MDI_Recv_Command, MDI_Conversion_Factor
from MDI_Library.mdi import MDI_INT, MDI_DOUBLE, MDI_CHAR

class MDIEngine():

    def __init__(self):

        self.molecule = psi4.core.get_active_molecule()
        self.energy = 0.0

        # lattice variables
        self.nlattice = 0 # number of lattice point charges
        self.clattice = [] # list of lattice coordinates
        self.lattice = [] # list of lattice charges
        self.lattice_field = psi4.QMMM() # Psi4 chargefield




def mdi():
    """ Begin functioning as an MDI engine

#    Arguments:
#        molecule: Initial molecule
#        serverdata: Configuration where to connect to ipi
#        options: LOT, multiplicity and charge
    """
    core.print_out("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n")
    pass
