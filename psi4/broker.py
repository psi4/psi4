# Copyright (c) 2013-2018 Thomas Spura <thomas.spura@gmail.com>
from __future__ import print_function

import string
import sys
import time
import numpy as np
from numpy.testing import assert_array_almost_equal as assert_equal

import psi4
try:
    from ipi.interfaces.sockets import Client
    ipi_available = True
except ImportError:
    ipi_available = False
    # Define Client to enable testing of the Broker in the unittests
    class Client(object): pass


class Broker(Client):
    def __init__(self, options, genbas=None, serverdata=False):
        self.serverdata = serverdata
        if not ipi_available:
            print("i-pi is not available for import:",
                  "The broker infrastructure will not be available!")
            super(Broker, self).__init__()
        elif serverdata:
            mode, address, port = serverdata.split(":")
            mode = string.lower(mode)
            super(Broker, self).__init__(address=address, port=port,
                                      mode=mode)
        else:
            super(Broker, self).__init__(_socket=False)
        self.options = options

        print("PSI4 options:")
        for item, value in self.options.items():
            print(item, value)
            if item not in ["LOT", "multiplicity", "charge"]:
                psi4.core.set_global_option(item, value)
        psi4.core.IO.set_default_namespace("xwrapper")

        self.timing = {}

    def run_frc(self, pos):
        """Fetch force, energy of PSI.
        """
        if len(pos.shape) == 1:
            pos = pos.reshape((len(pos) / 3, 3))
        self.frc, self.pot = self.callback(pos)
        return self.frc, self.pot

    def get_molecule(self, pos):
        ret = "".join(
                    ["%s %f %f %f\n" % (self.atoms_list[i],
                                        pos[i][0], pos[i][1],
                                        pos[i][2])
                                        for i in range(len(pos))]
                    + ["units bohr\n"]
                    + ["no_reorient\n"]
                    #+ ["no_com\n"]
                    + ["symmetry c1\n"]
                    + ["%d %d\n" % (self.options["charge"],
                                    self.options["multiplicity"])]
                )
        return ret

    def callback(self, pos):
        """Initialize psi with new positions and calculate force.
        """

        mol = self.get_molecule(pos)
        psi4.geometry(mol, "xwrapper")

        self.calc_full(self.options["LOT"])

        self.pot = psi4.core.get_variable('CURRENT ENERGY')
        self.frc = -np.array(self.grd)

        return self.frc, np.float64(self.pot)

    def calc_full(self, LOT, bypass_scf=False):
        """Calculate the gradient with @LOT.

        When bypass_scf=True a hf energy calculation has been done before.
        """
        start = time.time()
        self.grd = psi4.gradient(LOT, bypass_scf=bypass_scf)
        time_needed = time.time() - start
        self.timing[LOT] = self.timing.get(LOT, []) + [time_needed]


def broker(serverdata=False, options=None):
    """ Run Broker to connect to i-pi

    Arguments:
        dryrun: Calculate forces with read in positions and the mirror image only.
    """
    b = Broker(serverdata=serverdata, options=options)

    mol = psi4.core.get_active_molecule()
    atoms = np.array(mol.geometry())
    names = [mol.symbol(i) for i in range(len(atoms))]

    print("Initial atoms", names)
    b.atoms_list = names

    try:
        if b.serverdata:
            b.run()
        else:
            print("ATOMS", atoms)
            if len(atoms.shape) == 1:
                ratio = len(atoms)/3
                atoms = atoms.reshape((ratio, 3))
                print("Calculating force for", atoms)
            print("FORCE:")
            frc, pot = b.run_frc(atoms)
            print(frc, pot)

            atoms *= -1.0
            frc2, pot2 = b.run_frc(atoms)
            print("FORCE MIRROR:")
            print(frc2, pot2)
            assert_equal(pot, pot2)
            assert_equal(frc, -1.0*frc2)

    except KeyboardInterrupt:
        print("Killing Broker")
        b.__del__()
        sys.exit(1)
    return b
