#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2017 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

from __future__ import absolute_import
from __future__ import print_function
#import os
#import re
#import math
#import copy
import itertools
from .molecule import Molecule
#from periodictable import *
#from physconst import *
from .vecutil import *
from .exceptions import *

#LINEAR_A_TOL = 1.0E-2  # When sin(a) is below this, we consider the angle to be linear
#DEFAULT_SYM_TOL = 1.0E-8
#FULL_PG_TOL = 1.0e-8
#ZERO = 1.0E-14
NOISY_ZERO = 1.0E-8
COORD_ZERO = 1.0E-5  # tolerance in coordinate alignment btwn qc programs


class OrientMols(object):
    """Class to encode a transformation between two molecular coordinate
    systems. After initializing with two qcdb.Molecule objects at the
    same geometry in possible different frames and orderings, class
    can apply the appropriate transformations to coordinate, gradient,
    Hessian, etc. arrays.

    """

    def __init__(self, molPermanent, molChangeable):
        """Stores the shift, rotation, axis exchange, axis inversion,
        and atom remapping necessary to bring the geometry of
        *molChangeable* into coincidence with the geometry of
        *molPermanent*. *molPermanent* and *molChangeable* must be
        :py:class:`qcdb.Molecule` and represent the same geometry.

        """
        # <<< Permanent (Psi4) >>>

        # Molecule
        self.Pmol = molPermanent
        # Vector to shift Pmol to center of mass
        self.Pshift = []
        # Matrix to rotate Pmol to inertial frame
        self.Protate = []

        # <<< Changeable (Cfour) >>>

        # Molecule
        self.Cmol = molChangeable
        # Vector to shift Cmol to center of mass
        self.Cshift = []
        # Matrix to rotate Cmol to inertial frame
        self.Crotate = []
        # Matrix to rotate Cmol to axis representation of Pmol
        self.Cexchflip = []
        # Vector to map Cmol to atom ordering of Pmol
        self.Catommap = []

        try:
            if ((self.Pmol.nallatom() == self.Cmol.nallatom()) and \
               (abs(self.Pmol.nuclear_repulsion_energy() - self.Cmol.nuclear_repulsion_energy()) < 1.0e-3)):
                self.create_orientation_from_molecules(self.Pmol, self.Cmol)
            else:
                print('qcdb.orient.__init__ debug info')
                self.Pmol.print_out()
                print('natom', self.Pmol.natom(), 'NRE', self.Pmol.nuclear_repulsion_energy(), 'rotor', self.Pmol.rotor_type())
                self.Cmol.print_out()
                print('natom', self.Cmol.natom(), 'NRE', self.Cmol.nuclear_repulsion_energy(), 'rotor', self.Cmol.rotor_type())
                raise ValidationError("""OrientMols Molecule arguments differ fatally.""")
        except AttributeError:
            raise ValidationError("""OrientMols must be instantiated with two qcdb.Molecule objects.""")

    def __str__(self):
        text = """  ==> qcdb OrientMols <==\n\n"""
        text += """   natom:     %d\n\n""" % (self.Pmol.natom())
        text += """   PNRE:       %16.8f\n""" % (self.Pmol.nuclear_repulsion_energy())
        text += """   Pshift:    %s\n""" % (self.Pshift)
        text += """   Protate:   %s\n""" % (self.Protate)
        text += """\n   CNRE:       %16.8f\n""" % (self.Cmol.nuclear_repulsion_energy())
        text += """   Cshift:    %s\n""" % (self.Cshift)
        text += """   Crotate:   %s\n""" % (self.Crotate)
        text += """   Cexchflip: %s\n""" % (self.Cexchflip)
        text += """   Catommap:  %s\n""" % (self.Catommap)
        return text

    def create_orientation_from_molecules(self, Pmol, Cmol):
        """Finds the shift, rotation, axis exchange, axis inversion,
        and atom remapping necessary to bring the geometry of *Cmol*
        into coincidence with the geometry of *Pmol*. *Pmol* and *Cmol*
        must be :py:class:`qcdb.Molecule` and represent the same
        geometry. Presently catches some errors of orientation that
        Cfour as *Cmol* should properly fulfill. These are unnecessary
        restrictions and can be relaxed later.

        """
        p4mol = Pmol.clone()
        c4mol = Cmol.clone()
        Nat = p4mol.natom()
        eye3 = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]

        # Find translation to CoM, straightforward
#        com = [ x if abs(x) > NOISY_ZERO else 0.0 for x in com]
        com = p4mol.center_of_charge()  # avoids trouble of ghost atoms w/center_of_mass
        p4mol.translate(scale(com, -1.0))
        self.Pshift = com

        com = c4mol.center_of_charge()
        c4mol.translate(scale(com, -1.0))
        self.Cshift = com

        # Find rotation to MoI frame, straightforward
#        frame = [[ x if abs(x) > NOISY_ZERO else 0.0 for x in ax] for ax in frame]
        moi, frame = p4mol.inertial_system(masswt=False, zero=NOISY_ZERO)  # avoids trouble of ghost atoms w/masswt'd
        Psort = sorted(range(3), key=lambda x: moi[x])
        p4mol.rotate(frame)
        self.Protate = frame

        moi, frame = c4mol.inertial_system(masswt=False, zero=NOISY_ZERO)
        Csort = sorted(range(3), key=lambda x: moi[x])
        c4mol.rotate(frame)
        self.Crotate = frame

        # Find degrees of freedom among axis exchanges
        rotor = p4mol.rotor_type()
#        if rotor != c4mol.rotor_type():
#            raise ValidationError("""molPermanent (%s) and molChangeable (%s) of different rotor types.""" % \
#                (rotor, c4mol.rotor_type()))
        # TODO: is this safe? differences in masses of ghost atoms can cause different rotor types
        if rotor == 'RT_ATOM':
            freebytop = []
        elif rotor == 'RT_LINEAR':                 # 0  <  IB == IC      inf > B == C
            freebytop = [1, 2]
        elif rotor == 'RT_SPHERICAL_TOP':          # IA == IB == IC       A == B == C
            freebytop = [0, 1, 2]
        elif rotor == 'RT_PROLATE_SYMMETRIC_TOP':  # IA <  IB == IC       A >  B == C
            freebytop = [1, 2]
        elif rotor == 'RT_OBLATE_SYMMETRIC_TOP':   # IA == IB <  IC       A == B >  C
            freebytop = [0, 1]
        elif rotor == 'RT_ASYMMETRIC_TOP':         # IA <  IB <  IC       A >  B >  C
            freebytop = []

        # Find possible mappings of axis exchange and flipping that brings Cgeom into coincidence with Pgeom
        axExch = [[], [], []]
        axPhse = [[], [], []]

        Pgeom = p4mol.geometry()
        Cgeom = c4mol.geometry()
        for Paxs in range(3):
            allowed = freebytop if Paxs in freebytop else [Paxs]

            for Caxs in allowed:
                PcolS = sorted([row[Psort[Paxs]] for row in Pgeom])
                CcolS = sorted([row[Csort[Caxs]] for row in Cgeom])
                CcolMS = sorted([-row[Csort[Caxs]] for row in Cgeom])

                if all([abs(CcolS[at] - PcolS[at]) < COORD_ZERO for at in range(Nat)]):
                    if all([abs(CcolS[at] - CcolMS[at]) < COORD_ZERO for at in range(Nat)]):
                        axPhse[Psort[Paxs]] = [1, -1]  # Indeterminate when symm P4 col == C4 col
                        axExch[Psort[Paxs]].append(Csort[Caxs])
                    else:
                        axPhse[Psort[Paxs]] = [1]  # Pos when asym P4 col == C4 col
                        axExch[Psort[Paxs]].append(Csort[Caxs])
                elif all([abs(CcolMS[at] - PcolS[at]) < COORD_ZERO for at in range(Nat)]):
                    axPhse[Psort[Paxs]] = [-1]  # Neg when P4 col == -C4 col
                    axExch[Psort[Paxs]].append(Csort[Caxs])

            if len(axPhse[Psort[Paxs]]) == 0:
                print('qcdb.orient.create_orientation_from_molecules debug info')
                print('\nrotor', rotor, 'Paxs', Paxs, 'Caxs', Caxs, \
                    'allowed', allowed, 'P(axs)', Psort.index(Paxs), 'C(Caxs)', Csort.index(Caxs), \
                    '\nPcolS: ', PcolS, '\nCcolS: ', CcolS, '\nCcolMS: ', CcolMS, \
                    'axExch', axExch, 'axPhse', axPhse)
                print('\nPgeom: ')
                for item in Pgeom:
                    print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))
                print('\nCgeom: ')
                for item in Cgeom:
                    print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))
                print(self)
                raise ValidationError("""Axis unreconcilable between QC programs.""")

        allowedExchflip = []
        for lee in itertools.product(*axExch):
            for lpp in itertools.product(*axPhse):
                if (sorted(lee) == [0, 1, 2]):
                    temp = zero(3, 3)
                    for ax in range(3):
                        temp[lee[ax]][ax] = lpp[ax]
                    allowedExchflip.append(temp)

        # Find mapping of atom exchange that brings Cgeom into coincidence with Pgeom
        for exfp in allowedExchflip:
            Cgeom = mult(c4mol.geometry(), exfp)

            mapMat = [None] * Nat
            Pwhite = list(range(Nat))
            Cwhite = list(range(Nat))

#            print '\nPgeom: '
#            for item in Pgeom:
#                print '       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2])
#            print '\nCgeom: '
#            for item in Cgeom:
#                print '       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2])

            while len(Pwhite) > 0:
                Patm = Pwhite[0]
                sameElem = [at for at in range(Nat) if ((c4mol.symbol(at) == p4mol.symbol(Patm)) or (c4mol.symbol(at) == 'GH'))]
                allowed = list(set(sameElem) & set(Cwhite))

                for Catm in allowed:
                    if all([abs(Cgeom[Catm][ax] - Pgeom[Patm][ax]) < COORD_ZERO for ax in range(3)]):
                        mapMat[Patm] = Catm
                        Pwhite.remove(Patm)
                        Cwhite.remove(Catm)
#                        print 'matchd on atom', 'Patm', Patm, 'Catm', Catm
                        break
                else:
#                    print 'failed on atom', 'Patm', Patm, 'rejecting', exfp
                    break
            else:
#                print 'accept exchflp', exfp, 'with map', mapMat
                break
        else:
            print('else of for', exfp, mapMat)
            print('qcdb.orient.create_orientation_from_molecules debug info')
            print('\nPatm', Patm, 'Catm', Catm, 'Pwhite', Pwhite, \
                'Cwhite', Cwhite, 'sameElem', sameElem, 'allowed', allowed, \
                '\nCgeom[Catm]: ', Cgeom[Catm], '\nPgeom[Patm]: ', Pgeom[Patm], \
                '\nmapMat', mapMat)
            print('\nPgeom: ')
            for item in Pgeom:
                print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))
            print('\nCgeom: ')
            for item in Cgeom:
                print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))
            print(self)
            raise ValidationError("""Atom unreconcilable between QC programs.""")

        self.Cexchflip = exfp
        self.Catommap = mapMat
        # Note that this is resetting the geom but not the atoms, so c4mol.print_out() is deceptive
        c4mol.rotate(exfp)
        new_geom = []
        for at in range(Nat):
            new_geom.append(Cgeom[mapMat[at]])
        c4mol.set_geometry(new_geom)

        # One last check that p4mol and c4mol align
        Pgeom = p4mol.geometry()
        Cgeom = c4mol.geometry()

        if not all([all([abs(Cgeom[at][ax] - Pgeom[at][ax]) < COORD_ZERO for ax in range(3)]) for at in range(Nat)]):
            raise ValidationError("""Geometries unreconcilable between QC programs:\n  P4 %s\n  C4 %s""" % (Pgeom, Cgeom))

    def transform_coordinates(self, coord):
        """

        """
        #print self
        print("Original")
        coord.print_out()

        coord.translate(scale(self.Cshift, -1.0))
        print("Shift")
        coord.print_out()

        coord.rotate(self.Crotate)
        print("Rotate")
        coord.print_out()

        coord.rotate(self.Cexchflip)
        print("ExchFlip")
        coord.print_out()

        geom = coord.geometry()
        new_geom = []
        for at in range(coord.natom()):
            new_geom.append(geom[self.Catommap[at]])
        coord.set_geometry(new_geom)
        print("AtomMap")
        coord.print_out()

        coord.rotate(transpose(self.Protate))
        print("P4 Rotate")
        coord.print_out()

        coord.translate(self.Pshift)
        print("P4 Shift")
        coord.print_out()

    def transform_coordinates2(self, coord):
        """

        """
        geom = coord.geometry()

#        print self
#        print "Original"
#        for item in geom:
#            print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))

        # ?????
        #coord.translate(scale(self.Cshift, -1.0))
        geom2 = []
        for item in geom:
            geom2.append(sub(item, self.Cshift))
        geom = geom2
#        print "Shift"
#        for item in geom:
#            print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))

        #coord.rotate(self.Crotate)
        geom2 = mult(geom, self.Crotate)
        geom = geom2
#        print "Rotate"
#        for item in geom:
#            print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))

        #coord.rotate(self.Cexchflip)
        geom2 = mult(geom, self.Cexchflip)
        geom = geom2
 #       print "ExchFlip"
 #       for item in geom:
 #           print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))

        geom2 = []
        for at in range(coord.natom()):
            geom2.append(geom[self.Catommap[at]])
        coord.set_geometry(geom2)
        geom = geom2
#        print "AtomMap"
#        for item in geom:
#            print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))

        #coord.rotate(transpose(self.Protate))
        geom2 = mult(geom, transpose(self.Protate))
        geom = geom2
#        print "P4 Rotate"
#        for item in geom:
#            print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))

        geom2 = []
        for item in geom:
            geom2.append(add(item, self.Pshift))
        geom = geom2
#        print "P4 Shift"
#        for item in geom:
#            print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))

        Pgeom = self.Pmol.geometry()
        Cgeom = geom
        Nat = len(geom)
        if not all([all([abs(Cgeom[at][ax] - Pgeom[at][ax]) < COORD_ZERO for ax in range(3)]) for at in range(Nat)]):
            raise ValidationError("""Geometries unreconcilable between QC programs:\n  P4 %s\n  C4 %s""" % (Pgeom, Cgeom))
        return geom

    def transform_gradient(self, arr):
        """Applies to *arr* the transformation appropriate to bring a
        gradient in *molChangeable* orientation into *molPermanent*
        orientation. In particular, applies a rotation to place it
        in the inertial frame, a column exchange and phasing to place
        it in the axis system, a row exchange to place it in the atom
        ordering, and a rotation to remove it from the inertial frame.

        """
        arr = mult(arr, self.Crotate)
#        print "Rotate"
#        for item in arr:
#            print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))

        arr = mult(arr, self.Cexchflip)
#        print "ExchFlip"
#        for item in arr:
#            print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))

        arr2 = []
        for at in range(len(arr)):
            arr2.append(arr[self.Catommap[at]])
        arr = arr2
#        print "AtomMap"
#        for item in arr:
#            print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))

        arr = mult(arr, transpose(self.Protate))
#        print "P4 Rotate"
#        for item in arr:
#            print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))

        return arr

    def transform_elementlist(self, elem):
        """Applies to *elem* the transformation appropriate to bring a
        list of atomic numbers in *molChangeable* orientation into
        *molPermanent* orientation. In particular, applies a row exchange
        to place it in the atom ordering.

        """
        elem2 = []
        for at in range(self.Pmol.natom()):
            elem2.append(elem[self.Catommap[at]])
        return elem2

    def transform_vector(self, vec):
        """Applies to *vec* the transformation appropriate to bring a
        vector in *molChangeable* orientation into *molPermanent*
        orientation. In particular, applies a rotation to place it
        in the inertial frame, a column exchange and phasing to place
        it in the axis system, and a rotation to remove it from the
        inertial frame.

        """
        vec = [vec]  # hack since vecutil handles matrices, not vectors

        #coord.rotate(self.Crotate)
        vec2 = mult(vec, self.Crotate)
        vec = vec2
        #print "Rotate"
        #for item in vec:
        #    print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))

        #coord.rotate(self.Cexchflip)
        vec2 = mult(vec, self.Cexchflip)
        vec = vec2
        #print "ExchFlip"
        #for item in vec:
        #    print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))

        #coord.rotate(transpose(self.Protate))
        vec2 = mult(vec, transpose(self.Protate))
        vec = vec2
        #print "P4 Rotate"
        #for item in vec:
        #    print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))

        return vec[0]
