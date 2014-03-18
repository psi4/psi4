#
#@BEGIN LICENSE
#
# PSI4: an ab initio quantum chemistry software package
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
#@END LICENSE
#

import math

import qcdb.exceptions
import qcformat
#import molpro_basissets


class Psi4In(qcformat.InputFormat):

    def __init__(self, mem, mtd, bas, mol, sys, cast):
        qcformat.InputFormat.__init__(self, mem, mtd, bas, mol, sys, cast)

        # memory in MB --> MW
        #self.memory = int(math.ceil(mem / 8.0))
        # auxiliary basis sets
        #[self.unaugbasis, self.augbasis, self.auxbasis] = self.corresponding_aux_basis()

    def format_global_parameters(self):
        text = ''

        #if self.method in ['mp2c', 'dft-sapt-shift', 'dft-sapt', 'dft-sapt-pbe0ac', 'dft-sapt-pbe0acalda']:
        #    text += """GTHRESH,ZERO=1.e-14,ONEINT=1.e-14,TWOINT=1.e-14,ENERGY=1.e-8,ORBITAL=1.e-8,GRID=1.e-8\n\n"""
        #elif self.method in ['b3lyp', 'b3lyp-d', 'df-b3lyp', 'df-b3lyp-d']:
        #    text += """GTHRESH,ZERO=1.e-14,ONEINT=1.e-14,TWOINT=1.e-14,ENERGY=1.e-8,ORBITAL=1.e-7,GRID=1.e-8\n\n"""
        #else:
        #    text += """GTHRESH,ZERO=1.e-14,ONEINT=1.e-14,TWOINT=1.e-14,ENERGY=1.e-9\n\n"""

        return text

    def format_basis(self):
        text = ''
        text += """set basis %s\n\n""" % (self.basis)
#        text += """basis={\n"""
#
#        try:
#            # jaxz, maxz, etc.
#            for line in molpro_basissets.altbasis[self.basis]:
#                text += """%s\n""" % (line)
#            text += '\n'
#        except KeyError:
#            # haxz
#            if self.basis.startswith('heavy-aug-'):
#                text += """set,orbital; default,%s,H=%s\n""" % (self.basis[6:], self.unaugbasis)
#            # xz, axz, 6-31g*
#            else:
#                text += """set,orbital; default,%s\n""" % (self.basis)
#
#        if ('df-' in self.method) or ('f12' in self.method) or (self.method in ['mp2c', 'dft-sapt', 'dft-sapt-pbe0acalda']):
#            if self.unaugbasis and self.auxbasis:
#
#                text += """set,jkfit;   default,%s/jkfit\n""" % (self.auxbasis)
#                text += """set,jkfitb;  default,%s/jkfit\n""" % (self.unaugbasis)
#                text += """set,mp2fit;  default,%s/mp2fit\n""" % (self.auxbasis)
#                text += """set,dflhf;   default,%s/jkfit\n""" % (self.auxbasis)
#            else:
#                raise qcdb.exceptions.ValidationError("""Auxiliary basis not predictable from orbital basis '%s'""" % (self.basis))
#
#        text += """}\n\n"""
        return text

    def format_infile_string(self):
        text = ''

        # format comment and memory
        text += """# %s %s\n""" % (self.index, self.molecule.tagline)
        text += """memory %d mb\n""" % (self.memory)

        # format molecule, incl. charges and dummy atoms
        text += self.molecule.format_molecule_for_psi4()

        # format global convergence directions
        #text += self.format_global_parameters()

        # format castup directions
        if self.castup is True:
            text += """set basis_guess sto-3g\n"""

        # format basis set
        text += self.format_basis()

        # format method
        for line in qcmtdIN[self.method]:
            text += """%s\n""" % (line)
        text += """print_variables()\n"""
        text += '\n'

        return text


qcmtdIN = {
'sapt0': [
    """energy('sapt0')"""],
'saptccd': [
    """set scf_type df""",
    """set guess sad""",
    """set freeze_core true""",
    """set nat_orbs_t2 true""",
    """set nat_orbs_t3 true""",
    """set nat_orbs_v4 true""",
    """energy('sapt2+3(ccd)')"""],

}
