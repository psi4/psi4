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

r"""Module to define a class :py:class:`~BasisFamily` that associates
fitting basis sets to an orbital basis and to provide functions to
query appropriate fitting bases for any orbital basis distributed
with Psi4.

"""
import os


basisfamily_list = []


class BasisFamily(object):
    """Class to associate with an orbital basis name *ornate*
    the gbs file names in which the orbital basis *orbital*
    (usually the coded form of *ornate*) and *jfit*, *jkfit*,
    *rifit*, and *dualfit* auxiliary bases can be found.

    """

    def __init__(self, ornate, orbital=None):
        """Constructor"""
        # literature name of orbital basis set, e.g., aug-cc-pVTZ or 6-31+G*
        self.ornate = ornate
        # gbs file name of orbital basis set, e.g., aug-cc-pvtz or 6-31pgs
        self.orbital = sanitize_basisname(ornate) if orbital is None else sanitize_basisname(orbital)
        # gbs file name of JKFIT designed for orbital basis
        self.jkfit = None
        # gbs friendly file name of JFIT designed for orbital basis
        self.jfit = None
        # gbs file name of CFIT designed for orbital basis
        self.rifit = None
        # gbs file name of DUAL designed for orbital basis
        self.dualfit = None
        # gbs file name of JKFIT default when self.jkfit unavailable
        #self.jkdef = jkdef
        # gbs file name of JFIT default when self.jfit unavailable
        #self.jdef = jdef
        # gbs file name of CFIT default when self.rifit unavailable
        #self.ridef = ridef

    def __str__(self):
        text = ''
        text += """  ==> %s Family <==\n\n""" % (self.ornate)
        text += """  Orbital basis:        %s\n""" % (self.orbital)
        #text += """  JK auxiliary basis:   %s\n""" % (self.jkfit)
        #text += """  MP2 auxiliary basis:  %s\n""" % (self.rifit)
        text += """  JK auxiliary basis:   %s  Def: %s\n""" % (self.jkfit, self.jkdef)
        text += """  J auxiliary basis:    %s  Def: %s\n""" % (self.jfit, self.jdef)
        text += """  MP2 auxiliary basis:  %s  Def: %s\n""" % (self.rifit, self.ridef)
        text += """  DUAL auxiliary basis: %s\n""" % (self.dualfit)
        text += """\n"""
        return text

    def name(self):
        """Function to return the ornate name of the orbital basis,
        e.g., 6-311++G** for 6-311ppgss.

        """
        return self.ornate

    def add_jkfit(self, fit):
        """Function to add basis *fit* as associated fitting basis
        member *jkfit* to a BasisFamily object.

        """
        self.jkfit = sanitize_basisname(fit)

    def add_rifit(self, fit):
        """Function to add basis *fit* as associated fitting basis
        member *rifit* to a BasisFamily object.

        """
        self.rifit = sanitize_basisname(fit)

    def add_dualfit(self, fit):
        """Function to add basis *fit* as associated helper basis
        member *dualfit* to a BasisFamily object.

        """
        self.dualfit = sanitize_basisname(fit)

    def add_jfit(self, fit):
        """Function to add basis *fit* as associated fitting basis
        member *jfit* to a BasisFamily object.

        """
        self.jfit = sanitize_basisname(fit)

    def add_jfit_default(self, fit):
        """Function to add basis *fit* as associated fitting basis
        member *jdef* to a BasisFamily object.

        """
        self.jdef = sanitize_basisname(fit)

    def add_jkfit_default(self, fit):
        """Function to add basis *fit* as associated fitting basis
        member *jkdef* to a BasisFamily object.

        """
        self.jkdef = sanitize_basisname(fit)

    def add_rifit_default(self, fit):
        """Function to add basis *fit* as associated fitting basis
        member *ridef* to a BasisFamily object.

        """
        self.ridef = sanitize_basisname(fit)


def sanitize_basisname(name):
    """Function to return *name* in coded form, stripped of
    characters that confuse filenames, characters into lowercase,
    ``+`` into ``p``, ``*`` into ``s``, and ``(``, ``)``, & ``,``
    into ``_``.

    """
    temp = name.lower()
    temp = temp.replace('+', 'p')
    temp = temp.replace('*', 's')
    temp = temp.replace('(', '_')
    temp = temp.replace(')', '_')
    temp = temp.replace(',', '_')
    return temp


def load_basis_families():
    """Function to load into the array ``basisfamily_list``
    BasisFamily objects for all Psi4's standard installed bases.

    """
    from basislistdunning import load_basfam_dunning
    from basislistother import load_basfam_other

    if len(basisfamily_list) == 0:
        load_basfam_dunning()
        load_basfam_other()
    return basisfamily_list


def print_basis_families():
    """Function to print to the output file a formatted summary
    of all the BasisFamily objects in ``basisfamily_list``, by
    default all Psi4's standard installed bases.

    """
    basisfamily_list = load_basis_families()

    text = ''
    for fam in basisfamily_list:
        text += '%s' % (fam)
    return text


def corresponding_basis(name, role='BASIS'):
    """Function to validate if the orbital basis *name* in coded or
    ornate form is in Psi4's standard installed bases list. ``None``
    is returned if the orbital basis is not found.

    """
    role = role.upper()
    basisfamily_list = load_basis_families()

    for fam in basisfamily_list:
        if sanitize_basisname(fam.ornate) == sanitize_basisname(name):
            if role == 'ORNATE':
                return fam.ornate
            elif role == 'BASIS' or role == 'ORBITAL':
                return fam.orbital
            elif role == 'JFIT':
                return fam.jfit
            elif role == 'JKFIT':
                return fam.jkfit
            elif role == 'RIFIT':
                return fam.rifit
            elif role == 'DUALFIT':
                return fam.dualfit
    return None
