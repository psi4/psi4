#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2016 The Psi4 Developers.
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
    def __init__(self, ornate, orbital=None, jk=None, ri=None, dual=None):
        self.ornate = ornate
        if orbital is None:
            self.orbital = sanitize_basisname(ornate)
        else:
            self.orbital = sanitize_basisname(orbital)
        self.jkfit = jk
        self.rifit = ri
        self.dualfit = dual

    def __str__(self):
        text = ''
        text += """  ==> %s Family <==\n\n""" % (self.ornate)
        text += """  Orbital basis:        %s\n""" % (self.orbital)
        text += """  JK auxiliary basis:   %s\n""" % (self.jkfit)
        text += """  MP2 auxiliary basis:  %s\n""" % (self.rifit)
        text += """  DUAL auxiliary basis: %s\n""" % (self.dualfit)
        text += """\n"""
        return text

    def name(self):
        """Function to return the ornate name of the orbital basis,
        e.g., 6-311++G** for 6-311ppgss.
        """
        return self.ornate

    def add_jfit(self, fit):
        """Function to add basis *fit* as associated fitting basis
        member *jfit* to a BasisFamily object.
        """
        self.jfit = sanitize_basisname(fit)

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
    from p4util.basislistdunning import load_basfam_dunning
    from p4util.basislistother import load_basfam_other

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


def corresponding_orbital(name):
    """Function to validate if the orbital basis *name* in coded or
    ornate form is in Psi4's standard installed bases list. ``None``
    is returned if the orbital basis is not found.
    """
    basisfamily_list = load_basis_families()

    for fam in basisfamily_list:
        if sanitize_basisname(fam.ornate) == sanitize_basisname(name):
            return fam.ornate
    return None

def corresponding_jfit(name):
    """Function to return an appropriate J fitting basis for
        the orbital basis *name* in coded or ornate form. ``None``
        is returned if no fitting basis is defined or if the
        orbital basis is not found.
        """
    basisfamily_list = load_basis_families()
    
    for fam in basisfamily_list:
        if sanitize_basisname(fam.ornate) == sanitize_basisname(name):
            return fam.jfit
    return None


def corresponding_jkfit(name):
    """Function to return an appropriate JK fitting basis for
    the orbital basis *name* in coded or ornate form. ``None``
    is returned if no fitting basis is defined or if the
    orbital basis is not found.
    """
    basisfamily_list = load_basis_families()

    for fam in basisfamily_list:
        if sanitize_basisname(fam.ornate) == sanitize_basisname(name):
            return fam.jkfit
    return None


def corresponding_rifit(name):
    """Function to return an appropriate RI fitting basis for
    the orbital basis *name* in coded or ornate form. ``None``
    is returned if no fitting basis is defined or if the
    orbital basis is not found.
    """
    basisfamily_list = load_basis_families()

    for fam in basisfamily_list:
        if sanitize_basisname(fam.ornate) == sanitize_basisname(name):
            return fam.rifit
    return None


def corresponding_dualfit(name):
    """Function to return an appropriate DUAL helper basis for
    the orbital basis *name* in coded or ornate form. ``None``
    is returned if no fitting basis is defined or if the
    orbital basis is not found.
    """
    basisfamily_list = load_basis_families()

    for fam in basisfamily_list:
        if sanitize_basisname(fam.ornate) == sanitize_basisname(name):
            return fam.dualfit
    return None