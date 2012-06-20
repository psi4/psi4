r"""Module to define a class :py:class:`~BasisFamily` that associates 
fitting basis sets to an orbital basis and to provide functions to
query appropriate fitting bases for any orbital basis distributed
with Psi4.

"""
import os
import PsiMod
from psiexceptions import *


basisfamily_list = []


class BasisFamily(object):
    """Class to associate with an orbital basis name *ornate*
    the gbs file names in which the orbital basis *orbital* 
    (usually the coded form of *ornate*) and *jkfit*, *rifit*, 
    and *dualfit* auxiliary bases can be found.

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
    temp = PsiMod.BasisSet.make_filename(name)
    return os.path.splitext(os.path.splitext(temp)[0])[0]


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

    for fam in basisfamily_list:
        PsiMod.print_out('%s' % fam)


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
