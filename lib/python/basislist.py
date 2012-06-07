r"""Module to define a class that associates
fitting basis sets to each orbital basis.
"""
import os
import PsiMod
from psiexceptions import *


basisfamily_list = []


class BasisFamily(object):

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
        return self.ornate

    def add_jkfit(self, fit):
        self.jkfit = sanitize_basisname(fit)

    def add_rifit(self, fit):
        self.rifit = sanitize_basisname(fit)

    def add_dualfit(self, fit):
        self.dualfit = sanitize_basisname(fit)


def sanitize_basisname(name):
    temp = PsiMod.BasisSet.make_filename(name)
    return os.path.splitext(os.path.splitext(temp)[0])[0]


def load_basis_families():
    from basislistdunning import load_basfam_dunning
    from basislistother import load_basfam_other

    if len(basisfamily_list) == 0:
        load_basfam_dunning()
        load_basfam_other()
    return basisfamily_list


def print_basis_families():
    basisfamily_list = load_basis_families()

    for fam in basisfamily_list:
        PsiMod.print_out('%s' % fam)


def corresponding_jkfit(name):
    basisfamily_list = load_basis_families()

    for fam in basisfamily_list:
        if sanitize_basisname(fam.ornate) == sanitize_basisname(name):
            return fam.jkfit
    return None


def corresponding_rifit(name):
    basisfamily_list = load_basis_families()

    for fam in basisfamily_list:
        if sanitize_basisname(fam.ornate) == sanitize_basisname(name):
            return fam.rifit
    return None


def corresponding_dualfit(name):
    basisfamily_list = load_basis_families()

    for fam in basisfamily_list:
        if sanitize_basisname(fam.ornate) == sanitize_basisname(name):
            return fam.dualfit
    return None
