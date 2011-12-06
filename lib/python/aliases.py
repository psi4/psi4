import PsiMod
import re
import os
import input
import math
import warnings
from driver import *
from wrappers import *
from molecule import *
from text import *

# Place in this file quickly defined procedures such as
#   (1) aliases for complex methods
#   (2) simple modifications to existing methods


def sherrillgroup_gold_standard(**kwargs):

    if not (kwargs.has_key('func_cbs')):  kwargs['func_cbs']                 = energy

    if not (kwargs.has_key('scf_basis')):  kwargs['scf_basis']               = 'aug-cc-pVQZ'
    if not (kwargs.has_key('scf_scheme')):  kwargs['scf_scheme']             = highest_1

    if not (kwargs.has_key('name')):  kwargs['name']                         = 'mp2'
    if not (kwargs.has_key('corl_basis')):  kwargs['corl_basis']             = 'aug-cc-pV[TQ]Z'
    if not (kwargs.has_key('corl_scheme')):  kwargs['corl_scheme']           = corl_xtpl_helgaker_2

    if not (kwargs.has_key('delta_wfn')):  kwargs['delta_wfn']               = 'ccsd(t)'
    if not (kwargs.has_key('delta_wfn_lesser')):  kwargs['delta_wfn_lesser'] = 'mp2'
    if not (kwargs.has_key('delta_basis')):  kwargs['delta_basis']           = 'aug-cc-pVTZ'
    if not (kwargs.has_key('delta_scheme')):  kwargs['delta_scheme']         = highest_1

    return cbs(**kwargs)

