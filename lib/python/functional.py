"""
Module to provide lightweight definitions of functionals and 
SuperFunctionals
"""
import PsiMod
import re
import os
import sys
from psiexceptions import *

## ==> Functionals <== ##

def build_s_x_functional(name):

    # Call this first
    fun = PsiMod.Functional.build_base('S_X')

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    fun.set_name('S_X')
    # Tab in, trailing newlines 
    fun.set_description('    Slater LSDA Exchange\n')
    # Tab in, trailing newlines 
    fun.set_citation('    J.C. Slater, Phys. Rev., 81(3):385-390, 1951\n')
    
    # These should be set by build_base, but prove that you know what's up
    fun.set_gga(False)
    fun.set_meta(False)
    fun.set_alpha(1.0)
    fun.set_omega(0.0)

    # Custom parameters

    # => End User-Customization <= #

    return fun

# Functional lookup table
functionals = {
        's_x' : build_s_x_functional
    }

def build_functional(alias):
    name = alias.lower()
    return functionals[name](name) 

def functional_list():
    val = []
    for key in functionals.keys():
        val.append(functionals[key](key))
    return val

## ==> SuperFunctionals <== ##

def build_s_x_superfunctional(name, npoints, deriv):

    # Call this first
    sup = PsiMod.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('S_X')
    # Tab in, trailing newlines 
    sup.set_description('    Slater LSDA Exchange\n')
    # Tab in, trailing newlines 
    sup.set_citation('    J.C. Slater, Phys. Rev., 81(3):385-390, 1951\n')

    # Add member functionals
    sup.add_x_functional(build_functional('S_X'))

    # Set GKS up after adding functionals
    sup.set_x_omega(0.0)
    sup.set_c_omega(0.0)
    sup.set_x_alpha(0.0)
    sup.set_c_alpha(0.0)

    # => End User-Customization <= #

    # Call this last
    sup.allocate()
    return sup 

# Superfunctional lookup table
superfunctionals = {
        's_x' : build_s_x_superfunctional
    }

def build_superfunctional(alias, npoints, deriv):
    name = alias.lower()
    return superfunctionals[name](name, npoints, deriv) 

def superfunctional_list():
    val = []
    for key in superfunctionals.keys():
        val.append(superfunctionals[key](key,1,1))
    return val
