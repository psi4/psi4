"""
Module to provide lightweight definitions of functionals and 
SuperFunctionals
"""
import PsiMod
import re
import os
import sys
import math
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

def build_b88_x_functional(name):

    # Call this first
    fun = PsiMod.Functional.build_base('B88_X')

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    fun.set_name('B88_X')
    # Tab in, trailing newlines 
    fun.set_description('    Becke88 GGA Exchange\n')
    # Tab in, trailing newlines 
    fun.set_citation('    A.D. Becke, Phys. Rev. A, 38(6):3098-3100, 1988\n')
    
    # These should be set by build_base, but prove that you know what's up
    fun.set_gga(True)
    fun.set_meta(False)
    fun.set_alpha(1.0)
    fun.set_omega(0.0)

    # Custom parameters
    fun.set_parameter('B88_d', 0.0042);
    fun.set_parameter('B88_a', 1.0000);

    # => End User-Customization <= #

    return fun

def build_b3_x_functional(name):

    # Call this first
    fun = PsiMod.Functional.build_base('B88_X')

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    fun.set_name('B3_X')
    # Tab in, trailing newlines 
    fun.set_description('    Becke88 GGA Exchange (B3LYP weighting)\n')
    # Tab in, trailing newlines 
    fun.set_citation('    P.J. Stephens et. al., J. Phys. Chem., 98, 11623-11627, 1994\n')
    
    # These should be set by build_base, but prove that you know what's up
    fun.set_gga(True)
    fun.set_meta(False)
    fun.set_alpha(1.0)
    fun.set_omega(0.0)

    # Custom parameters
    fun.set_parameter('B88_d', 0.0042);
    fun.set_parameter('B88_a', 0.9000);

    # => End User-Customization <= #

    return fun

def build_pbe_x_functional(name):

    # Call this first
    fun = PsiMod.Functional.build_base('PBE_X')

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    fun.set_name('PBE_X')
    # Tab in, trailing newlines 
    fun.set_description('    PBE GGA Exchange Hole (Parameter Free)\n')
    # Tab in, trailing newlines 
    fun.set_citation('    J.P. Perdew et. al., Phys. Rev. Lett., 77(18), 3865-3868, 1996\n')
    
    # These should be set by build_base, but prove that you know what's up
    fun.set_gga(True)
    fun.set_meta(False)
    fun.set_alpha(1.0)
    fun.set_omega(0.0)

    # Custom parameters

    # => End User-Customization <= #

    return fun

def build_pw91_x_functional(name):

    # Call this first
    fun = PsiMod.Functional.build_base('PW91_X')

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    fun.set_name('PW91_X')
    # Tab in, trailing newlines 
    fun.set_description('    PW91 Parameterized GGA Exchange\n')
    # Tab in, trailing newlines 
    fun.set_citation('    J.P. Perdew et. al., Phys. Rev. B., 46(11), 6671-6687, 1992\n')
    
    # These should be set by build_base, but prove that you know what's up
    fun.set_gga(True)
    fun.set_meta(False)
    fun.set_alpha(1.0)
    fun.set_omega(0.0)

    # Custom parameters
    k01 = math.pow(6.0 * math.pi * math.pi, 1.0/3.0);
    k02 = k01 * k01
    k04 = k02 * k02
    fun.set_parameter('PW91_a1', 0.19645 / ( 2.0 * k01));
    fun.set_parameter('PW91_a2', 7.79560 / ( 2.0 * k01));
    fun.set_parameter('PW91_a3', 0.27430 / ( 4.0 * k02));
    fun.set_parameter('PW91_a4', 0.15080 / ( 4.0 * k02));
    fun.set_parameter('PW91_a5', 100.000 / ( 4.0 * k02));
    fun.set_parameter('PW91_a6', 0.00400 / (16.0 * k04));

    # => End User-Customization <= #

    return fun

def build_ws_x_functional(name):

    # Call this first
    fun = PsiMod.Functional.build_base('wS_X')

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    fun.set_name('wS_X')
    # Tab in, trailing newlines 
    fun.set_description('    Slater Short-Range LSDA Exchange\n')
    # Tab in, trailing newlines 
    fun.set_citation('    Adamson et. al., J. Comput. Chem., 20(9), 921-927, 1999\n')
    
    # These should be set by build_base, but prove that you know what's up
    fun.set_gga(False)
    fun.set_meta(False)
    fun.set_alpha(1.0)
    fun.set_omega(0.3)

    # Custom parameters

    # => End User-Customization <= #

    return fun

def build_wpbe_x_functional(name):

    # Call this first
    fun = PsiMod.Functional.build_base('wPBE_X')

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    fun.set_name('wPBE_X')
    # Tab in, trailing newlines 
    fun.set_description('    PBE Short-Range GGA Exchange (Exact SR-Hole)\n')
    # Tab in, trailing newlines 
    fun.set_citation('    Henderson et. al., J. Chem. Phys., 128, 194105, 2008\n')
    
    # These should be set by build_base, but prove that you know what's up
    fun.set_gga(True)
    fun.set_meta(False)
    fun.set_alpha(1.0)
    fun.set_omega(0.3)

    # Custom parameters

    # => End User-Customization <= #

    return fun

# Functional lookup table
functionals = {
        's_x'         : build_s_x_functional,
        'b88_x'       : build_b88_x_functional,
        'b3_x'        : build_b3_x_functional,
        'pbe_x'       : build_pbe_x_functional,
        'pw91_x'      : build_pw91_x_functional,
        'ws_x'        : build_ws_x_functional,
        'wpbe_x'      : build_wpbe_x_functional,
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

def build_b88_x_superfunctional(name, npoints, deriv):

    # Call this first
    sup = PsiMod.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('B88_X')
    # Tab in, trailing newlines 
    sup.set_description('    Becke88 GGA Exchange\n')
    # Tab in, trailing newlines 
    sup.set_citation('    A.D. Becke, Phys. Rev. A, 38(6):3098-3100, 1988\n')

    # Add member functionals
    sup.add_x_functional(build_functional('B88_X'))

    # Set GKS up after adding functionals
    sup.set_x_omega(0.0)
    sup.set_c_omega(0.0)
    sup.set_x_alpha(0.0)
    sup.set_c_alpha(0.0)

    # => End User-Customization <= #

    # Call this last
    sup.allocate()
    return sup 

def build_b3_x_superfunctional(name, npoints, deriv):

    # Call this first
    sup = PsiMod.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('B3_X')
    # Tab in, trailing newlines 
    sup.set_description('    Becke88 GGA Exchange (B3LYP weighting)\n')
    # Tab in, trailing newlines 
    sup.set_citation('    P.J. Stephens et. al., J. Phys. Chem., 98, 11623-11627, 1994\n')

    # Add member functionals
    sup.add_x_functional(build_functional('B3_X'))

    # Set GKS up after adding functionals
    sup.set_x_omega(0.0)
    sup.set_c_omega(0.0)
    sup.set_x_alpha(0.2)
    sup.set_c_alpha(0.0)

    # => End User-Customization <= #

    # Call this last
    sup.allocate()
    return sup 

def build_pbe_x_superfunctional(name, npoints, deriv):

    # Call this first
    sup = PsiMod.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('PBE_X')
    # Tab in, trailing newlines 
    sup.set_description('    PBE GGA Exchange Hole (Parameter Free)\n')
    # Tab in, trailing newlines 
    sup.set_citation('    J.P. Perdew et. al., Phys. Rev. Lett., 77(18), 3865-3868, 1996\n')

    # Add member functionals
    sup.add_x_functional(build_functional('PBE_X'))

    # Set GKS up after adding functionals
    sup.set_x_omega(0.0)
    sup.set_c_omega(0.0)
    sup.set_x_alpha(0.0)
    sup.set_c_alpha(0.0)

    # => End User-Customization <= #

    # Call this last
    sup.allocate()
    return sup 

def build_pw91_x_superfunctional(name, npoints, deriv):

    # Call this first
    sup = PsiMod.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('PW91_X')
    # Tab in, trailing newlines 
    sup.set_description('    PW91 Parameterized GGA Exchange\n')
    # Tab in, trailing newlines 
    sup.set_citation('    J.P. Perdew et. al., Phys. Rev. B., 46(11), 6671-6687, 1992\n')

    # Add member functionals
    sup.add_x_functional(build_functional('PW91_X'))

    # Set GKS up after adding functionals
    sup.set_x_omega(0.0)
    sup.set_c_omega(0.0)
    sup.set_x_alpha(0.0)
    sup.set_c_alpha(0.0)

    # => End User-Customization <= #

    # Call this last
    sup.allocate()
    return sup 

def build_ws_x_superfunctional(name, npoints, deriv):

    # Call this first
    sup = PsiMod.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('wS_X')
    # Tab in, trailing newlines 
    sup.set_description('    Slater Short-Range LSDA Exchange\n')
    # Tab in, trailing newlines 
    sup.set_citation('    Adamson et. al., J. Comput. Chem., 20(9), 921-927, 1999\n')

    # Add member functionals
    sup.add_x_functional(build_functional('wS_X'))

    # Set GKS up after adding functionals
    sup.set_x_omega(0.3)
    sup.set_c_omega(0.0)
    sup.set_x_alpha(0.0)
    sup.set_c_alpha(0.0)

    # => End User-Customization <= #

    # Call this last
    sup.allocate()
    return sup 

def build_wpbe_x_superfunctional(name, npoints, deriv):

    # Call this first
    sup = PsiMod.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('wPBE_X')
    # Tab in, trailing newlines 
    sup.set_description('    PBE Short-Range GGA Exchange (Exact SR-Hole)\n')
    # Tab in, trailing newlines 
    sup.set_citation('    Henderson et. al., J. Chem. Phys., 128, 194105, 2008\n')

    # Add member functionals
    sup.add_x_functional(build_functional('wPBE_X'))

    # Set GKS up after adding functionals
    sup.set_x_omega(0.3)
    sup.set_c_omega(0.0)
    sup.set_x_alpha(0.0)
    sup.set_c_alpha(0.0)

    # => End User-Customization <= #

    # Call this last
    sup.allocate()
    return sup 

# Superfunctional lookup table
superfunctionals = {
        's_x'    : build_s_x_superfunctional,
        'b88_x'  : build_b88_x_superfunctional,
        'b3_x'   : build_b3_x_superfunctional,
        'pbe_x'  : build_pbe_x_superfunctional,
        'pw91_x' : build_pw91_x_superfunctional,
        'ws_x'   : build_ws_x_superfunctional,
        'wpbe_x' : build_wpbe_x_superfunctional,
    }

def build_superfunctional(alias, npoints, deriv):
    name = alias.lower()
    return superfunctionals[name](name, npoints, deriv) 

def superfunctional_list():
    val = []
    for key in superfunctionals.keys():
        val.append(superfunctionals[key](key,1,1))
    return val
