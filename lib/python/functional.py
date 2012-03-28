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

def build_lyp_c_functional(name):

    # Call this first
    fun = PsiMod.Functional.build_base('LYP_C')

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    fun.set_name('LYP_C')
    # Tab in, trailing newlines 
    fun.set_description('    LYP GGA Correlation Functional\n')
    # Tab in, trailing newlines 
    fun.set_citation('    B. Miehlich et. al., Chem. Phys. Lett., 157(3), 200-206 1989\n')
    
    # These should be set by build_base, but prove that you know what's up
    fun.set_gga(True)
    fun.set_meta(False)
    fun.set_alpha(1.0)
    fun.set_omega(0.0)

    # Custom parameters
    # Always built-in for this functional

    # => End User-Customization <= #

    return fun

def build_primitive_functional(name):

    # Call this first
    key = name.upper()
    if (key[0] == 'W'):
        key[0] = 'w'
    fun = PsiMod.Functional.build_base(key)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    fun.set_name(key)
    # Tab in, trailing newlines 
    fun.set_description(fun.description())
    # Tab in, trailing newlines 
    fun.set_citation(fun.citation())
    
    # These should be set by build_base, but prove that you know what's up
    fun.set_gga(fun.is_gga())
    fun.set_meta(fun.is_meta())
    fun.set_alpha(fun.alpha())
    fun.set_omega(fun.omega())

    # Custom parameters
    # Always built-in for this functional

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

def build_lyp_c_superfunctional(name, npoints, deriv):

    # Call this first
    sup = PsiMod.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('LYP_C')
    # Tab in, trailing newlines 
    sup.set_description('    LYP GGA Correlation Functional\n')
    # Tab in, trailing newlines 
    sup.set_citation('    B. Miehlich et. al., Chem. Phys. Lett., 157(3), 200-206 1989\n')

    # Add member functionals
    sup.add_c_functional(build_functional('LYP_C'))

    # Set GKS up after adding functionals
    sup.set_x_omega(0.0)
    sup.set_c_omega(0.0)
    sup.set_x_alpha(0.0)
    sup.set_c_alpha(0.0)

    # => End User-Customization <= #

    # Call this last
    sup.allocate()
    return sup 

def build_blyp_superfunctional(name, npoints, deriv):

    # Call this first
    sup = PsiMod.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    # No spaces, keep it short and according to convention
    sup.set_name('BLYP')
    # Tab in, trailing newlines 
    sup.set_description('    BLYP GGA Exchange-Correlation Functional\n')
    # Tab in, trailing newlines 
    sup.set_citation('    P.J. Stephens et. al., J. Phys. Chem., 98, 11623-11627, 1994\n    B. Miehlich et. al., Chem. Phys. Lett., 157(3), 200-206 1989\n')

    # Add member functionals
    sup.add_x_functional(build_functional('B88_X'))
    sup.add_c_functional(build_functional('LYP_C'))

    # Set GKS up after adding functionals
    sup.set_x_omega(0.0)
    sup.set_c_omega(0.0)
    sup.set_x_alpha(0.0)
    sup.set_c_alpha(0.0)

    # => End User-Customization <= #

    # Call this last
    sup.allocate()
    return sup 

def build_primitive_superfunctional(name, npoints, deriv):

    # Call this first
    sup = PsiMod.SuperFunctional.blank()
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)

    # => User-Customization <= #

    key = name.upper()
    if (key[0] == 'W'):
        key[0] = 'w'
    fun = build_functional(key)

    # No spaces, keep it short and according to convention
    sup.set_name(key)
    # Tab in, trailing newlines 
    sup.set_description(fun.description())
    # Tab in, trailing newlines 
    sup.set_citation(fun.citation())

    # Add member functionals
    
    if (key[-1] == 'X'):
        sup.add_x_functional(fun)
    else:
        sup.add_c_functional(fun)

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
        's_x'    : build_s_x_superfunctional,
        'b88_x'  : build_b88_x_superfunctional,
        'b3_x'   : build_b3_x_superfunctional,
        'pbe_x'  : build_pbe_x_superfunctional,
        'pw91_x' : build_pw91_x_superfunctional,
        'ws_x'   : build_ws_x_superfunctional,
        'wpbe_x' : build_wpbe_x_superfunctional,
        'blyp'   : build_blyp_superfunctional,
    }

def build_superfunctional(alias, npoints, deriv):
    name = alias.lower()
    return superfunctionals[name](name, npoints, deriv) 

def superfunctional_list():
    val = []
    for key in superfunctionals.keys():
        val.append(superfunctionals[key](key,1,1))
    return val

def test_ccl_functional(functional, ccl_functional):

    if (not os.path.exists('data_pt_%s.html' %(ccl_functional))):
        os.system('wget ftp://ftp.dl.ac.uk/qcg/dft_library/data_pt_%s.html' % ccl_functional)
    fh = open('data_pt_%s.html' %(ccl_functional))
    lines = fh.readlines()
    fh.close()
    
    points = []
    point = {}
    
    rho_line = re.compile(r'^\s*rhoa=\s*(-?\d+\.\d+E[+-]\d+)\s*rhob=\s*(-?\d+\.\d+E[+-]\d+)\s*sigmaaa=\s*(-?\d+\.\d+E[+-]\d+)\s*sigmaab=\s*(-?\d+\.\d+E[+-]\d+)\s*sigmabb=\s*(-?\d+\.\d+E[+-]\d+)\s*')
    val_line = re.compile(r'^\s*(\w*)\s*=\s*(-?\d+\.\d+E[+-]\d+)')
    
    aliases = { 'zk'            : 'v',
                'vrhoa'         : 'v_rho_a',
                'vrhob'         : 'v_rho_b',
                'vsigmaaa'      : 'v_gamma_aa',
                'vsigmaab'      : 'v_gamma_ab',
                'vsigmabb'      : 'v_gamma_bb',
                'v2rhoa2'       : 'v_rho_a_rho_a',
                'v2rhoab'       : 'v_rho_a_rho_b',
                'v2rhob2'       : 'v_rho_b_rho_b',
                'v2rhoasigmaaa' : 'v_rho_a_gamma_aa',
                'v2rhoasigmaab' : 'v_rho_a_gamma_ab',
                'v2rhoasigmabb' : 'v_rho_a_gamma_bb',
                'v2rhobsigmaaa' : 'v_rho_b_gamma_aa',
                'v2rhobsigmaab' : 'v_rho_b_gamma_ab',
                'v2rhobsigmabb' : 'v_rho_b_gamma_bb',
                'v2sigmaaa2'    : 'v_gamma_aa_gamma_aa',
                'v2sigmaaaab'   : 'v_gamma_aa_gamma_ab',
                'v2sigmaaabb'   : 'v_gamma_aa_gamma_bb',
                'v2sigmaab2'    : 'v_gamma_ab_gamma_ab',
                'v2sigmaabbb'   : 'v_gamma_ab_gamma_bb',
                'v2sigmabb2'    : 'v_gamma_bb_gamma_bb',
              }
    
    for line in lines:
    
        mobj = re.match(rho_line, line)
        if (mobj):
    
            if len(point):
                points.append(point)
                point = {}    
    
            point['rho_a']    = float(mobj.group(1))
            point['rho_b']    = float(mobj.group(2))
            point['gamma_aa'] = float(mobj.group(3))
            point['gamma_ab'] = float(mobj.group(4))
            point['gamma_bb'] = float(mobj.group(5))
    
            continue
    
        mobj = re.match(val_line, line)
        if (mobj):
            point[aliases[mobj.group(1)]] = float(mobj.group(2))
    
    points.append(point)
            
    N = len(points)
    rho_a = PsiMod.Vector(N)
    rho_b = PsiMod.Vector(N)
    gamma_aa = PsiMod.Vector(N)
    gamma_ab = PsiMod.Vector(N)
    gamma_bb = PsiMod.Vector(N)
    tau_a = PsiMod.Vector(N)
    tau_b = PsiMod.Vector(N)
    
    index = 0;
    for point in points:
        rho_a[index] = point['rho_a']
        rho_b[index] = point['rho_b']
        gamma_aa[index] = point['gamma_aa']
        gamma_ab[index] = point['gamma_ab']
        gamma_bb[index] = point['gamma_bb']
        index = index + 1
    
    super = build_superfunctional(functional, N, 1) 
    super.test_functional(rho_a,rho_b,gamma_aa,gamma_ab,gamma_bb,tau_a,tau_b)
    
    v = super.value('V')
    v_rho_a = super.value('V_RHO_A')
    v_rho_b = super.value('V_RHO_B')
    v_gamma_aa = super.value('V_GAMMA_AA')
    v_gamma_ab = super.value('V_GAMMA_AB')
    v_gamma_bb = super.value('V_GAMMA_BB')
    
    if not v_gamma_aa:
        v_gamma_aa = tau_a
        v_gamma_ab = tau_a
        v_gamma_bb = tau_a
    
    tasks = ['v', 'v_rho_a', 'v_rho_b', 'v_gamma_aa', 'v_gamma_ab', 'v_gamma_bb']
    mapping = {
            'v' : v,
            'v_rho_a' : v_rho_a,
            'v_rho_b' : v_rho_b,
            'v_gamma_aa' : v_gamma_aa,
            'v_gamma_ab' : v_gamma_ab,
            'v_gamma_bb' : v_gamma_bb,
        }
    
    super.print_detail(3)
    index = 0;
    for point in points:
        PsiMod.print_out('rho_a= %11.3E, rho_b= %11.3E, gamma_aa= %11.3E, gamma_ab= %11.3E, gamma_bb= %11.3E\n' %(rho_a[index], rho_b[index], gamma_aa[index], gamma_ab[index], gamma_bb[index]))
    
        for task in tasks:
            v_ref = point[task]
            v_obs = mapping[task][index]
            delta = v_obs - v_ref
            if (v_ref == 0.0):
                epsilon = 0.0
            else:
                epsilon = abs(delta / v_ref)
            if (epsilon < 1.0E-11):
                passed = 'PASSED'
            else:
                passed = 'FAILED'
    
            PsiMod.print_out('\t%-15s %24.16E %24.16E %24.16E %24.16E %6s\n' % (task, v_ref, v_obs, delta, epsilon, passed))
    
        index = index + 1

    PsiMod.print_out('\n')
