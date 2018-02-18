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
# This file is part of Psi4.
#
# Psi4 is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3.
#
# Psi4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with Psi4; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#
"""
Builder function
"""
from psi4 import core
from psi4.driver.p4util.exceptions import *

from . import dict_lda_funcs
from . import dict_gga_funcs
from . import dict_mgga_funcs
from . import dict_hyb_funcs
from . import dict_dh_funcs

dict_functionals = {}
dict_functionals.update(dict_lda_funcs.functional_list)
dict_functionals.update(dict_gga_funcs.functional_list)
dict_functionals.update(dict_mgga_funcs.functional_list)
dict_functionals.update(dict_hyb_funcs.functional_list)
dict_functionals.update(dict_dh_funcs.functional_list)

def check_components(components):
    # sanity checks definition of x_functionals and c_functionals:
    # only 1 component can have "use_libxc": True
    use_libxc = 0
    for item in components.keys():
        if "use_libxc" in components[item].keys() and components[item]["use_libxc"]:
            use_libxc += 1
    if use_libxc > 1:
        return(False)
    # require explicit X_HF or C_MP2 with defined alpha for complex definitions
    #if len(components) > 1 and "X_HF" not in components.keys() and "C_MP2" not in components.keys():
    #    return(False)
    if "X_HF" in components.keys() and "alpha" not in components["X_HF"].keys():
        return(False)
    if "C_MP2" in components.keys() and "alpha" not in components["C_MP2"].keys():
        return(False)
    return(True)

def build_superfunctional_from_dictionary(name, npoints, deriv, restricted):
    if "xc_functionals" in dict_functionals[name].keys():
        for xc_key in dict_functionals[name]["xc_functionals"].keys():
            xc_name = "XC_" + xc_key
        sup = core.SuperFunctional.XC_build(xc_name, restricted)
        descr = "    " + name.upper() + " "
        if sup.is_gga():
            if sup.x_alpha() > 0:
                descr += "Hyb-GGA "
            else:
                descr += "GGA "
        descr += "Exchange-Correlation Functional\n"
        sup.set_description(descr)
    else:
        sup = core.SuperFunctional.blank()
        descr = []
        citation = []
        # first sanity-check input
        if check_components(dict_functionals[name]["x_functionals"]):
            x_HF = {"ALPHA": 0.0, "OMEGA": 0.0, "BETA": 0.0}
            # then process each functional
            for x_key in dict_functionals[name]["x_functionals"].keys():
                if x_key != "X_HF":
                    x_name = "XC_" + x_key
                    x_func = core.LibXCFunctional(x_name, restricted)
                    x_params = dict_functionals[name]["x_functionals"][x_key]
                    if "use_libxc" in x_params and x_params["use_libxc"]:
                        x_HF = x_func.query_libxc("XC_HYB_CAM_COEF")
                        x_func.set_alpha(1.0)
                    if "tweak" in x_params.keys():
                        x_func.set_tweak(x_params["tweak"])
                    if "alpha" in x_params.keys():
                        x_func.set_alpha(x_params["alpha"])
                    if "omega" in x_params.keys():
                        x_func.set_omega(x_params["omega"])
                    sup.add_x_functional(x_func)
                    if x_func.citation() not in citation:
                        citation.append(x_func.citation())
                    if x_func.description() not in descr:
                        descr.append(x_func.description())
            # finally set HF exchange
            if "X_HF" in dict_functionals[name]["x_functionals"].keys():
                x_params = dict_functionals[name]["x_functionals"]["X_HF"]
                
                if "alpha" in x_params.keys():
                    sup.set_x_alpha(x_params["alpha"])
                else:
                    x_params["alpha"] = 0.0
                if "beta" in x_params.keys():
                    sup.set_x_beta(x_params["beta"])
                if "omega" in x_params.keys():
                    sup.set_x_omega(x_params["omega"])
            else:
                sup.set_x_alpha(x_HF["ALPHA"] + x_HF["BETA"])
                sup.set_x_beta(x_HF["ALPHA"]-x_HF["BETA"])
                sup.set_x_omega(x_HF["OMEGA"])
        else:
            raise ValidationError("SCF: Incorrect specification of exchange in functional %s." % (name))
        if check_components(dict_functionals[name]["c_functionals"]):
            for c_key in dict_functionals[name]["c_functionals"].keys():
                if c_key != "C_MP2":
                    c_name = "XC_" + c_key
                    c_func = core.LibXCFunctional(c_name, restricted)
                    c_params = dict_functionals[name]["c_functionals"][c_key]
                    if "tweak" in c_params.keys():
                        c_func.set_tweak(c_params["tweak"])
                    if "alpha" in c_params.keys():
                        c_func.set_alpha(c_params["alpha"])
                    else:
                        c_func.set_alpha(1.0)
                    sup.add_c_functional(c_func)
                    if c_func.citation() not in citation:
                        citation.append(c_func.citation())
                    if c_func.description() not in descr:
                        descr.append(c_func.description())
            if "C_MP2" in dict_functionals[name]["c_functionals"].keys():
                c_params = dict_functionals[name]["c_functionals"]["C_MP2"]
                sup.set_c_alpha(c_params["alpha"])
                if "ss" in c_params:
                    sup.set_c_ss_alpha(c_params["ss"])
                if "os" in c_params:
                    sup.set_c_os_alpha(c_params["os"])
        else:
            raise ValidationError("SCF: Incorrect specification of correlation in functional %s." % (name))        
            
        descr = "\n".join(descr)
        citation = "\n".join(citation)
        sup.set_citation(citation)
        sup.set_description(descr)
    
    if "citation" in dict_functionals[name].keys():
        sup.set_citation(dict_functionals[name]["citation"])
    if "description" in dict_functionals[name].keys():
        sup.set_description(dict_functionals[name]["description"])
    
    dispersion = False
    if "dispersion" in dict_functionals[name].keys():
        d_params = dict_functionals[name]["dispersion"]
        if d_params["type"].lower() in ["d","d2","d2p4"]:
            d_params["type"] = "d2p4"
        elif d_params["type"].lower() in ["d3","d3(0)","d3zero"]:
            d_params["type"] = "d3zero"
        elif d_params["type"].lower() in ["d3bj","d3(bj)"]:
            d_params["type"] = "d3bj"
        elif d_params["type"].lower() in ["d3m","d3m(0)","d3m0","d3mzero"]:
            d_params["type"] = "d3mzero"
        elif d_params["type"].lower() in ["d3mbj","d3m(bj)"]:
            d_params["type"] = "d3mbj"
        elif d_params["type"].lower() in ["chg", "das2009", "das2010"]:
            pass
        else:
            raise ValidationError("SCF: Impossible to classify dispersion type %s." % (d_params["type"]))
        if "citation" not in d_params.keys():
            d_params["citation"] = False
        dispersion = ("custom", d_params["type"], d_params["params"], d_params["citation"])
        
    
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)
    sup.set_name(name.upper())
    sup.allocate()
    return (sup, dispersion)
