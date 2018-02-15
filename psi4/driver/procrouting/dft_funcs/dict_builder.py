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

from . import dict_xc_funcs
from . import dict_gga_funcs
from . import dict_hyb_funcs
from . import dict_dh_funcs

dict_functionals = {}
dict_functionals.update(dict_xc_funcs.functional_list)
dict_functionals.update(dict_gga_funcs.functional_list)
dict_functionals.update(dict_hyb_funcs.functional_list)
dict_functionals.update(dict_dh_funcs.functional_list)

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
        x_HF_alpha = 1.0
        x_HF = {"ALPHA": 1.0, "OMEGA": 0.0, "BETA": 1.0}
        for x_key in dict_functionals[name]["x_functionals"].keys():
            if x_key != "X_HF":
                x_name = "XC_" + x_key
                x_func = core.LibXCFunctional(x_name, restricted)
                x_params = dict_functionals[name]["x_functionals"][x_key]
                x_HF = x_func.query_libxc("XC_HYB_CAM_COEF")
                x_exx_params = x_func.query_libxc("XC_HYB_EXX_COEF")
                #if x_func.is_gga():
                if "alpha" in x_params.keys():
                    x_func.set_alpha(x_params["alpha"])
                    x_HF["ALPHA"] = 1.0 - x_params["alpha"]
                else:
                    x_func.set_alpha(1.0)
                if "omega" in x_params.keys():
                    x_HF["OMEGA"] = x_params["omega"]
                if "beta" in x_params.keys():
                    x_HF["BETA"] = x_params["beta"]
                sup.add_x_functional(x_func)
                if x_func.citation() not in citation:
                    citation.append(x_func.citation())
                if x_func.description() not in descr:
                    descr.append(x_func.description())
        if "X_HF" in dict_functionals[name]["x_functionals"].keys():
            x_params = dict_functionals[name]["x_functionals"]["X_HF"]
            if "alpha" in x_params.keys():
                x_HF["ALPHA"] = x_params["alpha"]
            if "beta" in x_params.keys():
                x_HF["BETA"] = x_params["beta"]
            if "omega" in x_params.keys():
                x_HF["OMEGA"] = x_params["omega"]
        sup.set_x_alpha(x_HF["ALPHA"] + x_HF["BETA"])
        sup.set_x_beta(-x_HF["BETA"])
        sup.set_x_omega(x_HF["OMEGA"])
        c_MP2_alpha = 1.0
        for c_key in dict_functionals[name]["c_functionals"].keys():
            if c_key != "C_MP2":
                c_name = "XC_" + c_key
                c_func = core.LibXCFunctional(c_name, restricted)
                c_params = dict_functionals[name]["c_functionals"][c_key]
                c_cam_params = c_func.query_libxc("XC_HYB_CAM_COEF")
                c_exx_params = c_func.query_libxc("XC_HYB_EXX_COEF")
                if "alpha" in c_params.keys():
                    c_func.set_alpha(c_params["alpha"])
                    c_MP2_alpha = c_MP2_alpha - c_params["alpha"]
                else:
                    c_func.set_alpha(1.0)
                    c_MP2_alpha = c_cam_params["ALPHA"]
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
            sup.set_c_alpha(c_MP2_alpha)
            
        descr = "\n".join(descr)
        citation = "\n".join(citation)
        sup.set_citation(citation)
        sup.set_description(descr)
    
    if "citation" in dict_functionals[name].keys():
        sup.set_citation(dict_functionals[name]["citation"])
    if "description" in dict_functionals[name].keys():
        sup.set_description(dict_functionals[name]["description"])
    
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)
    sup.set_name(name.upper())
    sup.allocate()
    return (sup, False)
