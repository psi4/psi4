#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2023 The Psi4 Developers.
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
import re

"""
Superfunctional builder function & handlers.
The new definition of functionals is based on a dictionary with the following structure
dict = {
           "name":  "",       name of the functional - matched against name.lower() in method lookup

          "alias":  [""],     alternative names for the method in lookup functions, processed with .lower()

       "citation":  "",       citation of the method in the standard indented format, printed in output

    "description":  "",       description of the method, printed in output

 "xc_functionals":  {         definition of a full XC functional from LibXC
      "XC_METHOD_NAME": {}      must match a LibXC method, see libxc_functionals.py for examples
     },                         if present, the x/c_functionals and x_hf/c_mp2 parameters are not read!

  "x_functionals":  {          definition of X contributions
       "X_METHOD_NAME":  {       must match a LibXC method
                 "alpha": 1.0,   coefficient for (global) GGA exchange, by default 1.0
                 "omega": 0.0,   range-separation parameter
             "use_libxc": False  whether "x_hf" parameters should be set from LibXC values for this method
                 "tweak": {},    tweak the underlying functional
     },

           "x_hf":  {          definition of HF exchange for hybrid functionals
               "alpha": 0.0,             coefficient for (global) HF exchange, by default none
                "beta": 0.0,             coefficient for long range HF exchange
               "omega": 0.0,             range separation parameters
           "use_libxc": "X_METHOD_NAME"  reads the above 3 values from specified X functional
     },

  "c_functionals":  {          definition of C contributions
       "C_METHOD_NAME":  {       must match a LibXC method
                 "alpha": 1.0,   coefficient for (global) GGA correlation, by default 1.0
                 "tweak": {},    tweak the underlying functional
    },

          "c_mp2":  {          definition of MP2 correlation double hybrid functionals
               "alpha": 0.0,     coefficient for MP2 correlation, by default none
                  "ss": 0.0,     coefficient for same spin correlation in SCS methods, forces alpha = 1.0
                  "os": 0.0,     coefficient for opposite spin correlation in SCS methods, forces alpha = 1.0
    },

     "dispersion":  {          definition of dispersion corrections
               "type": "",       dispersion type - "d2", "d3zero", "d3bj" etc., see empirical_dispersion.py
             "params": {},       parameters for the dispersion correction
                "nlc": False     (optional) logical switch to turn off nlc (e.g. VV10) correction defined by LibXC
           "citation": "",       special reference for the dispersion correction parameters, appended to output
                                   (if defined in driver, not if defined in input file)
    },
}
"""
import collections
import copy

from qcengine.programs.empirical_dispersion_resources import dashcoeff, get_dispersion_aliases

from psi4 import core

from ...p4util.exceptions import ValidationError
from . import dh_functionals, gga_functionals, hyb_functionals, lda_functionals, libxc_functionals, mgga_functionals

dict_functionals = {}
dict_functionals.update(libxc_functionals.functional_list)
dict_functionals.update(lda_functionals.functional_list)
dict_functionals.update(gga_functionals.functional_list)
dict_functionals.update(mgga_functionals.functional_list)
dict_functionals.update(hyb_functionals.functional_list)
dict_functionals.update(dh_functionals.functional_list)


def get_functional_aliases(functional_dict):
    if "alias" in functional_dict:
        aliases = [each.lower() for each in functional_dict["alias"]]
        aliases.append(functional_dict["name"].lower())
    else:
        aliases = [functional_dict["name"].lower()]
    return aliases


_dispersion_aliases = get_dispersion_aliases()

functionals = {}
dashcoeff_supplement = collections.defaultdict(lambda: collections.defaultdict(dict))

for functional_name in dict_functionals:
    functional_aliases = get_functional_aliases(dict_functionals[functional_name])

    # first create copies for aliases of parent functional
    for alias in functional_aliases:
        functionals[alias] = dict_functionals[functional_name]

    # if the parent functional is already dispersion corrected:
    if "dispersion" in dict_functionals[functional_name]:
        disp = dict_functionals[functional_name]['dispersion']
        for formal in functional_aliases:
            # "bless" the original functional dft/*_functionals dispersion definition including aliases
            dashcoeff_supplement[disp['type']]['definitions'][formal] = disp
            # add omitted default parameters of the dispersion correction
            for p,val in dashcoeff[disp['type']]['default'].items():
                if p not in dashcoeff_supplement[disp['type']]['definitions'][formal]['params'].keys():
                    dashcoeff_supplement[disp['type']]['definitions'][formal]['params'][p]=val
            # generate dispersion aliases for every functional alias
            for nominal_dispersion_level, resolved_dispersion_level in _dispersion_aliases.items():
                if resolved_dispersion_level == disp["type"]:
                    alias = formal.replace(disp["type"], nominal_dispersion_level.lower())
                    if alias not in functionals:
                        dashcoeff_supplement[disp['type']]['definitions'][formal] = disp
                        functionals[alias] = dict_functionals[functional_name]
        continue

    # else loop through dispersion types in empirical_dispersion_resources.dashcoeff (also considering aliases)
    #   and build dispersion corrected version (applies also for aliases)
    for nominal_dispersion_level, resolved_dispersion_level in _dispersion_aliases.items():
        # first check whether there is a pre-defined dispersion corrected functional
        # of the same resolved_dispersion_level type
        if functional_name + "-" + resolved_dispersion_level in dict_functionals:
            formal = functional_name + "-" + resolved_dispersion_level
            for alias in functional_aliases:
                alias += "-" + nominal_dispersion_level.lower()
                if alias not in functionals:
                    functionals[alias] = dict_functionals[formal]
        # if not, build it from dashparam logic if possible
        else:
            for dispersion_functional in dashcoeff[resolved_dispersion_level]['definitions']:
                if dispersion_functional.lower() in functional_aliases:
                    func = copy.deepcopy(dict_functionals[functional_name])
                    func["name"] += "-" + resolved_dispersion_level
                    func["dispersion"] = copy.deepcopy(
                        dashcoeff[resolved_dispersion_level]['definitions'][dispersion_functional])
                    func["dispersion"]["type"] = resolved_dispersion_level

                    # this ensures that M06-2X-D3, M06-2X-D3ZERO, M062X-D3 or M062X-D3ZERO
                    #   all point to the same method (M06-2X-D3ZERO)
                    for alias in functional_aliases:
                        alias += "-" + nominal_dispersion_level.lower()
                        functionals[alias] = func


def check_consistency(func_dictionary):
    """
    This checks the consistency of the definitions of exchange and correlation components
    of the functional, including detecting duplicate requests for LibXC params, inconsistent
    requests for HF exchange and missing correlation. It also makes sure that names of methods
    passed in using dft_functional={} syntax have a non-implemented name.
    """
    # 0a) make sure method name is set:
    if "name" not in func_dictionary:
        raise ValidationError("SCF: No method name was specified in functional dictionary.")
    else:
        name = func_dictionary["name"]
        # 0b) make sure provided name is unique:
        if (name.lower() in functionals.keys()) and (func_dictionary not in functionals.values()):
            raise ValidationError(
                "SCF: Provided name for a custom dft_functional matches an already defined one: %s." % (name))

    # 1a) sanity checks definition of xc_functionals
    if "xc_functionals" in func_dictionary:
        if "x_functionals" in func_dictionary or "x_hf" in func_dictionary:
            raise ValidationError("SCF: Duplicate specification of exchange (XC + X) in functional %s." % (name))
        elif "c_functionals" in func_dictionary or "c_mp2" in func_dictionary:
            raise ValidationError("SCF: Duplicate specification of correlation (XC + C) in functional %s." % (name))

    # 1b) require at least an empty exchange functional entry or X_HF
    elif "x_functionals" not in func_dictionary and "x_hf" not in func_dictionary:
        raise ValidationError("SCF: No exchange specified in functional %s." % (name))

    # 1c) require at least an empty correlation functional entry or C_MP2
    elif "c_functionals" not in func_dictionary and "c_mp2" not in func_dictionary:
        raise ValidationError("SCF: No correlation specified in functional %s." % (name))

    # 2) use_libxc handling:
    use_libxc = 0
    if "x_functionals" in func_dictionary:
        for item in func_dictionary["x_functionals"]:
            if "use_libxc" in func_dictionary["x_functionals"][item] and \
            func_dictionary["x_functionals"][item]["use_libxc"]:
                use_libxc += 1

    # 2a) only 1 component in x_functionals can have "use_libxc": True to prevent libxc conflicts
    if use_libxc > 1:
        raise ValidationError("SCF: Duplicate request for libxc exchange parameters in functional %s." % (name))

    # 2b) if "use_libxc" is defined in x_functionals, there shouldn't be an "x_hf" key
    elif use_libxc == 1 and "x_hf" in func_dictionary:
        raise ValidationError("SCF: Inconsistent definition of exchange in functional %s." % (name))

    # 2c) ensure libxc params requested in "x_hf" are for a functional that is included in "x_functionals"
    elif "x_hf" in func_dictionary and "use_libxc" in func_dictionary["x_hf"] \
    and func_dictionary["x_hf"]["use_libxc"] not in func_dictionary["x_functionals"]:
        raise ValidationError(
            "SCF: Libxc parameters requested for an exchange functional not defined as a component of %s." % (name))

    # 3) checks would be caught at runtime or involve only formatting.
    #    included here to preempt driver definition problems, if specific fctl not in tests.
    # 3a) check formatting for citation
    if "citation" in func_dictionary:
        cit = func_dictionary["citation"]
        if cit and not (cit.startswith('    ') and cit.endswith('\n')):
            raise ValidationError(
                f"SCF: All citations should have the form '    A. Student, B. Prof, J. Goodstuff Vol, Page, Year\n', not : {cit}"
            )
    if "dispersion" in func_dictionary:
        disp = func_dictionary["dispersion"]
        # 3b) check dispersion type present and known
        if "type" not in disp or disp["type"] not in _dispersion_aliases:
            raise ValidationError(
                f"SCF: Dispersion type ({disp['type']}) should be among ({_dispersion_aliases.keys()})")
    # 3c) check dispersion params complete
        allowed_params = sorted(dashcoeff[_dispersion_aliases[disp["type"]]]["default"].keys())
        if "params" not in disp or sorted(disp["params"].keys()) != allowed_params:
            raise ValidationError(
                f"SCF: Dispersion params for {name} ({list(disp['params'].keys())}) must include all ({allowed_params})")
    # 3d) check formatting for dispersion citation
        if "citation" in disp:
            cit = disp["citation"]
            if cit and not ((cit.startswith('    ') and cit.endswith('\n')) or re.match(r"^10.\d{4,9}/[-._;()/:A-Z0-9]+$", cit)):
                raise ValidationError(
                    f"SCF: All citations should have the form '    A. Student, B. Prof, J. Goodstuff Vol, Page, Year\n', not : {cit}"
                )


def build_superfunctional_from_dictionary(func_dictionary, npoints, deriv, restricted):
    """
    This returns a (core.SuperFunctional, dispersion) tuple based on the requested name.
    The npoints, deriv and restricted parameters are also respected.
    """

    # Sanity check first, raises ValidationError if something is wrong
    check_consistency(func_dictionary)

    # Either process the "xc_functionals" special case
    if "xc_functionals" in func_dictionary:
        for xc_key in func_dictionary["xc_functionals"]:
            xc_name = ("XC_" + xc_key).upper()
            xc_params = func_dictionary["xc_functionals"][xc_key]
            if "tweak" in xc_params:
                sup = core.SuperFunctional.XC_build(xc_name, restricted, xc_params["tweak"])
            else:
                sup = core.SuperFunctional.XC_build(xc_name, restricted)
        descr = "    " + func_dictionary["name"] + " "
        if sup.is_gga():
            if sup.x_alpha() > 0:
                descr += "Hyb-GGA "
            else:
                descr += "GGA "
        descr += "Exchange-Correlation Functional\n"
        sup.set_description(descr)

    # or combine X and C contributions into a blank SuperFunctional
    else:
        sup = core.SuperFunctional.blank()
        descr = []
        citation = []

        # Exchange processing - first the GGA part:
        # LibXC uses capital labels for the CAM coefficients, by default we're not using LibXC parms
        x_HF = {"ALPHA": 0.0, "OMEGA": 0.0, "BETA": 0.0, "used": False}
        if "x_functionals" in func_dictionary:
            x_funcs = func_dictionary["x_functionals"]
            for x_key in x_funcs:

                # Lookup the functional in LibXC
                x_name = ("XC_" + x_key).upper()
                x_func = core.LibXCFunctional(x_name, restricted)
                x_params = x_funcs[x_key]

                # If we're told to use libxc parameters for x_hf from this GGA, do so and set flag
                if "use_libxc" in x_params and x_params["use_libxc"]:
                    x_HF.update(x_func.query_libxc("XC_HYB_CAM_COEF"))
                    x_HF["used"] = True
                    x_func.set_alpha(1.0)

                if "tweak" in x_params:
                    x_func.set_tweak(x_params["tweak"], quiet=True)
                if "alpha" in x_params:
                    x_func.set_alpha(x_params["alpha"])
                if "omega" in x_params:
                    x_func.set_omega(x_params["omega"])
                sup.add_x_functional(x_func)

                # This ensures there is at least some citation for the method
                if x_func.citation() not in citation:
                    citation.append(x_func.citation())
                if x_func.description() not in descr:
                    descr.append(x_func.description())

        # Exchange processing - HF part:
        # x_HF contains zeroes or "use_libxc" params from a GGA above
        if "x_hf" in func_dictionary:
            x_params = func_dictionary["x_hf"]

            # if "use_libxc" specified here, fetch parameters and set flag
            # Duplicate definition of "use_libxc" caught in check_consistency.
            if "use_libxc" in x_params:
                x_name = ("XC_" + x_params["use_libxc"]).upper()
                x_HF.update(core.LibXCFunctional(x_name, restricted).query_libxc("XC_HYB_CAM_COEF"))
                x_HF["used"] = True

            if "alpha" in x_params:
                sup.set_x_alpha(x_params["alpha"])
            else:
                x_params["alpha"] = 0.0
            if "beta" in x_params:
                sup.set_x_beta(x_params["beta"])
            if "omega" in x_params:
                sup.set_x_omega(x_params["omega"])

        # Set LibXC parameters if requested above.
        # LibXC uses different nomenclature:
        # we need to shuffle the long and short range contributions around
        # by default, all 3 are 0.0 - different values are set only if "use_libxc" is specified
        if x_HF["used"]:
            sup.set_x_alpha(x_HF["ALPHA"])
            sup.set_x_beta(x_HF["BETA"])
            sup.set_x_omega(x_HF["OMEGA"])

        # Correlation processing - GGA part, generally same as above.
        if "c_functionals" in func_dictionary:
            c_funcs = func_dictionary["c_functionals"]
            for c_key in c_funcs:
                c_name = ("XC_" + c_key).upper()
                c_func = core.LibXCFunctional(c_name, restricted)
                c_params = func_dictionary["c_functionals"][c_key]
                if "tweak" in c_params:
                    c_func.set_tweak(c_params["tweak"], quiet=True)
                if "alpha" in c_params:
                    c_func.set_alpha(c_params["alpha"])
                else:
                    c_func.set_alpha(1.0)
                sup.add_c_functional(c_func)
                if c_func.citation() not in citation:
                    citation.append(c_func.citation())
                if c_func.description() not in descr:
                    descr.append(c_func.description())

        # Correlation processing - MP2 part
        if "c_mp2" in func_dictionary:
            c_params = func_dictionary["c_mp2"]
            if "alpha" in c_params:
                sup.set_c_alpha(c_params["alpha"])
            else:
                sup.set_c_alpha(0.0)

            # The value of alpha is locked to 1.0 C++-side when SCS is detected
            if "ss" in c_params:
                sup.set_c_ss_alpha(c_params["ss"])
                sup.set_c_alpha(1.0)
            if "os" in c_params:
                sup.set_c_os_alpha(c_params["os"])
                sup.set_c_alpha(1.0)

        # Merge descriptions and citations from above components obtained from LibXC as a fallback.
        descr = "\n".join(descr)
        citation = "\n".join(citation)
        sup.set_citation(citation)
        sup.set_description(descr)

    # Here, the above joining is usually overwritten by the proper reference.
    if "citation" in func_dictionary:
        sup.set_citation(func_dictionary["citation"])
    if "description" in func_dictionary:
        if "doi" in func_dictionary:
            sup.set_description(func_dictionary["description"].replace("\n", "") + "  (" + func_dictionary["doi"].lstrip() + ")")
        else:
            sup.set_description(func_dictionary["description"])

    # Dispersion handling for tuple assembly
    dispersion = False
    if "dispersion" in func_dictionary:
        d_params = func_dictionary["dispersion"]
        if "citation" not in d_params:
            d_params["citation"] = False
        if "nlc" in d_params:
            sup.set_vv10_b(-1.0)
            sup.set_do_vv10(d_params["nlc"])
        if d_params["type"] == 'nl':
            sup.set_vv10_b(d_params["params"]["b"])
            sup.set_vv10_c(d_params["params"]["c"])
        dispersion = d_params

    sup.set_xclib_description(core.LibXCFunctional.xclib_description())
    sup.set_max_points(npoints)
    sup.set_deriv(deriv)
    sup.set_name(func_dictionary["name"].upper())
    sup.allocate()
    return (sup, dispersion)
