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
Runs a JSON input psi file.
"""


from psi4.driver import driver
from psi4.driver import molutil
from psi4 import core
import psi4

import json
import uuid
import copy
import os


methods_dict = {
    'energy': driver.energy,
    'gradient': driver.gradient,
    'property': driver.properties,
    'optimize': driver.optimize,
    'hessian': driver.hessian,
    'frequency': driver.frequency,
}


def run_json(json_data):
    """
    Runs and updates the input JSON data.

    Parameters
    ----------
    json_data : JSON input
        Required input fields:
            - molecule : str
                A string representation of the molecule. Any valid psi4 synatx is valid.
            - driver : str
                The driver method to use valid arguments are "energy", "gradient", "property"
            - args : str
                Input arguments to the driver function
            - kwargs : dict
                Input dictionary to the driver function

        Optional input fields:
            - kwargs : dict
                Additional kwargs to be passed to the driver.
            - options : dict
                Global options to set in the Psi4 run.
            - scratch_location : str
                Overide the default scratch location.
            - return_output : bool
                Return the full string reprsentation of the output or not.

        Output fields:
            - return_value : float, psi4.core.Vector, psi4.core.Matrix
                The return value of the input function.
            - variables : dict
                The list of Psi4 variables generated from the run.
            - success : bool
                Indicates if the run was successful or not.
            - error : str
                If an error is raised the result is returned here.
            - raw_output : str
                The full Psi4 output if requested.


    Notes
    -----
    !Warning! This function is experimental and likely to change in the future.
    Please report any suggestions or uses of this function on github.com/psi4/psi4.

    Examples
    --------

    # CP corrected Helium dimer energy in the STO-3G basis.
    >>> json_data = {}

    >>> json_data["molecule"] = "He 0 0 0\n--\nHe 0 0 1"
    >>> json_data["driver"] = "energy"
    >>> json_data["method"] = 'SCF'
    >>> json_data["kwargs"] = {"bsse_type": "cp"}
    >>> json_data["options"] = {"BASIS": "STO-3G"}

    >>> run_json(json_data)
    {
        "raw_output": "Output storing was not requested.",
        "options": {
            "BASIS": "STO-3G"
        },
        "driver": "energy",
        "molecule": "He 0 0 0\n--\nHe 0 0 1",
        "method": "SCF",
        "variables": {
            "SCF N ITERS": 2.0,
            "SCF DIPOLE Y": 0.0,
            "CURRENT DIPOLE Y": 0.0,
            "CP-CORRECTED 2-BODY INTERACTION ENERGY": 0.1839360538612116,
            "HF TOTAL ENERGY": -5.433191881443323,
            "SCF TOTAL ENERGY": -5.433191881443323,
            "TWO-ELECTRON ENERGY": 4.124089347186247,
            "SCF ITERATION ENERGY": -5.433191881443323,
            "CURRENT DIPOLE X": 0.0,
            "CURRENT DIPOLE Z": 0.0,
            "CURRENT REFERENCE ENERGY": -5.433191881443323,
            "CURRENT ENERGY": 0.1839360538612116,
            "COUNTERPOISE CORRECTED TOTAL ENERGY": -5.433191881443323,
            "SCF DIPOLE Z": 0.0,
            "COUNTERPOISE CORRECTED INTERACTION ENERGY": 0.1839360538612116,
            "NUCLEAR REPULSION ENERGY": 2.11670883436,
            "SCF DIPOLE X": 0.0,
            "ONE-ELECTRON ENERGY": -11.67399006298957
        },
        "return_value": 0.1839360538612116,
        "error": "",
        "success": true,
        "provenance": {
            "creator": "Psi4",
            "routine": "psi4.run_json",
            "version": "1.1a1"
        },
        "kwargs": {
            "bsse_type": "cp"
        }
    }
    """

    # Set a few variables
    json_data["error"] = ""
    json_data["raw_output"] = "Output storing was not requested."

    # Check input
    for check in ["driver", "method", "molecule"]:
        if check not in list(json_data):
            json_data["error"] = "Minimum input requires the %s field." % check
            return False

    if json_data["driver"] not in list(methods_dict):
        json_data["error"] = "Driver parameters '%s' not recognized" % str(json_data["driver"])
        return False


    if "scratch_location" in list(json_data):
        psi4_io = core.IOManager.shared_object()
        psi4_io.set_default_path(json_data["scratch_location"])

    # Do we return the output?
    return_output = json_data.pop("return_output", False)
    if return_output:
        outfile = str(uuid.uuid4()) + ".json_out"
        core.set_output_file(outfile, False)
        json_data["raw_output"] = "Not yet run."

    try:
        # Set options
        if "options" in json_data.keys() and (json_data["options"] is not None):
            for k, v in json_data["options"].items():
                core.set_global_option(k, v)

        # Rework args
        args = json_data["method"]
        if not isinstance(args, (list, tuple)):
            args = [args]

        # Deep copy kwargs
        if "kwargs" in list(json_data):
            kwargs = copy.deepcopy(json_data["kwargs"])
            kwargs["molecule"] = molutil.geometry(json_data["molecule"])
        else:
            kwargs = {}

        # Full driver call
        kwargs["return_wfn"] = True

        val, wfn = methods_dict[json_data["driver"]](*args, **kwargs)

        if isinstance(val, (float)):
            json_data["return_value"] = val
        elif isinstance(val, (core.Matrix, core.Vector)):
            json_data["return_value"] = val.to_serial()
        else:
            json_data["error"] += "Unrecognized return value of type %s\n" % type(val)
            json_data["return_value"] = val

        json_data["variables"] = core.get_variables()
        json_data["success"] = True

        prov = {}
        prov["version"] = psi4.__version__
        prov["routine"] = "psi4.run_json"
        prov["creator"] = "Psi4"
        json_data["provenance"] = prov

    except Exception as error:
        json_data["error"] += repr(error)
        json_data["success"] = False

    if return_output:
        with open(outfile, 'r') as f:
            json_data["raw_output"] = f.read()
        os.unlink(outfile)

    return json_data

