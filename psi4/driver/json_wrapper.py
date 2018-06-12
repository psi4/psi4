#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2018 The Psi4 Developers.
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

import os
import copy
import uuid
import json

import psi4
from psi4.driver import driver
from psi4.driver import molutil
from psi4.driver import p4util
from psi4 import core

methods_dict_ = {
    'energy': driver.energy,
    'gradient': driver.gradient,
    'properties': driver.properties,
    'hessian': driver.hessian,
    'frequency': driver.frequency,
}
can_do_properties_ = {
    "dipole", "quadrupole", "mulliken_charges", "lowdin_charges", "wiberg_lowdin_indices", "mayer_indices"
}

def _clean_psi_environ(do_clean):
    if do_clean:
        psi4.core.clean_variables()
        psi4.core.clean_options()
        core.clean()

def run_json(json_data, clean=True):

    # Set scratch
    if "scratch_location" in json_data:
        psi4_io = core.IOManager.shared_object()
        psi4_io.set_default_path(json_data["scratch_location"])

    # Direct output
    outfile = os.path.join(core.IOManager.shared_object().get_default_path(), str(uuid.uuid4()) + ".json_out")
    core.set_output_file(outfile, False)

    # Set memory
    if "memory" in json_data:
        psi4.set_memory(json_data["memory"])


    # Do we return the output?
    return_output = json_data.pop("return_output", False)
    if return_output:
        json_data["raw_output"] = "Not yet run."

    # Set a few flags
    json_data["raw_output"] = None
    json_data["success"] = False

    # Attempt to run the computer
    try:
        if json_data.get("schema_name", "").startswith("qc_schema"):
            # qc_schema should be copied
            json_data = run_json_qc_schema(copy.deepcopy(json_data), clean)
        else:
            # Original run updates inplace
            run_json_original_v1_1(json_data, clean)

    except Exception as error:
        json_data["error"] = repr(error)
        json_data["success"] = False

    if return_output:
        with open(outfile, 'r') as f:
            json_data["raw_output"] = f.read()
        os.unlink(outfile)

    return json_data


def run_json_qc_schema(json_data, clean):
    """
    An implementation of the QC JSON Schema (molssi-qc-schema.readthedocs.io/en/latest/index.html#) implementation in Psi4.


    Parameters
    ----------
    json_data : JSON
        Please see molssi-qc-schema.readthedocs.io/en/latest/spec_components.html for further details.

    Notes
    -----
    !Warning! This function is experimental and likely to change in the future.
    Please report any suggestions or uses of this function on github.com/MolSSI/QC_JSON_Schema.

    Examples
    --------

    """

    # Clean a few things
    _clean_psi_environ(clean)

    # This is currently a forced override
    if json_data["schema_name"] != "qc_schema_input":
        raise KeyError("Schema name of '{}' not understood".format(json_data["schema_name"]))

    if json_data["schema_version"] != 1:
        raise KeyError("Schema version of '{}' not understood".format(json_data["schema_version"]))

    if json_data.get("nthreads", False) is not False:
        psi4.set_num_threads(json_data["nthreads"], quiet=True)

    json_data["provenance"] = {"creator": "Psi4", "version": psi4.__version__, "routine": "psi4.json.run_json"}

    # Build molecule
    mol = core.Molecule.from_schema(json_data)

    # Update molecule geometry as we orient and fix_com
    json_data["molecule"]["geometry"] = mol.geometry().np.ravel().tolist()


    # Set options
    for k, v in json_data["keywords"].items():
        core.set_global_option(k, v)

    # Setup the computation
    method = json_data["model"]["method"]
    core.set_global_option("BASIS", json_data["model"]["basis"])
    kwargs = {"return_wfn": True, "molecule": mol}

    # Handle special properties case
    if json_data["driver"] == "properties":
        if "properties" in json_data["model"]:
            kwargs["properties"] = [x.lower() for x in json_data["model"]["properties"]]

            extra = set(kwargs["properties"]) - can_do_properties_
            if len(extra):
                raise KeyError("Did not understand property key %s." % kwargs["properties"])
        else:
            kwargs["properties"] = list(can_do_properties_)


    # Actual driver run
    val, wfn = methods_dict_[json_data["driver"]](method, **kwargs)

    # Pull out a standard set of SCF properties
    psi_props = psi4.core.get_variables()
    json_data["psi4:qcvars"] = psi_props

    # Handle the return result
    if json_data["driver"] == "energy":
        json_data["return_result"] = val
    elif json_data["driver"] in ["gradient", "hessian"]:
        json_data["return_result"] = val.np.ravel().tolist()
    elif json_data["driver"] == "properties":
        ret = {}
        mtd = json_data["model"]["method"].upper()
        if "dipole" in kwargs["properties"]:
            ret["dipole"] = [psi_props[mtd + " DIPOLE " + x] for x in ["X", "Y", "Z"]]
        if "quadrupole" in kwargs["properties"]:
            ret["quadrupole"] = [psi_props[mtd + " QUADRUPOLE " + x] for x in ["XX", "XY", "XZ", "YY", "YZ", "ZZ"]]
        if "mulliken_charges" in kwargs["properties"]:
            ret["mulliken_charges"] = wfn.get_array("MULLIKEN_CHARGES").np.ravel().tolist()
        if "lowdin_charges" in kwargs["properties"]:
            ret["lowdin_charges"] = wfn.get_array("LOWDIN_CHARGES").np.ravel().tolist()
        if "wiberg_lowdin_indices" in kwargs["properties"]:
            ret["wiberg_lowdin_indices"] = wfn.get_array("WIBERG_LOWDIN_INDICES").np.ravel().tolist()
        if "mayer_indices" in kwargs["properties"]:
            ret["mayer_indices"] = wfn.get_array("MAYER_INDICES").np.ravel().tolist()

        json_data["return_result"] = ret
    else:
        raise KeyError("Did not understand Driver key %s." % json_data["driver"])


    props = {
        "calcinfo_nbasis": wfn.nso(),
        "calcinfo_nmo": wfn.nmo(),
        "calcinfo_nalpha": wfn.nalpha(),
        "calcinfo_nbeta": wfn.nbeta(),
        "calcinfo_natom": mol.geometry().shape[0],
        "scf_one_electron_energy": psi_props["ONE-ELECTRON ENERGY"],
        "scf_two_electron_energy": psi_props["TWO-ELECTRON ENERGY"],
        "nuclear_repulsion_energy": psi_props["NUCLEAR REPULSION ENERGY"],
        "scf_dipole_moment": [psi_props[x] for x in ["SCF DIPOLE X", "SCF DIPOLE Y", "SCF DIPOLE Z"]],
        "scf_iterations": int(psi_props["SCF ITERATIONS"]),
        "scf_total_energy": psi_props["SCF TOTAL ENERGY"],
        "return_energy": psi_props["CURRENT ENERGY"],
    }

    # Pull out optional SCF keywords
    other_scf = [("DFT VV10 ENERGY", "scf_vv10_energy"), ("DFT XC ENERGY", "scf_xc_energy"),
                 ("DISPERSION CORRECTION ENERGY", "scf_dispersion_correction_energy")]
    for pkey, skey in other_scf:
        if (pkey in psi_props) and (psi_props[pkey] != 0):
            props[skey] = psi_props[pkey]

    # Write out MP2 keywords
    if "MP2 CORRELATION ENERGY" in psi_props:
        props["mp2_same_spin_correlation_energy"] = psi_props["MP2 SAME-SPIN CORRELATION ENERGY"]
        props["mp2_opposite_spin_correlation_energy"] = psi_props["MP2 OPPOSITE-SPIN CORRELATION ENERGY"]
        props["mp2_singles_energy"] = 0.0
        props["mp2_doubles_energy"] = psi_props["MP2 CORRELATION ENERGY"]
        props["mp2_total_correlation_energy"] = psi_props["MP2 CORRELATION ENERGY"]
        props["mp2_total_energy"] = psi_props["MP2 TOTAL ENERGY"]

    json_data["properties"] = props
    json_data["success"] = True

    # Reset state
    _clean_psi_environ(clean)

    json_data["schema_name"] = "qc_schema_output"

    return json_data


def run_json_original_v1_1(json_data, clean):
    """
    Runs and updates the input JSON data.

    This was a trial specification introduced in Psi4 v1.1. This will be deprecated in Psi4 v1.3 in favour of the QC JSON format
    found here: http://molssi-qc-schema.readthedocs.io/en/latest/index.html#

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

    # Clean a few things
    _clean_psi_environ(clean)

    # Set a few variables
    json_data["error"] = ""
    json_data["success"] = False
    json_data["raw_output"] = None
    json_data[
        "warning"] = "Warning! This format will be deprecated in Psi4 v1.3. Please switch the standard QC JSON format."

    # Add the provenance data
    prov = {}
    prov["version"] = psi4.__version__
    prov["routine"] = "psi4.run_json"
    prov["creator"] = "Psi4"
    json_data["provenance"] = prov

    # Check input
    for check in ["driver", "method", "molecule"]:
        if check not in list(json_data):
            json_data["error"] = "Minimum input requires the %s field." % check
            return json_data

    # Check driver call
    if json_data["driver"] not in list(methods_dict_):
        json_data["error"] = "Driver parameters '%s' not recognized" % str(json_data["driver"])
        return json_data

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
    else:
        kwargs = {}

    kwargs["molecule"] = molutil.geometry(json_data["molecule"])

    # Full driver call
    kwargs["return_wfn"] = True

    if json_data["driver"] not in methods_dict_:
        raise KeyError("Driver type '%s' not recognized." % str(json_data["driver"]))

    val, wfn = methods_dict_[json_data["driver"]](*args, **kwargs)

    if isinstance(val, (float, int)):
        json_data["return_value"] = val
    elif isinstance(val, (core.Matrix, core.Vector)):
        json_data["return_value"] = val.to_serial()
    else:
        raise TypeError("Unrecognized return value of type %s\n" % type(val))

    json_data["variables"] = core.get_variables()
    json_data["success"] = True

    # Reset state
    _clean_psi_environ(clean)

    return json_data
