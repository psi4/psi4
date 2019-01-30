#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2019 The Psi4 Developers.
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

import atexit
import copy
import json
import numpy as np
import os
import uuid

import psi4
from psi4.driver import driver
from psi4.driver import molutil
from psi4.driver import p4util
from psi4 import core

## Methods and properties blocks

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

## QCSchema translation blocks

_qcschema_translation = {

    # Generics
    "generics": {
        "return_energy": {"variables": "CURRENT ENERGY"},
        "nuclear_repulsion_energy": {"variables": "NUCLEAR REPULSION ENERGY"},
    },

    # Properties
    "properties": {
        "mulliken_charges": {"variables": "MULLIKEN_CHARGES", "skip_null": True},
        "lowdin_charges": {"variables": "LOWDIN_CHARGES", "skip_null": True},
        "wiberg_lowdin_indices": {"variables": "WIBERG_LOWDIN_INDICES", "skip_null": True},
        "mayer_indices": {"variables": "MAYER_INDICES", "skip_null": True},
    },

    # SCF variables
    "scf": {
        "scf_one_electron_energy": {"variables": "ONE-ELECTRON ENERGY"},
        "scf_two_electron_energy": {"variables": "TWO-ELECTRON ENERGY"},
        "scf_dipole_moment": {"variables": ["SCF DIPOLE X", "SCF DIPOLE Y", "SCF DIPOLE Z"]},
        "scf_iterations": {"variables": "SCF ITERATIONS", "cast": int},
        "scf_total_energy": {"variables": "SCF TOTAL ENERGY"},
        "scf_vv10_energy": {"variables": "DFT VV10 ENERGY", "skip_zero": True},
        "scf_xc_energy": {"variables": "DFT XC ENERGY", "skip_zero": True},
        "scf_dispersion_correction_energy": {"variables": "DISPERSION CORRECTION ENERGY", "skip_zero": True},

        # SCF Properties (experimental)
        # "scf_quadrupole_moment": {"variables": ["SCF QUADRUPOLE " + x for x in ["XX", "XY", "XZ", "YY", "YZ", "ZZ"]], "skip_null": True},
        # "scf_mulliken_charges": {"variables": "MULLIKEN_CHARGES", "skip_null": True},
        # "scf_lowdin_charges": {"variables": "LOWDIN_CHARGES", "skip_null": True},
        # "scf_wiberg_lowdin_indices": {"variables": "WIBERG_LOWDIN_INDICES", "skip_null": True},
        # "scf_mayer_indices": {"variables": "MAYER_INDICES", "skip_null": True},
    },

    # MP2 variables
    "mp2": {
        "mp2_same_spin_correlation_energy": {"variables": "MP2 SAME-SPIN CORRELATION ENERGY"},
        "mp2_opposite_spin_correlation_energy": {"variables": "MP2 OPPOSITE-SPIN CORRELATION ENERGY"},
        "mp2_singles_energy": {"variables": "NYI", "default": 0.0},
        "mp2_doubles_energy": {"variables": "MP2 CORRELATION ENERGY"},
        "mp2_total_correlation_energy": {"variables": "MP2 CORRELATION ENERGY"},
        "mp2_total_energy": {"variables": "MP2 TOTAL ENERGY"},
    },

#    "": {"variables": },

} # yapf: disable

def _json_translation(value):
    """
    Translates from Psi4 to JSON data types
    """

    if isinstance(value, (psi4.core.Matrix, psi4.core.Vector)):
        value = value.np.ravel().tolist()
    elif isinstance(value, (psi4.core.Dimension)):
        value = value.to_tuple()
    elif isinstance(value, np.ndarray):
        value = value.ravel().tolist()

    return value

def _convert_variables(data, context=None):
    """
    Converts dictionaries of variables based on translation metadata
    """

    # Build the correct translation units
    if context is None:
        needed_vars = {}
        for v in _qcschema_translation.values():
            needed_vars.update(v)
    else:
        needed_vars = _qcschema_translation[context]

    ret = {}
    for key, var in needed_vars.items():

        # Get the actual variables
        if isinstance(var["variables"], str):
            value = data.get(var["variables"], None)
        elif isinstance(var["variables"], (list, tuple)):
            value = [data.get(x, None) for x in var["variables"]]
            if not any(value):
                value = None
        else:
            raise TypeError("variables type not understood.")

        # Handle skips
        if var.get("skip_zero", False) and (value == 0):
            continue

        if (var.get("skip_zero") or var.get("skip_null", False)) and (value is None):
            continue

        # Add defaults
        if (value is None) and ("default" in var):
            value = var["default"]

        # Cast if called
        if "cast" in var:
            value = var["cast"](value)

        ret[key] = _json_translation(value)

    return ret


## Execution functions


def _clean_psi_environ(do_clean):
    if do_clean:
        psi4.core.clean_variables()
        psi4.core.clean_options()
        psi4.core.clean()


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
        # qc_schema should be copied
        json_data = run_json_qc_schema(copy.deepcopy(json_data), clean)

    except Exception as error:
        json_data["error"] = repr(error)
        json_data["success"] = False

    if return_output:
        with open(outfile, 'r') as f:
            json_data["raw_output"] = f.read()
        atexit.register(os.unlink, outfile)

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

    # Pull out a standard set of Psi variables
    psi_props = psi4.core.scalar_variables()
    json_data["psi4:qcvars"] = psi_props

    # Still a bit of a mess at the moment add in local vars as well.
    for k, v in wfn.variables().items():
        if k not in json_data["psi4:qcvars"]:
            json_data["psi4:qcvars"][k] = _json_translation(v)

    # Handle the return result
    if json_data["driver"] == "energy":
        json_data["return_result"] = val
    elif json_data["driver"] in ["gradient", "hessian"]:
        json_data["return_result"] = val.np.ravel().tolist()
    elif json_data["driver"] == "properties":
        ret = {}
        mtd = json_data["model"]["method"].upper()

        # Dipole/quadrupole still special case
        if "dipole" in kwargs["properties"]:
            ret["dipole"] = [psi_props[mtd + " DIPOLE " + x] for x in ["X", "Y", "Z"]]
        if "quadrupole" in kwargs["properties"]:
            ret["quadrupole"] = [psi_props[mtd + " QUADRUPOLE " + x] for x in ["XX", "XY", "XZ", "YY", "YZ", "ZZ"]]
        ret.update(_convert_variables(wfn.variables(), context="properties"))

        json_data["return_result"] = ret
    else:
        raise KeyError("Did not understand Driver key %s." % json_data["driver"])

    props = {
        "calcinfo_nbasis": wfn.nso(),
        "calcinfo_nmo": wfn.nmo(),
        "calcinfo_nalpha": wfn.nalpha(),
        "calcinfo_nbeta": wfn.nbeta(),
        "calcinfo_natom": mol.geometry().shape[0],
    }
    props.update(_convert_variables(psi_props, context="generics"))
    props.update(_convert_variables(psi_props, context="scf"))

    # Write out MP2 keywords
    if "MP2 CORRELATION ENERGY" in psi_props:
        props.update(_convert_variables(psi_props, context="mp2"))

    json_data["properties"] = props
    json_data["success"] = True

    # Reset state
    _clean_psi_environ(clean)

    json_data["schema_name"] = "qc_schema_output"

    return json_data
