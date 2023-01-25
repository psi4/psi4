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
"""
Runs a JSON input psi file.
"""

__all__ = [
    "run_json",
    "run_qcschema",
]

import atexit
import copy
import datetime
import json
import os
import pprint
import sys
import traceback
import uuid
import warnings
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, Union

import numpy as np
import qcelemental as qcel
import qcengine as qcng
from psi4.extras import exit_printing
from psi4.header import print_header
from psi4.metadata import __version__
from psi4.driver import driver, p4util

from psi4 import core

pp = pprint.PrettyPrinter(width=120, compact=True, indent=1)




## Methods and properties blocks

methods_dict_ = {
    'energy': driver.energy,
    'gradient': driver.gradient,
    'properties': driver.properties,
    'hessian': driver.hessian,
    'frequency': driver.frequency,
}
default_properties_ = {
    "dipole", "quadrupole", "mulliken_charges", "lowdin_charges", "wiberg_lowdin_indices", "mayer_indices"
}

## QCSchema translation blocks

_qcschema_translation = {

    # Generics
    "generics": {
        "return_energy": {"variables": "CURRENT ENERGY"},
        "return_gradient": {"variables": "CURRENT GRADIENT"},
        "return_hessian": {"variables": "CURRENT HESSIAN"},
        # "nuclear_repulsion_energy": {"variables": "NUCLEAR REPULSION ENERGY"},  # use mol instead
    },

    # Properties
    "properties": {
        "mulliken_charges": {"variables": "MULLIKEN CHARGES", "skip_null": True},
        "lowdin_charges": {"variables": "LOWDIN CHARGES", "skip_null": True},
        "wiberg_lowdin_indices": {"variables": "WIBERG LOWDIN INDICES", "skip_null": True},
        "mayer_indices": {"variables": "MAYER INDICES", "skip_null": True},
    },

    # SCF variables
    "scf": {
        "scf_one_electron_energy": {"variables": "ONE-ELECTRON ENERGY"},
        "scf_two_electron_energy": {"variables": "TWO-ELECTRON ENERGY"},
        "scf_dipole_moment": {"variables": "SCF DIPOLE"},
        "scf_iterations": {"variables": "SCF ITERATIONS", "cast": int},
        "scf_total_energy": {"variables": "SCF TOTAL ENERGY"},
        "scf_vv10_energy": {"variables": "DFT VV10 ENERGY", "skip_zero": True},
        "scf_xc_energy": {"variables": "DFT XC ENERGY", "skip_zero": True},
        "scf_dispersion_correction_energy": {"variables": "DISPERSION CORRECTION ENERGY", "skip_zero": True},
        "scf_total_gradient": {"variables": "SCF TOTAL GRADIENT"},
        "scf_total_hessian": {"variables": "SCF TOTAL HESSIAN"},

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
        "mp2_singles_energy": {"variables": "MP2 SINGLES ENERGY"},
        "mp2_doubles_energy": {"variables": "MP2 DOUBLES ENERGY"},
        "mp2_correlation_energy": {"variables": "MP2 CORRELATION ENERGY"},
        "mp2_total_energy": {"variables": "MP2 TOTAL ENERGY"},
        "mp2_dipole_moment": {"variables": "NYI"}
    },

    "ccsd": {
        "ccsd_same_spin_correlation_energy": {"variables": "CCSD SAME-SPIN CORRELATION ENERGY"},
        "ccsd_opposite_spin_correlation_energy": {"variables": "CCSD OPPOSITE-SPIN CORRELATION ENERGY"},
        "ccsd_singles_energy": {"variables": "CCSD SINGLES ENERGY"},
        "ccsd_doubles_energy": {"variables": "CCSD DOUBLES ENERGY"},
        "ccsd_correlation_energy": {"variables": "CCSD CORRELATION ENERGY"},
        "ccsd_total_energy": {"variables": "CCSD TOTAL ENERGY"},
        "ccsd_dipole_moment": {"variables": "NYI"},
        "ccsd_iterations": {"variables": "CCSD ITERATIONS", "cast": int},
    },

    "ccsd(t)": {
        "ccsd_prt_pr_correlation_energy": {"variables": "CCSD(T) CORRELATION ENERGY"},
        "ccsd_prt_pr_total_energy": {"variables": "CCSD(T) TOTAL ENERGY"},
        "ccsd_prt_pr_dipole_moment": {"variables": "NYI"},
    }

#        "ccsd_singles_energy": {"variables": "NYI", "default": 0.0},
#    "": {"variables": },

} # yapf: disable


def _serial_translation(value, json=False):
    """
    Translates from Psi4 to JSON data types
    """
    if isinstance(value, (core.Dimension)):
        return value.to_tuple()

    if json:
        if isinstance(value, (core.Matrix, core.Vector)):
            return value.np.ravel().tolist()
        elif isinstance(value, np.ndarray):
            return value.ravel().tolist()
    else:
        if isinstance(value, (core.Matrix, core.Vector)):
            return value.np

    return value


def _convert_variables(data, context=None, json=False):
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

        conversion_factor = var.get("conversion_factor", None)

        # Get the actual variables
        if isinstance(var["variables"], str):
            value = data.get(var["variables"], None)
            if value is not None and conversion_factor:
                if isinstance(value, (int, float, np.ndarray)):
                    value *= conversion_factor
                elif isinstance(value, (core.Matrix, core.Vector)):
                    value.scale(conversion_factor)
        elif isinstance(var["variables"], (list, tuple)):
            value = [data.get(x, None) for x in var["variables"]]
            if not any(value):
                value = None
            elif conversion_factor:
                value = [x * conversion_factor for x in value]
        else:
            raise TypeError(f"Type of variable not understood: {key} {var}")

        # Handle skips
        if var.get("skip_zero", False) and (value == 0):
            continue

        if (var.get("skip_zero") or var.get("skip_null", False)) and (value is None):
            continue

        # Add defaults
        if (value is None) and ("default" in var):
            value = var["default"]

        # Cast if called
        if (value is not None) and ("cast" in var):
            value = var["cast"](value)

        ret[key] = _serial_translation(value, json=json)

    return ret


def _convert_basis(basis):
    """Converts a Psi4 basis object to a QCElemental basis.
    """
    centers = []
    symbols = []

    # Loop over centers
    for c in range(basis.molecule().natom()):
        center_shells = []
        symbols.append(basis.molecule().symbol(c).title())

        # Loop over shells *on* a center
        for s in range(basis.nshell_on_center(c)):
            shell = basis.shell(basis.shell_on_center(c, s))
            if shell.is_pure():
                htype = "spherical"
            else:
                htype = "cartesian"

            # Build the shell
            coefs = [[shell.coef(x) for x in range(shell.nprimitive)]]
            exps = [shell.exp(x) for x in range(shell.nprimitive)]
            qshell = qcel.models.basis.ElectronShell(angular_momentum=[shell.am],
                                                     harmonic_type=htype,
                                                     exponents=exps,
                                                     coefficients=coefs)
            center_shells.append(qshell)

        centers.append(qcel.models.basis.BasisCenter(electron_shells=center_shells))

    # Take unique to prevent duplicate data, doesn't matter too much
    hashes = [hash(json.dumps(centers[x].dict(), sort_keys=True)) for x in range(len(centers))]

    uniques = {k: v for k, v in zip(hashes, centers)}
    name_map = {}
    counter = defaultdict(int)

    # Generate reasonable names
    for symbol, h in zip(symbols, hashes):
        if h in name_map:
            continue

        counter[symbol] += 1

        name_map[h] = f"{basis.name()}_{symbol}{counter[symbol]}"

    center_data = {name_map[k]: v for k, v in uniques.items()}
    atom_map = [name_map[x] for x in hashes]

    ret = qcel.models.BasisSet(name=basis.name(), center_data=center_data, atom_map=atom_map)
    return ret


def _convert_wavefunction(wfn, context=None):

    basis = _convert_basis(wfn.basisset())
    # We expect CCA ordering.
    # Psi4 Cartesian is CCA (nothing to do)
    # Psi4 Spherical is in "Gaussian" reorder
    # Spherical Map: 0 -1 +1, ... -> -1, ..., 0, ..., +1

    spherical_maps = {}
    for L in range(wfn.basisset().max_am() + 1):
        mapper = list(range(L * 2 - 1, 0, -2)) + [0] + list(range(2, L * 2 + 1, 2))
        spherical_maps[L] = np.array(mapper)

    # Build a flat index that we can transform the AO quantities
    reorder = True
    ao_map = []
    cnt = 0
    for atom in basis.atom_map:
        center = basis.center_data[atom]
        for shell in center.electron_shells:
            if shell.harmonic_type == "cartesian":
                ao_map.append(np.arange(cnt, cnt + shell.nfunctions()))
            else:
                smap = spherical_maps[shell.angular_momentum[0]]
                ao_map.append(smap + cnt)
                reorder = True

            cnt += shell.nfunctions()

    ao_map = np.hstack(ao_map)

    # Build remap functions
    def re2d(mat, both=True):
        arr = np.array(mat)
        if reorder:
            if both:
                arr = arr[ao_map[:, None], ao_map]
            else:
                arr = arr[ao_map[:, None]]
        return arr

    # get occupations in orbital-energy ordering
    def sort_occs(noccpi, epsilon):
        occs = []
        for irrep, nocc in enumerate(noccpi):
            for i, e in enumerate(epsilon[irrep]):
                occs.append((e, int(i < nocc)))

        occs.sort(key = lambda x : x[0])
        return np.array([occ[1] for occ in occs])

    # Map back out what we can
    ret = {
        "basis": basis,
        "restricted": (wfn.same_a_b_orbs() and wfn.same_a_b_dens()),

        # Return results
        "orbitals_a": "scf_orbitals_a",
        "orbitals_b": "scf_orbitals_b",
        "density_a": "scf_density_a",
        "density_b": "scf_density_b",
        "fock_a": "scf_fock_a",
        "fock_b": "scf_fock_b",
        "eigenvalues_a": "scf_eigenvalues_a",
        "eigenvalues_b": "scf_eigenvalues_b",
        "occupations_a": "scf_occupations_a",
        "occupations_b": "scf_occupations_b",
        # "h_effective_a": re2d(wfn.H()),
        # "h_effective_b": re2d(wfn.H()),

        # SCF quantities
        "scf_orbitals_a": re2d(wfn.Ca_subset("AO", "ALL"), both=False),
        "scf_orbitals_b": re2d(wfn.Cb_subset("AO", "ALL"), both=False),
        "scf_density_a": re2d(wfn.Da_subset("AO")),
        "scf_density_b": re2d(wfn.Db_subset("AO")),
        "scf_fock_a": re2d(wfn.Fa_subset("AO")),
        "scf_fock_b": re2d(wfn.Fa_subset("AO")),
        "scf_eigenvalues_a": wfn.epsilon_a_subset("AO", "ALL"),
        "scf_eigenvalues_b": wfn.epsilon_b_subset("AO", "ALL"),
        "scf_occupations_a": sort_occs((wfn.doccpi() + wfn.soccpi()).to_tuple(), wfn.epsilon_a().nph),
        "scf_occupations_b": sort_occs(wfn.doccpi().to_tuple(), wfn.epsilon_b().nph),
    }

    return ret


## Execution functions


def _clean_psi_environ(do_clean: bool):
    """Reset work environment to new Psi4 instance state.
    This includes global variables (P::e.globals, P::e.arrays, P::e.options) and any
    non-explicitly-retained PSIO-managed scratch files.

    """
    if do_clean:
        core.clean_variables()
        core.clean_options()
        core.clean()


def _clean_psi_output(do_clean: bool, outfile: str):
    """Reset primary output file depending on run mode.
    When ``True``, remove the output file; when ``False``, close the output file.
    The latter is appropriate when calling :py:func:`psi4.run_qcschema` *from a Psi4 Python session*
    as otherwise, the parent session's own files could get cleaned away.

    """
    if do_clean:
        atexit.register(_quiet_remove, outfile)
    else:
        core.close_outfile()


def _read_output(outfile):
    try:
        with open(outfile, 'r') as f:
            output = f.read()

        return output
    except OSError:
        return "Could not read output file."


def _quiet_remove(filename):
    """
    Destroy the created at file at exit, pass if error.
    """
    try:
        os.unlink(filename)
    except OSError:
        pass


def run_qcschema(
    input_data: Union[Dict[str, Any], qcel.models.AtomicInput],
    clean: bool = True,
    postclean: bool = True,
) -> Union[qcel.models.AtomicResult, qcel.models.FailedOperation]:
    """Run a quantum chemistry job specified by :py:class:`qcelemental.models.AtomicInput` **input_data** in |PSIfour|.

    Parameters
    ----------
    input_data
        Quantum chemistry job in either AtomicInput class or dictionary form.
    clean
        Reset global QCVariables, options, and scratch files to default state.
    postclean
        When ``False``, *remove* the output file since absorbed into AtomicResult.
        When ``True``, simply *close* the output file. True is useful when calling
        from a Psi4 session to avoid removing the parent Psi4's output file.

    Returns
    -------
    qcelemental.models.AtomicResult
        Full record of quantum chemistry calculation, including output text. Returned upon job success.
    qcelemental.models.FailedOperation
        Record to diagnose calculation failure, including output text and input specification. Returned upon job failure.

    """
    outfile = os.path.join(core.IOManager.shared_object().get_default_path(), str(uuid.uuid4()) + ".qcschema_tmpout")
    core.set_output_file(outfile, False)
    print_header()

    start_time = datetime.datetime.now()

    try:
        input_model = qcng.util.model_wrapper(input_data, qcel.models.AtomicInput)

        # Echo the infile on the outfile
        core.print_out("\n  ==> Input QCSchema <==\n")
        core.print_out("\n--------------------------------------------------------------------------\n")
        core.print_out(pp.pformat(json.loads(input_model.json())))
        core.print_out("\n--------------------------------------------------------------------------\n")

        keep_wfn = input_model.protocols.wavefunction != 'none'

        # qcschema should be copied
        ret_data = run_json_qcschema(input_model.dict(), clean, False, keep_wfn=keep_wfn)
        ret_data["provenance"].update({
            "creator": "Psi4",
            "version": __version__,
            "routine": "psi4.schema_runner.run_qcschema"
        })
        ret_data["native_files"]["input"] = json.dumps(json.loads(input_model.json()), indent=1)

        exit_printing(start_time=start_time, success=True)

        ret = qcel.models.AtomicResult(**ret_data, stdout=_read_output(outfile))

    except Exception as exc:

        if not isinstance(input_data, dict):
            input_data = input_data.dict()

        input_data = input_data.copy()
        input_data["stdout"] = _read_output(outfile)
        ret = qcel.models.FailedOperation(input_data=input_data,
                                          success=False,
                                          error={
                                              'error_type': type(exc).__name__,
                                              'error_message': input_data["stdout"] + ''.join(traceback.format_exception(*sys.exc_info())),
                                          })

    _clean_psi_output(postclean, outfile)

    return ret


def run_json(json_data: Dict[str, Any], clean: bool = True) -> Dict[str, Any]:

    warnings.warn(
        "Using `psi4.schema_wrapper.run_json` or `psi4.json_wrapper.run_json` instead of `psi4.schema_wrapper.run_qcschema` is deprecated, and as soon as 1.5 it will stop working\n",
        category=FutureWarning)

    # Set scratch
    if "scratch_location" in json_data:
        psi4_io = core.IOManager.shared_object()
        psi4_io.set_default_path(json_data["scratch_location"])

    # Direct output
    outfile = os.path.join(core.IOManager.shared_object().get_default_path(), str(uuid.uuid4()) + ".json_out")
    core.set_output_file(outfile, False)

    # Set memory
    if "memory" in json_data:
        p4util.set_memory(json_data["memory"])

    # Do we return the output?
    return_output = json_data.pop("return_output", False)
    if return_output:
        json_data["raw_output"] = "Not yet run."

    # Set a few flags
    json_data["raw_output"] = None
    json_data["success"] = False
    json_data["provenance"] = {"creator": "Psi4", "version": __version__, "routine": "psi4.json_wrapper.run_json"}

    # Attempt to run the computer
    try:
        # qcschema should be copied
        json_data = run_json_qcschema(copy.deepcopy(json_data), clean, True)

    except Exception as exc:
        json_data["error"] = {
            'error_type': type(exc).__name__,
            'error_message': ''.join(traceback.format_exception(*sys.exc_info())),
        }
        json_data["success"] = False

        json_data["raw_output"] = _read_output(outfile)

    if return_output:
        json_data["raw_output"] = _read_output(outfile)

    atexit.register(_quiet_remove, outfile)

    return json_data


def run_json_qcschema(json_data, clean, json_serialization, keep_wfn=False):
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
    if json_data["schema_name"] in ["qc_schema_input", "qcschema_input"]:
        json_data["schema_name"] = "qcschema_input"
    else:
        raise KeyError("Schema name of '{}' not understood".format(json_data["schema_name"]))

    if json_data["schema_version"] != 1:
        raise KeyError("Schema version of '{}' not understood".format(json_data["schema_version"]))

    if json_data.get("nthreads", False) is not False:
        core.set_num_threads(json_data["nthreads"], quiet=True)

    # Build molecule
    if "schema_name" in json_data["molecule"]:
        molschemus = json_data["molecule"]  # dtype >=2
    else:
        molschemus = json_data  # dtype =1
    mol = core.Molecule.from_schema(molschemus, nonphysical=True)

    # Update molecule geometry as we orient and fix_com
    json_data["molecule"]["geometry"] = mol.geometry().np.ravel().tolist()

    # Set options
    ## The Forte plugin needs special treatment.
    try:
        import forte # Needed for Forte options to run.
    except ImportError:
        pass
    else:
        # Initialization tasks with Psi options.
        psi_options = core.get_options()
        current_module = psi_options.get_current_module()
        # Get the current Forte options from Forte
        forte_options = forte.ForteOptions()
        forte.register_forte_options(forte_options)
        psi_options.set_current_module('FORTE')
        forte_options.push_options_to_psi4(psi_options)
        # Restore current module
        psi_options.set_current_module(current_module)

    kwargs = json_data["keywords"].pop("function_kwargs", {})
    p4util.set_options(json_data["keywords"])

    # Setup the computation
    method = json_data["model"]["method"]
    core.set_global_option("BASIS", json_data["model"]["basis"])
    kwargs.update({"return_wfn": True, "molecule": mol})

    # Handle special properties case
    if json_data["driver"] == "properties":
        if "properties" not in kwargs:
            kwargs["properties"] = list(default_properties_)

    # Actual driver run
    val, wfn = methods_dict_[json_data["driver"]](method, **kwargs)

    # Pull out a standard set of SCF properties
    if "extras" not in json_data:
        json_data["extras"] = {}
    json_data["extras"]["qcvars"] = {}

    current_qcvars_only = json_data["extras"].get("current_qcvars_only", False)
    if json_data["extras"].get("wfn_qcvars_only", False):
        psi_props = wfn.variables(include_deprecated_keys=(not current_qcvars_only))
    else:
        psi_props = core.variables(include_deprecated_keys=(not current_qcvars_only))
        for k, v in psi_props.items():
            if k not in json_data["extras"]["qcvars"]:
                json_data["extras"]["qcvars"][k] = _serial_translation(v, json=json_serialization)

    # Still a bit of a mess at the moment add in local vars as well.
    for k, v in wfn.variables().items():
        if k not in json_data["extras"]["qcvars"]:
            # interpreting wfn_qcvars_only as no deprecated qcvars either
            if not (json_data["extras"].get("wfn_qcvars_only", False) and (
                any([k.upper().endswith(" DIPOLE " + cart) for cart in ["X", "Y", "Z"]])
                or any([k.upper().endswith(" QUADRUPOLE " + cart) for cart in ["XX", "YY", "ZZ", "XY", "XZ", "YZ"]])
                or k.upper()
                in [
                    "SOS-MP2 CORRELATION ENERGY",
                    "SOS-MP2 TOTAL ENERGY",
                    "SOS-PI-MP2 CORRELATION ENERGY",
                    "SOS-PI-MP2 TOTAL ENERGY",
                    "SCS-MP3 CORRELATION ENERGY",
                    "SCS-MP3 TOTAL ENERGY",
                ]
            )):
                json_data["extras"]["qcvars"][k] = _serial_translation(v, json=json_serialization)

    # Add in handling of matrix arguments which need to be obtained by a
    # a function call.

    if json_data["model"]["method"].lower() in ["ccsd"] :
        if json_data["extras"].get("psi4:save_tamps", False):
            if type(wfn.reference_wavefunction()) is core.RHF :
                json_data["extras"]["psi4:tamps"] = {}
                json_data["extras"]["psi4:tamps"]["tIjAb"] = wfn.get_amplitudes()["tIjAb"].to_array().tolist()
                json_data["extras"]["psi4:tamps"]["tIA"] = wfn.get_amplitudes()["tIA"].to_array().tolist()
                json_data["extras"]["psi4:tamps"]["Da"] = wfn.Da().to_array().tolist()

        
    # Handle the return result
    if json_data["driver"] == "energy":
        json_data["return_result"] = val
    elif json_data["driver"] in ["gradient", "hessian"]:
        json_data["return_result"] = _serial_translation(val, json=json_serialization)
    elif json_data["driver"] == "properties":
        ret = {}
        mtd = json_data["model"]["method"].upper()

        # Dipole/quadrupole still special case
        if "dipole" in kwargs["properties"]:
            ret["dipole"] = _serial_translation(psi_props[mtd + " DIPOLE"], json=json_serialization)
        if "quadrupole" in kwargs["properties"]:
            ret["quadrupole"] = _serial_translation(psi_props[mtd + " QUADRUPOLE"], json=json_serialization)
        ret.update(_convert_variables(wfn.variables(), context="properties", json=json_serialization))

        json_data["return_result"] = ret
    else:
        raise KeyError("Did not understand Driver key %s." % json_data["driver"])

    props = {
        "calcinfo_nbasis": wfn.nso(),
        "calcinfo_nmo": wfn.nmo(),
        "calcinfo_nalpha": wfn.nalpha(),
        "calcinfo_nbeta": wfn.nbeta(),
        "calcinfo_natom": mol.geometry().shape[0],
        "nuclear_repulsion_energy": mol.nuclear_repulsion_energy(),  # use this b/c psivar is monomer for SAPT
    }
    props.update(_convert_variables(psi_props, context="generics", json=json_serialization))
    if not list(set(['CBS NUMBER', 'NBODY NUMBER', 'FINDIF NUMBER']) & set(json_data["extras"]["qcvars"].keys())):
        props.update(_convert_variables(psi_props, context="scf", json=json_serialization))

    # Write out post-SCF keywords
    if "MP2 CORRELATION ENERGY" in psi_props:
        props.update(_convert_variables(psi_props, context="mp2", json=json_serialization))

    if "CCSD CORRELATION ENERGY" in psi_props:
        props.update(_convert_variables(psi_props, context="ccsd", json=json_serialization))

    if "CCSD(T) CORRELATION ENERGY" in psi_props:
        props.update(_convert_variables(psi_props, context="ccsd(t)", json=json_serialization))

    json_data["properties"] = props
    json_data["success"] = True

    json_data["provenance"]["module"] = wfn.module()

    if keep_wfn:
        json_data["wavefunction"] = _convert_wavefunction(wfn)

    files = {
        "psi4.grad": Path(core.get_writer_file_prefix(wfn.molecule().name()) + ".grad"),
        "psi4.hess": Path(core.get_writer_file_prefix(wfn.molecule().name()) + ".hess"),
        # binary "psi4.180.npy": Path(core.get_writer_file_prefix(wfn.molecule().name()) + ".180.npy"),
        "timer.dat": Path("timer.dat"),  # ok for `psi4 --qcschema` but no file collected for `qcengine.run_program(..., "psi4")`
    }
    json_data["native_files"] = {fl: flpath.read_text() for fl, flpath in files.items() if flpath.exists()}

    # Reset state
    _clean_psi_environ(clean)

    json_data["schema_name"] = "qcschema_output"

    return json_data
