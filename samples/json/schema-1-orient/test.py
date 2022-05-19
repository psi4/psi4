#! test QCSchema mol orientation

import numpy as np
import psi4
import json
import copy

# Generate JSON data
json_data = {
    "schema_name": "qc_schema_input",
    "schema_version": 1,
    "molecule": {
        "schema_name": "qcschema_molecule",
        "schema_version": 2,
        "geometry": [
            0.0,
            0.0,
            0.0,
            0.0,
            1.732,
            0.0,
        ],
        "symbols": ["F", "H"]
    },
    "driver": "energy",
    "model": {
        "method": "SCF",
        "basis": "cc-pVDZ"
    },
    "keywords": {
        "scf_type": "df"
    }
}

noorient_data = copy.deepcopy(json_data)
noorient_data["molecule"]["fix_orientation"] = True
noorient_data["molecule"]["fix_com"] = True

# Write expected output
expected_return_result = -100.0194177509218
linear_dipole = 0.7667930938

json_ret = psi4.schema_wrapper.run_qcschema(json_data)

# Orients to Z axis
psi4.compare_integers(True, json_ret.success, "JSON Success")  #TEST
psi4.compare_values(expected_return_result, json_ret.return_result, 5, "Return Value")  #TEST
psi4.compare_values(0.0, json_ret.properties.scf_dipole_moment[0], 3, "DIPOLE X")  #TEST
psi4.compare_values(0.0, json_ret.properties.scf_dipole_moment[1], 3, "DIPOLE Y")  #TEST
psi4.compare_values(linear_dipole, json_ret.properties.scf_dipole_moment[2], 3, "DIPOLE Z")  #TEST

dist = np.linalg.norm(json_ret.molecule.geometry[0] - json_ret.molecule.geometry[1])  #TEST
psi4.compare_values(1.732, dist, 4, "HF Bond Distance")  #TEST


json_ret = psi4.schema_wrapper.run_qcschema(noorient_data)

# Orients to Z axis
psi4.compare_integers(True, json_ret.success, "JSON Success")  #TEST
psi4.compare_values(expected_return_result, json_ret.return_result, 5, "Return Value")  #TEST
psi4.compare_values(0.0, json_ret.properties.scf_dipole_moment[0], 3, "DIPOLE X")  #TEST
psi4.compare_values(linear_dipole, json_ret.properties.scf_dipole_moment[1], 3, "DIPOLE Y")  #TEST
psi4.compare_values(0.0, json_ret.properties.scf_dipole_moment[2], 3, "DIPOLE Z")  #TEST
psi4.compare_arrays([0.0, 0.0, 0.0], json_ret.molecule.geometry[0], 3, "H Position")  #TEST
psi4.compare_arrays([0.0, 1.732, 0.0], json_ret.molecule.geometry[1], 3, "F Position")  #TEST

