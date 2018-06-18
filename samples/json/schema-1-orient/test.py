#! test QC_JSON Schema mol orientation

import numpy as np
import psi4
import json
import copy

# Generate JSON data
json_data = {
    "schema_name": "qc_schema_input",
    "schema_version": 1,
    "molecule": {
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
linear_dipole = 1.948993625469663 

json_ret = psi4.json_wrapper.run_json(json_data)

# Orients to Z axis
psi4.compare_integers(True, json_ret["success"], "JSON Success")  #TEST
psi4.compare_values(expected_return_result, json_ret["return_result"], 5, "Return Value")  #TEST
psi4.compare_values(0.0, json_ret["properties"]["scf_dipole_moment"][0], 3, "DIPOLE X")  #TEST
psi4.compare_values(0.0, json_ret["properties"]["scf_dipole_moment"][1], 3, "DIPOLE Y")  #TEST
psi4.compare_values(linear_dipole, json_ret["properties"]["scf_dipole_moment"][2], 3, "DIPOLE Z")  #TEST

dist = np.linalg.norm(np.array(json_ret["molecule"]["geometry"])[:3] - np.array(np.array(json_ret["molecule"]["geometry"])[3:])) #TEST
psi4.compare_values(1.732, dist, 4, "HF Bond Distance")  #TEST


json_ret = psi4.json_wrapper.run_json(noorient_data)

# Orients to Z axis
psi4.compare_integers(True, json_ret["success"], "JSON Success")  #TEST
psi4.compare_values(expected_return_result, json_ret["return_result"], 5, "Return Value")  #TEST
psi4.compare_values(0.0, json_ret["properties"]["scf_dipole_moment"][0], 3, "DIPOLE X")  #TEST
psi4.compare_values(linear_dipole, json_ret["properties"]["scf_dipole_moment"][1], 3, "DIPOLE Y")  #TEST
psi4.compare_values(0.0, json_ret["properties"]["scf_dipole_moment"][2], 3, "DIPOLE Z")  #TEST
psi4.compare_arrays([0.0, 0.0, 0.0], json_ret["molecule"]["geometry"][:3], 3, "H Position")  #TEST
psi4.compare_arrays([0.0, 1.732, 0.0], json_ret["molecule"]["geometry"][3:], 3, "F Position")  #TEST

