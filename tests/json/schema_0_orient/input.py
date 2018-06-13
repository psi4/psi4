import numpy as np
import psi4
import json
import copy

# Generate JSON data
json_data = {
    "schema_name": "QC_JSON",
    "schema_version": 0,
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

psi4.json_wrapper.run_json(json_data)

# Orients to Z axis
psi4.compare_integers(True, json_data["success"], "JSON Success")  #TEST
psi4.compare_values(expected_return_result, json_data["return_result"], 5, "Return Value")  #TEST
psi4.compare_values(0.0, json_data["properties"]["scf_dipole_moment"][0], 5, "DIPOLE X")  #TEST
psi4.compare_values(0.0, json_data["properties"]["scf_dipole_moment"][1], 5, "DIPOLE Y")  #TEST
psi4.compare_values(linear_dipole, json_data["properties"]["scf_dipole_moment"][2], 5, "DIPOLE Z")  #TEST

dist = np.linalg.norm(np.array(json_data["molecule"]["geometry"])[:3] - np.array(np.array(json_data["molecule"]["geometry"])[3:])) #TEST
psi4.compare_values(1.732, dist, 4, "HF Bond Distance")  #TEST


psi4.json_wrapper.run_json(noorient_data)

# Orients to Z axis
psi4.compare_integers(True, noorient_data["success"], "JSON Success")  #TEST
psi4.compare_values(expected_return_result, noorient_data["return_result"], 5, "Return Value")  #TEST
psi4.compare_values(0.0, noorient_data["properties"]["scf_dipole_moment"][0], 5, "DIPOLE X")  #TEST
psi4.compare_values(linear_dipole, noorient_data["properties"]["scf_dipole_moment"][1], 5, "DIPOLE Y")  #TEST
psi4.compare_values(0.0, noorient_data["properties"]["scf_dipole_moment"][2], 5, "DIPOLE Z")  #TEST
psi4.compare_arrays([0.0, 0.0, 0.0], noorient_data["molecule"]["geometry"][:3], 5, "H Position")  #TEST
psi4.compare_arrays([0.0, 1.732, 0.0], noorient_data["molecule"]["geometry"][3:], 5, "F Position")  #TEST

