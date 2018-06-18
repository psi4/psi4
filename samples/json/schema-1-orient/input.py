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



json_ret = psi4.json_wrapper.run_json(noorient_data)

# Orients to Z axis

