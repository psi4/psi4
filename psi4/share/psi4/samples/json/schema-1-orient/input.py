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



json_ret = psi4.schema_wrapper.run_qcschema(noorient_data)

# Orients to Z axis

