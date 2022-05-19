#! test QCSchema noncontiguous mol

import numpy as np
import psi4
import json

# Generate JSON data
json_data = {
    "schema_name": "qc_schema_input",
    "schema_version": 1,
    "molecule": {
        "schema_name": "qcschema_molecule",
        "schema_version": 2,
        "geometry": [
            0.0, 0.0, -0.1294769411935893,
            0.0, -1.494187339479985, 1.0274465079245698,
            0.0,  1.494187339479985, 1.0274465079245698
        ],
        "symbols": ["O", "H", "H"],
        "fragments": [[2, 0, 1]]
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

# Check non-contiguous fragment throws
json_ret = psi4.schema_wrapper.run_qcschema(json_data)

psi4.compare_integers(False, json_ret.success, "JSON Failure")                           #TEST
psi4.compare_integers("non-contiguous frag" in json_ret.error.error_message, True, "Contiguous Fragment Error")  #TEST

# Check symbol length errors
del json_data["molecule"]["fragments"]
json_data["molecule"]["symbols"] = ["O", "H"]
json_ret = psi4.schema_wrapper.run_qcschema(json_data)

psi4.compare_integers(False, json_ret.success, "JSON Failure")                           #TEST
psi4.compare_integers("dropped atoms!" in json_ret.error.error_message, True, "Symbol Error")  #TEST

# Check keyword errors
json_data["molecule"]["symbols"] = ["O", "H", "H"]
json_data["model"] = {"method": "SCF", "basis": "sto-3g"}
json_data["keywords"] = {"scf_type": "super_df"}
json_ret = psi4.schema_wrapper.run_qcschema(json_data)

psi4.compare_integers(False, json_ret.success, "JSON Failure")                           #TEST
psi4.compare_integers("valid choice" in json_ret.error.error_message, True, "Keyword Error")  #TEST
