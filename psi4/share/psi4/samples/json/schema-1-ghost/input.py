#! test QCSchema with ghost atoms

import numpy as np
import psi4
import json

# Generate JSON data
json_data = {
    "schema_name": "qc_schema_input",
    "schema_version": 1,
    "molecule": {
        "geometry": [
            0.0,
            0.0,
            -5.0,
            0.0,
            0.0,
            5.0,
        ],
        "symbols": ["He", "He"],
        "real": [True, False]
    },
    "driver": "energy",
    "model": {
        "method": "SCF",
        "basis": "cc-pVDZ"
    },
    "keywords": {
        "scf_type": "df"
    },
}
psi4.set_num_threads(1)
psi4.set_memory(1024 * 1024 * 1024)

# Write expected output
expected_return_result = -2.85518836280515
expected_properties = {
    'calcinfo_nbasis': 10,
    'calcinfo_nmo': 10,
    'calcinfo_nalpha': 1,
    'calcinfo_nbeta': 1,
    'calcinfo_natom': 2,
    'scf_one_electron_energy': -3.8820496359492576,
    'scf_two_electron_energy': 1.0268612731441076,
    'nuclear_repulsion_energy': 0.0,
    'scf_total_energy': -2.85518836280515,
    'return_energy': -2.85518836280515
}

json_ret = psi4.schema_wrapper.run_qcschema(json_data)




