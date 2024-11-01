#! test QCSchema for gradient

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
      -0.1294769411935893,
      0.0,
      -1.494187339479985,
      1.0274465079245698,
      0.0,
      1.494187339479985,
      1.0274465079245698
    ],
    "symbols": [
      "O",
      "H",
      "H"
    ]
  },
  "driver": "gradient",
  "model": {
    "method": "HF",
    "basis": "cc-pVDZ"
  },
  "keywords": {"scf_type": "df"}
}

# Write expected output
expected_return_result = [
  0.0,
  0.0,
  -0.05959774096119619,
  0.0,
  -0.043039786289375104,
  0.02979887048056895,
  0.0,
  0.043039786289375104,
  0.02979887048056895
]
expected_properties = {
  "calcinfo_nbasis": 24,
  "calcinfo_nmo": 24,
  "calcinfo_nalpha": 5,
  "calcinfo_nbeta": 5,
  "calcinfo_natom": 3,
  "scf_one_electron_energy": -122.4452968291507,
  "scf_two_electron_energy": 37.62243738251799,
  "nuclear_repulsion_energy": 8.80146206062943,
  "scf_total_energy": -76.02139738600329,
  "return_energy": -76.02139738600329
}

json_ret = psi4.schema_wrapper.run_qcschema(json_data)

with open("output.json", "w") as ofile:
    json.dump(json_ret.json(), ofile, indent=2)


