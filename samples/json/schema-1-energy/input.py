#! test QC_JSON Schema for energy

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
    ],
    "connectivity" : [
    (0, 1, 1.0),
    (0, 2, 1.0)
    ]
  },
  "driver": "energy",
  "model": {
    "method": "MP2",
    "basis": "cc-pVDZ"
  },
  "keywords": {"scf_type": "df",
               "mp2_type": "df",
               "scf_properties": ["mayer_indices"]}
}

# Write expected output
expected_return_result = -76.22831410222477
expected_properties = {
  "calcinfo_nbasis": 24,
  "calcinfo_nmo": 24,
  "calcinfo_nalpha": 5,
  "calcinfo_nbeta": 5,
  "calcinfo_natom": 3,
  "scf_one_electron_energy": -122.44529682915068,
  "scf_two_electron_energy": 37.622437382517965,
  "nuclear_repulsion_energy": 8.80146206062943,
  "scf_total_energy": -76.02139738600329,
  "mp2_same_spin_correlation_energy": -0.05202760538221721,
  "mp2_opposite_spin_correlation_energy": -0.1548891108392641,
  "mp2_singles_energy": 0.0,
  "mp2_doubles_energy": -0.20691671622148142,
  "mp2_total_correlation_energy": -0.20691671622148142,
  "mp2_total_energy": -76.22831410222477,
  "return_energy": expected_return_result
}

json_ret = psi4.json_wrapper.run_json(json_data)




# Expected output with exact MP2
expected_return_result = -76.2283674281634
expected_properties = {
  "calcinfo_nbasis": 24,
  "calcinfo_nmo": 24,
  "calcinfo_nalpha": 5,
  "calcinfo_nbeta": 5,
  "calcinfo_natom": 3,
  "scf_one_electron_energy": -122.44534537436829,
  "scf_two_electron_energy": 37.62246494646352,
  "nuclear_repulsion_energy": 8.80146206062943,
  "scf_total_energy": -76.02141836727533,
  "mp2_same_spin_correlation_energy": -0.051980792907589016,
  "mp2_opposite_spin_correlation_energy": -0.1549682679804691,
  "mp2_singles_energy": 0.0,
  "mp2_doubles_energy": -0.2069490608880642,
  "mp2_total_correlation_energy": -0.2069490608880642,
  "mp2_total_energy": -76.2283674281634,
  "return_energy": expected_return_result
}

# Switch run to exact MP2
json_data["keywords"]["scf_type"] = "pk"
json_data["keywords"]["mp2_type"] = "conv"
json_ret = psi4.json_wrapper.run_json(json_data)


# print(json.dumps(json_ret, indent=2))


