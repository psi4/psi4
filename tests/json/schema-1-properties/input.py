#! test QC_JSON Schema for properties

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
    "no_com": True,
    "no_reorient": True
  },
  "driver": "properties",
  "model": {
    "method": "HF",
    "basis": "6-31G",
    "properties": [
      "dipole",
      "quadrupole",
      "mulliken_charges",
      "lowdin_charges",
      "wiberg_lowdin_indices",
      "mayer_indices",
    ]
  },
  "keywords": {"scf_type": "df",
               "mp2_type": "df",
               "e_convergence": 9}
}

# Write expected output
expected_return_result = {
  "dipole": [
    0.0,
    0.0,
    2.6443634497158492
  ],
  "quadrupole": [
    -7.300687696691922,
    0.0,
    0.0,
    -4.136264661490291,
    0.0,
    -5.872491231624151
  ],
  "mulliken_charges": [
    -0.7967275695997689,
    0.3983637847998823,
    0.3983637847998822
  ],
  "lowdin_charges": [
    -0.5945105406840803,
    0.29725527034203636,
    0.29725527034203636
  ],
  "wiberg_lowdin_indices": [
    0.0,
    0.9237385044125344,
    0.9237385044125329,
    0.9237385044125344,
    0.0,
    0.006992650019441531,
    0.9237385044125329,
    0.006992650019441531,
    0.0
  ],
  "mayer_indices": [
    0.0,
    0.802064044935596,
    0.8020640449355959,
    0.802064044935596,
    0.0,
    0.003020025778524219,
    0.8020640449355959,
    0.003020025778524219,
    0.0
  ]
}

expected_properties = {
  "calcinfo_nbasis": 13,
  "calcinfo_nmo": 13,
  "calcinfo_nalpha": 5,
  "calcinfo_nbeta": 5,
  "calcinfo_natom": 3,
  "scf_one_electron_energy": -122.27509111304202,
  "scf_two_electron_energy": 37.49348718008625,
  "nuclear_repulsion_energy": 8.80146206062943,
  "scf_total_energy": -75.98014187232634,
  "return_energy": -75.98014187232634
}


json_ret = psi4.json_wrapper.run_json(json_data)

with open("output.json", "w") as ofile:                                                     #TEST
    json.dump(json_ret, ofile, indent=2)                                                   #TEST

psi4.compare_integers(True, json_ret["success"], "JSON Success")                           #TEST
psi4.compare_strings("qc_schema_output", json_ret["schema_name"], "Schema Name")           #TEST
for k in expected_return_result.keys():                                                     #TEST
    psi4.compare_arrays(expected_return_result[k], json_ret["return_result"][k], 5, "Result: " + k.upper())  #TEST

for k in expected_properties.keys():                                                       #TEST
    psi4.compare_values(expected_properties[k], json_ret["properties"][k], 5, k.upper())   #TEST


