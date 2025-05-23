#! test QCSchema for properties

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
import copy
json_data_copy = copy.deepcopy(json_data)

# Write expected output (dipole & quadrupole in au)
expected_return_result = {
  "dipole": [
    0.0,
    0.0,
    1.04037263345
  ],
  "quadrupole": [
    -5.42788218076,
    0.0,
    0.0,
    0.0,
    -3.07521129293,
    0.0,
    0.0,
    0.0,
    -4.36605314966
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


## with deprecated `run_json`

json_ret = psi4.json_wrapper.run_json(json_data)





## with current `run_qcschema`

# first, modernize input and output models
json_data = json_data_copy
json_data["schema_name"] = "qcschema_input"
json_data["molecule"]["fix_com"] = True
json_data["molecule"]["fix_orientation"] = True
json_data["molecule"].pop("no_com")
json_data["molecule"].pop("no_reorient")
expected_return_result = copy.deepcopy(expected_return_result)
expected_return_result["dipole"] = np.array(expected_return_result["dipole"]).reshape((3,))
expected_return_result["quadrupole"] = np.array(expected_return_result["quadrupole"]).reshape((3, 3))
expected_return_result["wiberg_lowdin_indices"] = np.array(expected_return_result["wiberg_lowdin_indices"]).reshape((3, 3))
expected_return_result["mayer_indices"] = np.array(expected_return_result["mayer_indices"]).reshape((3, 3))

json_ret = psi4.schema_wrapper.run_qcschema(json_data).dict()

# can't write msgpack arrays to json


