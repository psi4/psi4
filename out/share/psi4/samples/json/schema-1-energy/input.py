#! test QC_JSON Schema for energy

import numpy as np
import psi4
import json

# Generate JSON data
json_data = {
  "schema_name": "qcschema_input",
  "schema_version": 1,
  "return_output": True,
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
    "method": "CCSD(T)",
    "basis": "6-31g"
  },
  "keywords": {"scf_type": "df",
               "mp2_type": "df",
               "cc_type": "df",
               "scf_properties": ["mayer_indices"]}
}

# Write expected output
expected_return_result = -76.12069540875206
expected_properties = {
  "calcinfo_nbasis": 13,
  "calcinfo_nmo": 13,
  "calcinfo_nalpha": 5,
  "calcinfo_nbeta": 5,
  "calcinfo_natom": 3,
  "return_energy": expected_return_result,
  "scf_total_energy": -76.02139738600329,
  "nuclear_repulsion_energy": 8.801462060629428,
  "scf_one_electron_energy": -122.27509488480224,
  "scf_two_electron_energy": 37.49349095184536,
  "scf_dipole_moment": [
    0.0,
    0.0,
    1.040372174058
  ],
  "scf_iterations": 10,
  "scf_total_energy": -75.98014187232745,
  "mp2_same_spin_correlation_energy": -0.031030063236104254,
  "mp2_opposite_spin_correlation_energy": -0.10168342161187537,
  "mp2_singles_energy": 0.0,
  "mp2_doubles_energy": -0.13271348484797962,
  "mp2_correlation_energy": -0.13271348484797962,
  "mp2_total_energy": -76.11285535717543,

  "ccsd_same_spin_correlation_energy": -0.024524221169926093,
  "ccsd_opposite_spin_correlation_energy": -0.11488521989019397,
  "ccsd_singles_energy": 0.0,
  "ccsd_doubles_energy": -0.13940944106012007,
  "ccsd_correlation_energy": -0.13940944106012007,
  "ccsd_total_energy": -76.11955131338757,
  "ccsd_iterations": 12,
  "ccsd_prt_pr_correlation_energy": -0.140553536424603,
  "ccsd_prt_pr_total_energy": expected_return_result,
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
  "mp2_correlation_energy": -0.2069490608880642,
  "mp2_total_energy": -76.2283674281634,
  "return_energy": expected_return_result
}

# Switch run to exact MP2
json_data["model"] = {
    "method": "MP2",
    "basis": "cc-pVDZ"
  }
json_data["keywords"]["scf_type"] = "pk"
json_data["keywords"]["mp2_type"] = "conv"
json_data["return_output"] = True
json_data["extras"] = {"current_qcvars_only": True}
json_ret = psi4.json_wrapper.run_json(json_data)


# print(json.dumps(json_ret, indent=2))


assert "Ugur Bozkaya" in json_ret["raw_output"]

