import numpy as np
import psi4
import json

#Generate JSON data
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
    "driver": "tamps",
    "model": {
        "method": "MP2",
        "basis" "6-31g"
    },
    "keywords": {
        "scf_type": "df",
        "mp2_type": "df"
    }
}
json_ret = psi4.json_wrapper.run_json(json_data)

with open("output.json", "w") as ofile:
    json.dump(json_ret, ofile, indent=2)

