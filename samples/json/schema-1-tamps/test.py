#! test QCSchema for CCSD amplitudes saving

import numpy as np
import psi4
import json

#Generate JSON data

method = "CCSD"
basis = "sto-3g"
json_data = {
    "schema_name": "qcschema_input",
    "schema_version": 1,
    "molecule": {
        "geometry": [
            0.0,
            0.0,
            -0.129476941212,
            0.0,
            -1.494187339480,
            1.027446507906,
            0.0,
            1.494187339480,
            1.027446507906
        ],
        "symbols": [
          "O",
          "H",
          "H"
        ],
        "connectivity" : [
            (0, 1, 1.0),
            (0, 2, 1.0)
        ],
        "fix_symmetry": "c1"
    },
    "driver": "energy",
    "model": {
        "method": method,
        "basis": basis
    },
    "keywords": {
        "scf_type": "df",
        "mp2_type": "df"
    },
    "extras": {
        "psi4:save_tamps": True
    }
}

h2o_test = psi4.geometry("""
0 1
units bohr
O 0.0 0.0 -0.129476941212
H 0.0 -1.494187339480 1.027446507906
H 0.0 1.494187339480 1.027446507906
symmetry c1
""")

h2o_schema = psi4.core.Molecule.from_schema(json_data)

psi4.compare_values(h2o_test.geometry().to_array(),
                    h2o_schema.geometry().to_array())

json_ret = psi4.schema_wrapper.run_qcschema(json_data, True)
psi4.set_options({"basis": basis, "scf_type": "df", "mp2_type": "df"})
e, wfn = psi4.energy(method, return_wfn = True)

tIJAB = wfn.get_amplitudes()["tIjAb"].to_array()
tIA = wfn.get_amplitudes()["tIA"].to_array()
Da = wfn.Da().to_array()

with open("output.json", "w") as ofile:
    json.dump(json_ret.json(), ofile, indent=2)

# Make sure the amplitudes compted from h2o_test are actually assigned.
if all(map(lambda x: x == 0, tIJAB.flatten())) or \
   all(map(lambda x: x == 0, tIA.flatten())) or \
   all(map(lambda x: x == 0, Da.flatten())) :
    raise Exception("Error in computing values to compare against. Could not compare.")
else :
    psi4.compare_values(tIJAB, np.array(json_ret.extras["psi4:tamps"]["tIjAb"]))
    psi4.compare_values(tIA, np.array(json_ret.extras["psi4:tamps"]["tIA"]))
    psi4.compare_values(Da, np.array(json_ret.extras["psi4:tamps"]["Da"]))

