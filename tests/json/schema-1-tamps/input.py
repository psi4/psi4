import numpy as np
import psi4
import json

#Generate JSON data

method = "CCSD"
basis = "sto-3g"
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
        "psi4:arrays": {
            "tIjAb": [],
            "tIA": [],
            "Da": []
        }
    }
}
json_ret = psi4.json_wrapper.run_json_qcschema(json_data, True, True)

h2o_test = psi4.geometry("""
0 1
O 0.0 0.0 0.1294769411935893
H 0.0 -1.494187339479985 1.0274465079245698
H 0.0 1.494187339479985 1.0274465079245698
symmetry c1
""")
e, wfn = psi4.energy(f"{method}/{basis}", return_wfn = True)
tIJAB = wfn.get_amplitudes()["tIjAb"].to_array()
tIA = wfn.get_amplitudes()["tIA"].to_array()
Da = wfn.Da().to_array()

with open("output.json", "w") as ofile:
    json.dump(json_ret, ofile, indent=2)

# Make sure the amplitudes compted from h2o_test are actually assigned.
if all(map(lambda x: x == 0, tIJAB.flatten())) or \
   all(map(lambda x: x == 0, tIA.flatten())) or \
   all(map(lambda x: x == 0, Da.flatten())) :
    print("Error in computing values to compare against. Could not compare.")
else :
    print("Comparing T2")
    psi4.compare_values(tIJAB, np.array(json_ret["extras"]["psi4:arrays"]["tIjAb"]))
    print("Comparing T1")
    psi4.compare_values(tIA, np.array(json_ret["extras"]["psi4:arrays"]["tIA"]))
    print("Comparing Densities")
    psi4.compare_values(Da, np.array(json_ret["extras"]["psi4:arrays"]["Da"]))

