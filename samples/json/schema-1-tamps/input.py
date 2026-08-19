#! test QCSchema for CCSD amplitudes saving

import numpy as np
import psi4
import sys

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

if sys.version_info >= (3, 14):
    json_ret = psi4.schema_wrapper.run_qcschema(json_data, True, return_dict=True)

    with open("output.json", "w", encoding="utf-8") as ofile:
        from qcelemental.models._v1v2 import AtomicResult

        json_model = AtomicResult(**json_ret)
        ofile.write(json_model.model_dump_json(indent=2))
        psi4.compare(True, True, "json-able")
else:
    json_ret = psi4.schema_wrapper.run_qcschema(json_data, True)

    with open("output.json", "w") as ofile:
        ofile.write(json_ret.json())
        psi4.compare(True, True, "json-able")

    json_ret = json_ret.dict()
psi4.set_options({"basis": basis, "scf_type": "df", "mp2_type": "df"})
e, wfn = psi4.energy(method, return_wfn = True)

tIJAB = wfn.get_amplitudes()["tIjAb"].to_array()
tIA = wfn.get_amplitudes()["tIA"].to_array()
Da = wfn.Da().to_array()

# Make sure the amplitudes compted from h2o_test are actually assigned.
if all(map(lambda x: x == 0, tIJAB.flatten())) or \
   all(map(lambda x: x == 0, tIA.flatten())) or \
   all(map(lambda x: x == 0, Da.flatten())) :
    raise Exception("Error in computing values to compare against. Could not compare.")
else :
    psi4.compare_values(tIJAB, np.array(json_ret["extras"]["psi4:tamps"]["tIjAb"]))
    psi4.compare_values(tIA, np.array(json_ret["extras"]["psi4:tamps"]["tIA"]))
    psi4.compare_values(Da, np.array(json_ret["extras"]["psi4:tamps"]["Da"]))
