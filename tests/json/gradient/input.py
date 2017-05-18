import numpy as np
import psi4
from psi4.driver import p4util

psi4.set_output_file("output.dat", False)

# Generate JSON data
json_data = {}
json_data["molecule"] = """He 0 0 0\n--\nHe 0 0 1"""
json_data["driver"] = "gradient"
json_data["method"] = 'SCF'
json_data["kwargs"] = {}
json_data["options"] = {"BASIS": "STO-3G"}
json_data["return_output"] = True

psi4.json_wrapper.run_json(json_data)

p4util.compare_strings("STO-3G", json_data["options"]["BASIS"], "Options test") # TEST
p4util.compare_integers(True, json_data["success"], "Success")                  # TEST


bench_energy = -5.433191881443323                                               # TEST
cenergy = json_data["variables"]["CURRENT ENERGY"]                              # TEST
p4util.compare_values(bench_energy, cenergy, 6, "SCF CURRENT ENERGY")           # TEST

bench_gradient = np.array([[  0.0 , 0.0 ,   0.4206844],
                           [  0.0 , 0.0 ,  -0.4206844]])
cgradient = psi4.core.Matrix.from_serial(json_data["return_value"])                  # TEST
p4util.compare_arrays(bench_gradient, cgradient.np, 4, "SCF RETURN_VALUE")      # TEST

return_wfn = "return_wfn" not in json_data["kwargs"]                            # TEST
p4util.compare_integers(True, return_wfn, "Immutable input")                    # TEST

with open("output.dat", "w") as f:
    f.write(json_data["raw_output"]) 
