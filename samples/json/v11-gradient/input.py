#! test c. v1.1 Schema for gradient

import numpy as np
import psi4

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




bench_gradient = np.array([[  0.0 , 0.0 ,   0.4206844],
                           [  0.0 , 0.0 ,  -0.4206844]])


with open("output.dat", "w") as f:
    f.write(json_data["raw_output"]) 
