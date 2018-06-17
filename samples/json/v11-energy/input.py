#! test c. v1.1 Schema

import numpy as np
import psi4

# Generate JSON data
json_data = {}
json_data["molecule"] = """He 0 0 0\n--\nHe 0 0 1"""
json_data["memory"] = "5GB"
json_data["driver"] = "energy"
json_data["method"] = 'SCF'
json_data["kwargs"] = {"bsse_type": "cp"}
json_data["options"] = {"BASIS": "STO-3G"}
json_data["return_output"] = True

psi4.json_wrapper.run_json(json_data)





with open("output.dat", "w") as f:
    f.write(json_data["raw_output"]) 
