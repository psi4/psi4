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
psi4.compare_strings("STO-3G", json_data["options"]["BASIS"], "Options test")        #TEST
psi4.compare_integers(True, len(json_data["raw_output"]) > 5000, "Output returned")  #TEST
psi4.compare_integers(True, json_data["success"], "Success")                         #TEST


bench_cp_energy = 0.183936053861                                                     #TEST
cenergy = json_data["variables"]["CURRENT ENERGY"]                                   #TEST
psi4.compare_values(bench_cp_energy, cenergy, 6, "SCF CURRENT ENERGY")               #TEST

cenergy = json_data["return_value"]                                                  #TEST
psi4.compare_values(bench_cp_energy, cenergy, 6, "SCF RETURN_VALUE")                 #TEST

return_wfn = "return_wfn" not in json_data["kwargs"]                                 #TEST
psi4.compare_integers(True, return_wfn, "Immutable input")                           #TEST

with open("output.dat", "w") as f:
    f.write(json_data["raw_output"]) 
