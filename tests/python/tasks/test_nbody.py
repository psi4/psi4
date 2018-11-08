import psi4
import json

from psi4.driver.task_base import SingleResult
from psi4.driver.driver_nbody import NBodyComputer

single = {
    "molecule": psi4.geometry("He"),
    "driver": "energy",
    "method": "HF",
    "basis": "sto-3g",
    "keywords": {}
}

r = SingleResult(**single)
# print(r)
# print(json.dumps(r.plan(), indent=2))
# print(json.dumps(r.compute(), indent=2))

He2 = psi4.geometry("He 0 0 -2\n--\nHe 0 0 2")
nb = NBodyComputer(molecule=He2, driver="energy", bsse_type=["nocp"])

data = {
    "driver": "energy",
    "method": "HF",
    "basis": "sto-3g",
    "keywords": {}
}
nb.build_tasks(SingleResult, bsse_type="nocp", **data) 

nb.compute()

print(json.dumps(nb.get_results(), indent=2))
