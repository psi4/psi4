import psi4
import json

from psi4.driver.task_base import SingleResult
from psi4.driver.driver_nbody import NBodyComputer
from psi4.driver.driver_cbs import CBSComputer

single = {
    "molecule": psi4.geometry("He"),
    "cbs_metadata": [{"wfn": "mp2", "basis": "cc-pV[D,T]Z"}],
    "driver": "energy",
    "keywords": {}
}

r = CBSComputer(**single)
# print(r)
print(json.dumps(r.plan(), indent=2))
# print(json.dumps(r.compute(), indent=2))

