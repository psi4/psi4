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

cbs = CBSComputer(**single)
print(cbs)
print(json.dumps(cbs.plan(), indent=2))
cbs.compute()
print(json.dumps(cbs.get_results(), indent=2))

