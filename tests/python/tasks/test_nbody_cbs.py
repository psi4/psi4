import psi4
import json

from psi4.driver.task_base import SingleResult
from psi4.driver.driver_nbody import NBodyComputer
from psi4.driver.driver_cbs import CBSComputer


molecule = psi4.geometry("H 0 0 -2\nH 0 0 -3\n--\nH 0 0 2\nH 0 0 3")

single = {
#    "cbs_metadata": [{"wfn": "mp2", "basis": "cc-pV[D,T]Z"}],
    "cbs_metadata": [{'wfn': 'mp2', 'basis': 'cc-pV[D,T]Z'}],
    "driver": "energy",
    "keywords": {}
}

nb = NBodyComputer(molecule=molecule, driver="energy", bsse_type=["nocp"])

nb.build_tasks(CBSComputer, bsse_type="nocp", **single) 

nb.compute()

print(nb.get_results())
