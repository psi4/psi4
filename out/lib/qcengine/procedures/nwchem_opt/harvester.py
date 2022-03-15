import re
from decimal import Decimal
from typing import List, Tuple

from qcelemental.models import AtomicResult, Molecule, OptimizationInput, Provenance
from qcelemental.util import unnp

from qcengine.programs.nwchem.harvester import harvest_outfile_pass
from qcengine.programs.qcvar_identities_resources import build_atomicproperties, build_out
from qcengine.programs.util import PreservingDict


def harvest_output(outtext: str) -> Tuple[List[PreservingDict], List[Molecule], List[list], str, str]:
    """Function to read an entire NWChem output file.

    Reads all of the different "line search" segments of a file and returns
    values from up to last segment for which a geometry was written.

    Args:
        outtext (str): Output written to stdout
    Returns:
        - Variables extracted from the output file in each optimization step
        - Molecules from the each complete step
        - Gradients from the each complete step
        - (str): Version string
        - (str): Error message, if any
    """

    # Loop over all steps
    # TODO (wardlt): Is it only necessary to read the last two steps?
    pass_psivar = []
    pass_coord = []
    pass_grad = []
    version = error = None
    for outpass in re.split(r"Step +\d", outtext, re.MULTILINE)[1:]:
        psivar, nwcoord, nwgrad, version, module, error = harvest_outfile_pass(outpass)
        pass_psivar.append(psivar)
        pass_coord.append(nwcoord)
        pass_grad.append(nwgrad)

    # Ensure last segment contains the last geometry
    if not pass_coord[-1]:
        pass_psivar.pop()
        pass_coord.pop()
        pass_grad.pop()

    return pass_psivar, pass_coord, pass_grad, version, module, error


def harvest_as_atomic_result(input_model: OptimizationInput, nwout: str) -> List[AtomicResult]:
    """Parse each step in the geometry relaxation as a separate AtomicResult

    Args:
        input_model: Input specification for the relaxation
        nwout: Standard out from the NWChem simulation
    Returns:
        A list of the results at each step
    """
    # Parse the files
    out_psivars, out_mols, out_grads, version, module, error = harvest_output(nwout)

    # Make atomic results
    results = []
    for qcvars, nwgrad, out_mol in zip(out_psivars, out_grads, out_mols):
        if nwgrad is not None:
            qcvars[f"{input_model.input_specification.model.method.upper()[4:]} TOTAL GRADIENT"] = nwgrad
            qcvars["CURRENT GRADIENT"] = nwgrad

        # Get the formatted properties
        build_out(qcvars)
        atprop = build_atomicproperties(qcvars)

        provenance = Provenance(creator="NWChem", version=version, routine="nwchem_opt").dict()
        if module is not None:
            provenance["module"] = module

        # Format them inout an output
        output_data = {
            "schema_version": 1,
            "molecule": out_mol,
            "driver": "gradient",
            "extras": input_model.extras.copy(),
            "model": input_model.input_specification.model,
            "keywords": input_model.input_specification.keywords,
            "properties": atprop,
            "provenance": provenance,
            "return_result": nwgrad,
            "success": True,
        }

        # got to even out who needs plump/flat/Decimal/float/ndarray/list
        # Decimal --> str preserves precision
        output_data["extras"]["qcvars"] = {
            k.upper(): str(v) if isinstance(v, Decimal) else v for k, v in unnp(qcvars, flat=True).items()
        }

        results.append(AtomicResult(**output_data))
    return results
