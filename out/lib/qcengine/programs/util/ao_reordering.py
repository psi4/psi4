import sys
import traceback
from typing import Any, Dict, List

import numpy as np
from qcelemental.models import BasisSet


def cca_ao_order_spherical(max_angular_momentum: int) -> Dict[int, List[int]]:
    ao_order = {}
    for ang_mom in range(max_angular_momentum):
        ao_order[ang_mom] = [x for x in range(-1 * ang_mom, ang_mom + 1)]
    return ao_order


def get_ao_conversion(new_ao_order: Dict[int, List[int]], old_ao_order: Dict[int, List[int]]) -> Dict[int, List[int]]:
    """
    This function determines the conversion of AO ordering for each angular momentum between the CCA standard and the
    provided program harness.
    """

    to_new_order = {}
    for ang_mom in old_ao_order.keys():
        old_order = old_ao_order[ang_mom]
        new_order = new_ao_order[ang_mom]
        sort_indices = []
        for element in new_order:
            sort_indices.append(old_order.index(element))

        to_new_order[ang_mom] = sort_indices

    return to_new_order


def reorder_column_ao_indices(
    matrix: np.ndarray, basis: BasisSet, to_new_ao_order: Dict[str, Dict[int, List[int]]]
) -> np.ndarray:
    """
    Reorder the column atomic orbital (AO) indices of matrix.
    """

    num_cols = len(matrix[0])
    num_rows = len(matrix)
    new_matrix = np.zeros([num_rows, num_cols])
    for row_index in range(0, num_rows):
        ao_shift = 0
        for atom in basis.atom_map:
            for electron_shell in basis.center_data[atom].electron_shells:
                to_cca_indices = to_new_ao_order[electron_shell.harmonic_type][electron_shell.angular_momentum[0]]
                for bas_index in range(0, len(to_cca_indices)):
                    new_matrix[row_index][ao_shift + to_cca_indices[bas_index]] = matrix[row_index][
                        ao_shift + bas_index
                    ]
                ao_shift += len(to_cca_indices)

    return new_matrix


def reorder_row_and_column_ao_indices(
    matrix: np.ndarray, basis: BasisSet, to_new_ao_order: Dict[str, Dict[int, List[int]]]
) -> np.ndarray:
    """
    Reorder both row and column atomic orbital (AO) indices of matrix. (e.g. Fock, density)
    """

    col_reordered_matrix = reorder_column_ao_indices(matrix, basis, to_new_ao_order)
    new_matrix = reorder_column_ao_indices(col_reordered_matrix.transpose(), basis, to_new_ao_order)

    return new_matrix


def mill_qcvars(mill: "AlignmentMill", qcvars: Dict[str, Any]) -> Dict[str, Any]:
    """Apply translation, rotation, atom shuffle defined in ``mill`` to the nonscalar quantities in ``qcvars``."""

    milled = {}
    for k, v in qcvars.items():
        if k.endswith("GRADIENT"):
            milled[k] = mill.align_gradient(v)
        elif k.endswith("HESSIAN"):
            milled[k] = mill.align_hessian(v)
        else:
            milled[k] = v

    return milled


def error_stamp(stdout: str = "", stderr: str = "", tb: str = None) -> str:
    """Return all useful information in error string."""

    if not tb:
        tb = traceback.format_exception(*sys.exc_info())
    return "STDOUT:\n" + stdout + "\nSTDERR:\n" + stderr + "\nTRACEBACK:\n" + "".join(tb)
