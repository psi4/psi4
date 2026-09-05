"""Tests for Python-side DIIS helpers."""

import numpy as np
import pytest

from psi4 import core
from psi4.driver.procrouting import diis


pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.quick]


def _vector(values):
    return core.Vector.from_array(np.asarray(values, dtype=float))


def test_chunked_vector_equivalence():
    """Multiple DIIS items must equal their concatenated logical vector."""
    # Two slots plus three entries also exercises replacement without retaining
    # both the displaced and replacement entries at once.
    chunked = diis.DIIS(2, "chunked test", storage_policy=diis.StoragePolicy.InCore)
    monolithic = diis.DIIS(2, "monolithic test", storage_policy=diis.StoragePolicy.InCore)

    chunk_template = [_vector(np.zeros(3)), _vector(np.zeros(2))]
    flat_template = _vector(np.zeros(5))
    chunked.set_error_vector_size_from_list(chunk_template)
    chunked.set_vector_size_from_list(chunk_template)
    monolithic.set_error_vector_size(flat_template)
    monolithic.set_vector_size(flat_template)

    residuals = (
        np.array([1.0, -0.2, 0.3, 0.7, -0.1]),
        np.array([0.4, -0.1, 0.2, 0.2, -0.05]),
        np.array([0.2, -0.04, 0.08, 0.1, -0.02]),
    )
    targets = (
        np.array([2.0, 3.0, -1.0, 0.5, 4.0]),
        np.array([2.2, 2.8, -0.8, 0.7, 3.7]),
        np.array([2.3, 2.7, -0.7, 0.8, 3.6]),
    )

    for residual, target in zip(residuals, targets):
        chunked.add_entry_from_lists(
            [_vector(residual[:3]), _vector(residual[3:])],
            [_vector(target[:3]), _vector(target[3:])],
        )
        monolithic.add_entry(_vector(residual), _vector(target))

    chunked_result = [_vector(np.zeros(3)), _vector(np.zeros(2))]
    flat_result = _vector(np.zeros(5))
    chunked.extrapolate_from_list(chunked_result)
    monolithic.extrapolate(flat_result)

    joined_result = np.concatenate([chunk.np for chunk in chunked_result])
    assert np.allclose(joined_result, flat_result.np, rtol=0.0, atol=1.0e-12)
