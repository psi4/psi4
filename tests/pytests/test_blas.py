import numpy as np
import pytest

import psi4

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.quick]


@pytest.mark.parametrize("values,expected", [
    ([1.0], 0),
    ([1.0, -3.0, 2.0], 1),
    ([-5.0, 1.0, 2.0], 0),
    ([1.0, 2.0, -7.5], 2),
    ([0.5, 0.5, 0.5, 0.5], 0),
])
def test_idamax(values, expected):
    """IDAMAX returns the zero-based position of the largest absolute value."""
    vec = psi4.core.Vector.from_array(np.array(values))
    assert psi4.core.IDAMAX(0, len(values), vec, 1) == expected


def test_idamax_random():
    rng = np.random.default_rng(20260910)
    for _ in range(50):
        arr = rng.uniform(-10.0, 10.0, size=rng.integers(1, 200))
        vec = psi4.core.Vector.from_array(arr)
        assert psi4.core.IDAMAX(0, arr.size, vec, 1) == int(np.argmax(np.abs(arr)))
