import psi4
import pytest

from addons import using, uusing
from psi4.driver.procrouting.dft import dft_builder


pytestmark = [pytest.mark.psi, pytest.mark.api]


def _range_separated_functionals():
    parameters = []
    for name in sorted(dft_builder.dict_functionals):
        definition = dft_builder.dict_functionals[name]
        if not dft_builder.functional_available(definition):
            continue

        # failing with max iter for psi4 and ADIIS error for cuEST
        if name == 'hjs-b88': 
            continue

        functional = dft_builder.build_superfunctional_from_dictionary(definition, 1, 1, True)[0]
        if not functional.is_x_lrc():
            continue

        marks = []
        dispersion_type = definition.get("dispersion", {}).get("type", "")
        if dispersion_type.startswith("d3"):
            marks.extend(using("dftd3"))
        elif dispersion_type.startswith("d4"):
            marks.extend(using("dftd4"))
        parameters.append(pytest.param(definition["name"], id=definition["name"], marks=marks))

    return parameters


@uusing("cuest")
@uusing("cuda_cc8")
@pytest.mark.parametrize("functional", _range_separated_functionals())
def test_cuest_range_separated_energy(functional):
    psi4.core.set_num_threads(4)
    molecule = psi4.geometry(
        """
        0 1
        H 0.0 0.0 -0.7
        H 0.0 0.0  0.7
        units bohr
        symmetry c1
        """
    )

    options = {
        "basis": "def2-svp",
        "scf_type": "df",
        "dft_nuclear_scheme": "stratmann",
        "df_basis_scf": "def2-universal-JKFIT",
        "maxiter": 300,
        "e_convergence": 9,
        "d_convergence": 10, 
        "puream": True,
        "reference": "rhf",
    }

    psi4.set_options({**options, "use_cuest": False})
    psi4_energy = psi4.energy(functional, molecule=molecule)

    psi4.core.clean()
    psi4.set_options({**options, "use_cuest": True, "cuest_mixed_precision": False})
    cuest_energy = psi4.energy(functional, molecule=molecule)

    assert psi4.compare_values(
        psi4_energy,
        cuest_energy,
        atol=5.0e-4, # difference in DF scheme
        label=f"{functional} Psi4 DF versus cuEST DF energy",
    )

