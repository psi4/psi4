import pytest

from utils import *

import psi4
import numpy as np
from pathlib import Path
import json
from addons import using


pytestmark = [pytest.mark.psi, pytest.mark.api]


with open(Path(__file__).parent / "test_adcc/adcc_reference_data.json") as f:
    reference_data = json.load(f)


pytestcases = []
for case in reference_data:
    config = case['config']
    marks = [*using('adcc')]
    if config['label'] == 'quick':
        marks.append(pytest.mark.quick)
    # for easier manual test selection
    casename = "_".join([
        config['molname'], config['basis'],
        config['method'].replace("(", "").replace(")", ""),
        config['kind'], str(config['n_states'])
    ])
    pytestcases.append(
        pytest.param(case, marks=marks, id=casename)
    )

@pytest.mark.parametrize('case', pytestcases)
def test_adcc_reference_data(case):
    conf = case["config"]
    psi4.core.clean()
    psi4.set_options({
        'reference': conf['reference'],
        'basis': conf['basis'],
        'guess': 'core',
        'scf_type': 'pk',
        'roots_per_irrep': [conf['n_states']],
        'kind': conf['kind'],
        'num_core_orbitals': conf['core_orbitals'],
        'freeze_core': conf['frozen_core'],
        'num_frozen_uocc': conf['frozen_virtual'],
        'e_convergence': 1e-10,
        'd_convergence': 1e-8,
        'gauge': conf['gauge'],
        'adc__r_convergence': 1e-8,
    })
    mol = psi4.core.Molecule.from_string(conf['molstring'])

    props = ['TRANSITION_DIPOLE', 'DIPOLE', 'OSCILLATOR_STRENGTH', 'ROTATIONAL_STRENGTH']
    en_adc, wfn = psi4.properties(conf['method'], properties=props, return_wfn=True, molecule=mol)

    np.testing.assert_allclose(en_adc, case['energy_gs'], atol=1e-7)
    np.testing.assert_allclose(wfn.variable("CURRENT ENERGY"), case['energy_gs'], atol=1e-7)

    method = conf['method'].upper()
    for i in range(conf['n_states']):
        root_index = i + 1
        prop_access_patterns = {
            'excitation_energy': ([
                f"{method} ROOT 0 -> ROOT {root_index} EXCITATION ENERGY",
                f"{method} ROOT 0 (IN A) -> ROOT {root_index} (IN A) EXCITATION ENERGY",
                f"{method} ROOT 0 (A) -> ROOT {root_index} (A) EXCITATION ENERGY",
                f"{method} ROOT 0 -> ROOT {root_index} EXCITATION ENERGY - A TRANSITION",
                f"ADC ROOT 0 -> ROOT {root_index} EXCITATION ENERGY",
                f"ADC ROOT 0 (IN A) -> ROOT {root_index} (IN A) EXCITATION ENERGY",
                f"ADC ROOT 0 (A) -> ROOT {root_index} (A) EXCITATION ENERGY",
                f"ADC ROOT 0 -> ROOT {root_index} EXCITATION ENERGY - A TRANSITION",
            ], 1e-6),
            'oscillator_strength': ([
                f"{method} ROOT 0 (IN A) -> ROOT {root_index} (IN A) OSCILLATOR STRENGTH (LEN)",
                f"{method} ROOT 0 (A) -> ROOT {root_index} (A) OSCILLATOR STRENGTH (LEN)",
                f"{method} ROOT 0 -> ROOT {root_index} OSCILLATOR STRENGTH (LEN)",
                f"{method} ROOT 0 -> ROOT {root_index} OSCILLATOR STRENGTH (LEN) - A TRANSITION",
                f"ADC ROOT 0 (IN A) -> ROOT {root_index} (IN A) OSCILLATOR STRENGTH (LEN)",
                f"ADC ROOT 0 (A) -> ROOT {root_index} (A) OSCILLATOR STRENGTH (LEN)",
                f"ADC ROOT 0 -> ROOT {root_index} OSCILLATOR STRENGTH (LEN)",
                f"ADC ROOT 0 -> ROOT {root_index} OSCILLATOR STRENGTH (LEN) - A TRANSITION",
            ], 1e-5),
            'oscillator_strength_velocity': ([
                f"{method} ROOT 0 (IN A) -> ROOT {root_index} (IN A) OSCILLATOR STRENGTH (VEL)",
                f"{method} ROOT 0 (A) -> ROOT {root_index} (A) OSCILLATOR STRENGTH (VEL)",
                f"{method} ROOT 0 -> ROOT {root_index} OSCILLATOR STRENGTH (VEL)",
                f"{method} ROOT 0 -> ROOT {root_index} OSCILLATOR STRENGTH (VEL) - A TRANSITION",
                f"ADC ROOT 0 (IN A) -> ROOT {root_index} (IN A) OSCILLATOR STRENGTH (VEL)",
                f"ADC ROOT 0 (A) -> ROOT {root_index} (A) OSCILLATOR STRENGTH (VEL)",
                f"ADC ROOT 0 -> ROOT {root_index} OSCILLATOR STRENGTH (VEL)",
                f"ADC ROOT 0 -> ROOT {root_index} OSCILLATOR STRENGTH (VEL) - A TRANSITION",
            ], 1e-5),
            'rotatory_strength': ([
                f"{method} ROOT 0 (IN A) -> ROOT {root_index} (IN A) ROTATORY STRENGTH (VEL)",
                f"{method} ROOT 0 (A) -> ROOT {root_index} (A) ROTATORY STRENGTH (VEL)",
                f"{method} ROOT 0 -> ROOT {root_index} ROTATORY STRENGTH (VEL)",
                f"{method} ROOT 0 -> ROOT {root_index} ROTATORY STRENGTH (VEL) - A TRANSITION",
                f"ADC ROOT 0 (IN A) -> ROOT {root_index} (IN A) ROTATORY STRENGTH (VEL)",
                f"ADC ROOT 0 (A) -> ROOT {root_index} (A) ROTATORY STRENGTH (VEL)",
                f"ADC ROOT 0 -> ROOT {root_index} ROTATORY STRENGTH (VEL)",
                f"ADC ROOT 0 -> ROOT {root_index} ROTATORY STRENGTH (VEL) - A TRANSITION",
            ], 1e-5),
            'state_dipole_moment': ([
                f"{method} ROOT {root_index} (IN A) DIPOLE MOMENT",
                f"{method} ROOT {root_index} (A) DIPOLE MOMENT",
                f"{method} ROOT {root_index} DIPOLE MOMENT",
                f"{method} ROOT {root_index} DIPOLE MOMENT - A TRANSITION",
                f"ADC ROOT {root_index} (IN A) DIPOLE MOMENT",
                f"ADC ROOT {root_index} (A) DIPOLE MOMENT",
                f"ADC ROOT {root_index} DIPOLE MOMENT",
                f"ADC ROOT {root_index} DIPOLE MOMENT - A TRANSITION",
            ], 5e-3),
        }
        for propname, (vars, atol) in prop_access_patterns.items():
            # Psi lets us select only one gauge, so skip length gauge oscillator strength
            # in case we have selected velocity gauge and vice versa
            if propname == "oscillator_strength" and conf['gauge'] == "velocity" or \
                propname == "oscillator_strength_velocity" and conf['gauge'] == "length":
                continue
            for v in vars:
                ret = wfn.variable(v)
                if isinstance(ret, psi4.core.Matrix):
                    ret = ret.np.flatten()
                np.testing.assert_allclose(ret, case[propname][i], atol=atol, err_msg=f"Result for {v} incorrect.")
