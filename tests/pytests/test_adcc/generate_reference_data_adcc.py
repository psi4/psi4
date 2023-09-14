import adcc
import psi4

import numpy as np
import json
from dataclasses import dataclass, asdict


molecules = {
    "cn": """
        0 2
        C 0 0 0
        N 0 0 2.2143810738114829
        symmetry c1
        units au
    """,
    "formaldehyde": """
        C 2.0092420208996 3.8300915804899 0.8199294419789
        O 2.1078857690998 2.0406638776593 2.1812021228452
        H 2.0682421748693 5.7438044586615 1.5798996515014
        H 1.8588483602149 3.6361694243085 -1.2192956060942
        symmetry c1
        units au
        no_reorient
        no_com
    """,
    "h2o": """
        O 0 0 0
        H 0 0 1.795239827225189
        H 1.693194615993441 0 -0.599043184453037
        symmetry c1
        units au
    """,
    "hf": """
        0 3
        H 0 0 0
        F 0 0 2.5
        symmetry c1
        units au
    """,
    "methyloxirane": """
        O	0.7971066654	0.9044360742	0.0836962049
        C	-0.1867183086	-0.0290724859	0.5536827176
        C	-1.4336843546	-0.1726679227	-0.2822214295
        C	1.1302222000	-0.4892393880	0.0894444115
        H	1.2197487995	-0.9517340291	-0.8946449424
        H	1.8923895176	-0.7869225283	0.8107731933
        H	-0.3474086480	0.0162374592	1.6337796505
        H	-2.0955293870	0.6891134744	-0.1384941617
        H	-1.9883466588	-1.0759327249	-0.0005360999
        H	-1.1805969868	-0.2349473270	-1.3455182514
        symmetry c1
    """,
}


@dataclass
class AdccConfig:
    method: str
    basis: str
    molstring: str
    molname: str
    n_states: int
    kind: str
    core_orbitals: int = 0
    reference: str = 'rhf'
    frozen_core: bool = False
    frozen_virtual: int = 0
    gauge: str = 'length'
    label: str = ''


if __name__ == "__main__":
    psi4.core.be_quiet()
    configs = [
        AdccConfig('adc(1)', 'cc-pvdz', molecules['cn'], 'cn', n_states=4, kind='any', reference='uhf', label='quick'),  # cn-adc1
        AdccConfig('adc(2)', 'cc-pvdz', molecules['cn'], 'cn', n_states=6, kind='any', reference='uhf'),  # cn-adc2
        AdccConfig('adc(1)', 'cc-pvdz', molecules['h2o'], 'h2o', n_states=10, kind='singlet', label='quick'),  # h2o-adc1
        AdccConfig('adc(2)', 'cc-pvdz', molecules['h2o'], 'h2o', n_states=10, kind='singlet'),  # h2o-adc2
        AdccConfig('adc(1)', 'cc-pvdz', molecules['h2o'], 'h2o', n_states=2, kind='singlet', label='quick'),  # h2o-adc1-quick
        AdccConfig('adc(2)', 'cc-pvdz', molecules['h2o'], 'h2o', n_states=6, kind='any'),  # h2o-adc2-any
        AdccConfig('adc(2)-x', 'cc-pvdz', molecules['h2o'], 'h2o', n_states=10, kind='singlet'),  # h2o-adc2x
        AdccConfig('adc(3)', 'cc-pvdz', molecules['h2o'], 'h2o', n_states=10, kind='singlet'),  # h2o-adc3
        AdccConfig('cvs-adc(1)', 'cc-pvdz', molecules['h2o'], 'h2o', n_states=10, kind='singlet', core_orbitals=1),  # h2o-cvs-adc1
        AdccConfig('cvs-adc(2)', 'cc-pvdz', molecules['h2o'], 'h2o', n_states=10, kind='singlet', core_orbitals=1),  # h2o-cvs-adc2
        AdccConfig('cvs-adc(2)-x', 'cc-pvdz', molecules['h2o'], 'h2o', n_states=10, kind='singlet', core_orbitals=1),  # h2o-cvs-adc2x
        AdccConfig('cvs-adc(2)-x', 'cc-pvdz', molecules['h2o'], 'h2o', n_states=7, kind='triplet', core_orbitals=1, label='quick'),  # h2o-cvs-adc2x-triplets
        AdccConfig('cvs-adc(3)', 'cc-pvdz', molecules['h2o'], 'h2o', n_states=10, kind='singlet', core_orbitals=1),  # h2o-cvs-adc3
        AdccConfig('adc(2)', 'cc-pvdz', molecules['h2o'], 'h2o', n_states=5, kind='singlet', frozen_core=True, label='quick'),  # h2o-fc-adc2
        AdccConfig('adc(2)-x', 'cc-pvdz', molecules['h2o'], 'h2o', n_states=5, kind='singlet', frozen_virtual=5, label='quick'),  # h2o-fv-adc2x
        AdccConfig('adc(2)', '6-31G', molecules['hf'], 'hf', n_states=5, kind='spin_flip', reference='uhf'),  # hf-adc2-spin-flip
        AdccConfig('adc(2)', 'sto-3g', molecules['methyloxirane'], 'methyloxirane', n_states=5, kind='singlet', gauge='velocity', label='quick')  # methyloxirane
    ]
    
    ret = []
    for conf in configs:
        case = {}
        case['config'] = asdict(conf)
        psi4.core.clean()
        
        psi4.set_options({
            'reference': conf.reference,
            'basis': conf.basis,
            'guess': 'core',
            'scf_type': 'pk',
            'freeze_core': conf.frozen_core,
            'num_frozen_uocc': conf.frozen_virtual,
            'e_convergence': 1e-10,
            'd_convergence': 1e-8,
            # 'roots_per_irrep': [conf.n_states],
            # 'kind': conf.kind,
        })
        mol = psi4.core.Molecule.from_string(conf.molstring)
        en_scf, wfn = psi4.energy('hf', molecule=mol, return_wfn=True)
        adc_name = conf.method.replace("(", "").replace(")", "")
        adc_name = adc_name.replace("-x", "x")
        
        scf_accuracy = max(psi4.core.get_option("SCF", "E_CONVERGENCE"),
                       psi4.core.get_option("SCF", "D_CONVERGENCE"))
        # if psi4.core.get_option("ADC", "R_CONVERGENCE") < 0:
        #     conv_tol = max(100 * scf_accuracy, 1e-6)
        # else:
        #     conv_tol = psi4.core.get_option("ADC", "R_CONVERGENCE")
        conv_tol = 1e-8
        
        fc = 0
        if wfn.frzcpi()[0] > 0:
            assert conf.frozen_core
            fc = wfn.frzcpi()[0]
        fv = 0
        if wfn.frzvpi()[0] > 0:
            fv = wfn.frzvpi()[0]
        state = adcc.run_adc(
            wfn, method=adc_name, n_states=conf.n_states,
            kind=conf.kind, core_orbitals=conf.core_orbitals,
            conv_tol=conv_tol, frozen_core=fc, frozen_virtual=fv
        )
        
        props = [
            'excitation_energy',
            'oscillator_strength',
            'oscillator_strength_velocity',
            'state_dipole_moment',
            'transition_dipole_moment',
            'rotatory_strength']
        for pr in props:
            case[pr] = getattr(state, pr).tolist()
        
        is_cvs_adc3 = state.method.level >= 3 and state.ground_state.has_core_occupied_space
        if is_cvs_adc3:
            energy_gs = state.ground_state.energy(2)
        else:
            energy_gs = state.ground_state.energy(state.method.level)
        case['energy_gs'] = energy_gs
        ret.append(case)
        
    with open("adcc_reference_data.json", "w") as outfile:
        json.dump(ret, outfile)