"""
Interaction-induced one-electron properties of dimers from SAPT0 and SAPT(DFT),
``psi4.properties("sapt0" | "sapt(dft)")``: the monomer values and the first-order
polarization and exchange terms.

The reference values come from propSAPT v0.5.0 (https://github.com/btyrcha/propSAPT),
whose first-order terms of each operator component were assembled into these
QCVariables; the port agrees with it to 1e-10 a.u. on the same Psi4.
"""

import numpy as np
import pytest

import psi4
from psi4 import compare_values

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.sapt]

HE_H2 = """
symmetry c1
no_com
no_reorient
units bohr
0 1
He  -2.500000000   0.000000000   0.000000000
--
0 1
H    1.775500000   0.000000000   0.000000000
H    3.224500000   0.000000000   0.000000000
"""

WATER_DIMER = """
symmetry c1
no_com
no_reorient
0 1
O   -1.41606479    0.00000000    0.00000000
H   -0.52004811    0.00000000   -0.08624091
H   -1.81884655    0.00000000   -0.89989100
--
0 1
O    1.41606479    0.00000000    0.00000000
H    1.71179607   -0.76501221    0.50083629
H    1.71179607    0.76501221    0.50083629
"""

CONVERGENCE = {"e_convergence": 1e-10, "d_convergence": 1e-11, "cphf_r_convergence": 1e-10}

# He-H2, SAPT0/aug-cc-pVDZ, SCF_TYPE DIRECT
REF_HE_H2_SAPT0 = {
    'SAPT MONOMER A DIPOLE': [-0.000666424666, -0.0, 0.0],
    'SAPT MONOMER B DIPOLE': [-3.7557148e-05, 0.0, 0.0],
    'SAPT POL10,R DIPOLE': [-0.000214727451, 0.0, -0.0],
    'SAPT POL10,R DIPOLE - A': [-0.004535069512, -0.0, -0.0],
    'SAPT POL10,R DIPOLE - B': [0.004320342061, 0.0, -0.0],
    'SAPT EXCH10(S^2) DIPOLE': [-0.016025621138, -0.0, 0.0],
    'SAPT EXCH10(S^2) DIPOLE - A': [0.005959150125, -0.0, 0.0],
    'SAPT EXCH10(S^2) DIPOLE - B': [-0.021984771263, -0.0, 0.0],
    'SAPT EXCH10 DIPOLE': [-0.016052648008, -0.0, 0.0],
    'SAPT EXCH10 DIPOLE - A': [0.005969362232, -0.0, 0.0],
    'SAPT EXCH10 DIPOLE - B': [-0.02202201024, -0.0, 0.0],
    'SAPT MONOMER A QUADRUPOLE': [[-0.793369038733, -0.0, -0.0], [-0.0, -0.798317558007, 0.0],
                                  [-0.0, 0.0, -0.798317558007]],
    'SAPT MONOMER B QUADRUPOLE': [[-1.075117252717, -0.0, 0.0], [-0.0, -1.564897242189, -0.0],
                                  [0.0, -0.0, -1.564897242189]],
    'SAPT POL10,R QUADRUPOLE': [[0.029805473662, 0.0, -0.0], [0.0, -7.17547e-05, 0.0], [-0.0, 0.0, -7.17547e-05]],
    'SAPT POL10,R QUADRUPOLE - A': [[0.020173323423, 0.0, 0.0], [0.0, 0.000155718774, -0.0],
                                    [0.0, -0.0, 0.000155718774]],
    'SAPT POL10,R QUADRUPOLE - B': [[0.009632150239, 0.0, -0.0], [0.0, -0.000227473474, 0.0],
                                    [-0.0, 0.0, -0.000227473474]],
    'SAPT EXCH10(S^2) QUADRUPOLE': [[-0.057776521313, -0.0, 0.0], [-0.0, 0.005465853993, 0.0],
                                    [0.0, 0.0, 0.005465853993]],
    'SAPT EXCH10(S^2) QUADRUPOLE - A': [[-0.01389998052, -0.0, 0.0], [-0.0, 0.002295703624, 0.0],
                                        [0.0, 0.0, 0.002295703624]],
    'SAPT EXCH10(S^2) QUADRUPOLE - B': [[-0.043876540793, -0.0, 0.0], [-0.0, 0.003170150369, -0.0],
                                        [0.0, -0.0, 0.003170150369]],
    'SAPT EXCH10 QUADRUPOLE': [[-0.057880706346, -0.0, 0.0], [-0.0, 0.005475807242, 0.0], [0.0, 0.0, 0.005475807242]],
    'SAPT EXCH10 QUADRUPOLE - A': [[-0.013925994723, -0.0, 0.0], [-0.0, 0.002299989748, 0.0],
                                   [0.0, 0.0, 0.002299989748]],
    'SAPT EXCH10 QUADRUPOLE - B': [[-0.043954711623, -0.0, 0.0], [-0.0, 0.003175817494, -0.0],
                                   [0.0, -0.0, 0.003175817494]],
    'SAPT MONOMER A OCTUPOLE': [[[5.953429697466, 0.0, 0.0], [0.0, 1.992559035438, 0.0], [0.0, 0.0, 1.992559035438]],
                                [[0.0, 1.992559035438, 0.0], [1.992559035438, -0.0, 0.0], [0.0, 0.0, -0.0]],
                                [[0.0, 0.0, 1.992559035438], [0.0, 0.0, -0.0], [1.992559035438, -0.0, 0.0]]],
    'SAPT MONOMER B OCTUPOLE': [[[-8.076795996831, -0.0, 0.0], [-0.0, -3.908630406444, 0.0],
                                 [0.0, 0.0, -3.908630406444]],
                                [[-0.0, -3.908630406444, 0.0], [-3.908630406444, 0.0, 0.0], [0.0, 0.0, 0.0]],
                                [[0.0, 0.0, -3.908630406444], [0.0, 0.0, 0.0], [-3.908630406444, 0.0, 0.0]]],
    'SAPT POL10,R OCTUPOLE': [[[-0.036560745934, 0.0, -0.0], [0.0, 0.000545288859, 0.0], [-0.0, 0.0, 0.000545288859]],
                              [[0.0, 0.000545288859, 0.0], [0.000545288859, 0.0, -0.0], [0.0, -0.0, 0.0]],
                              [[-0.0, 0.0, 0.000545288859], [0.0, -0.0, 0.0], [0.000545288859, 0.0, -0.0]]],
    'SAPT POL10,R OCTUPOLE - A': [[[-0.082388170843, -0.0, -0.0], [-0.0, -0.003893261094, -0.0],
                                   [-0.0, -0.0, -0.003893261094]],
                                  [[-0.0, -0.003893261094, -0.0], [-0.003893261094, -0.0, -0.0], [-0.0, -0.0, -0.0]],
                                  [[-0.0, -0.0, -0.003893261094], [-0.0, -0.0, -0.0], [-0.003893261094, -0.0, -0.0]]],
    'SAPT POL10,R OCTUPOLE - B': [[[0.04582742491, 0.0, -0.0], [0.0, 0.004438549954, 0.0], [-0.0, 0.0,
                                                                                            0.004438549954]],
                                  [[0.0, 0.004438549954, 0.0], [0.004438549954, 0.0, -0.0], [0.0, -0.0, 0.0]],
                                  [[-0.0, 0.0, 0.004438549954], [0.0, -0.0, 0.0], [0.004438549954, 0.0, -0.0]]],
    'SAPT EXCH10(S^2) OCTUPOLE': [[[-0.175573981562, -0.0, 0.0], [-0.0, -0.022172037679, 0.0],
                                   [0.0, 0.0, -0.022172037679]],
                                  [[-0.0, -0.022172037679, 0.0], [-0.022172037679, -0.0, 0.0], [0.0, 0.0, -0.0]],
                                  [[0.0, 0.0, -0.022172037679], [0.0, 0.0, -0.0], [-0.022172037679, -0.0, 0.0]]],
    'SAPT EXCH10(S^2) OCTUPOLE - A': [[[0.05268817706, -0.0, 0.0], [-0.0, 0.001985893714, 0.0],
                                       [0.0, 0.0, 0.001985893714]],
                                      [[-0.0, 0.001985893714, 0.0], [0.001985893714, -0.0, 0.0], [0.0, 0.0, -0.0]],
                                      [[0.0, 0.0, 0.001985893714], [0.0, 0.0, -0.0], [0.001985893714, -0.0, 0.0]]],
    'SAPT EXCH10(S^2) OCTUPOLE - B': [[[-0.228262158622, -0.0, 0.0], [-0.0, -0.024157931393, 0.0],
                                       [0.0, 0.0, -0.024157931393]],
                                      [[-0.0, -0.024157931393, 0.0], [-0.024157931393, -0.0, 0.0], [0.0, 0.0, -0.0]],
                                      [[0.0, 0.0, -0.024157931393], [0.0, 0.0, -0.0], [-0.024157931393, -0.0, 0.0]]],
    'SAPT EXCH10 OCTUPOLE': [[[-0.175864578797, -0.0, 0.0], [-0.0, -0.022210318452, 0.0], [0.0, 0.0, -0.022210318452]],
                             [[-0.0, -0.022210318452, 0.0], [-0.022210318452, -0.0, 0.0], [0.0, 0.0, -0.0]],
                             [[0.0, 0.0, -0.022210318452], [0.0, 0.0, -0.0], [-0.022210318452, -0.0, 0.0]]],
    'SAPT EXCH10 OCTUPOLE - A': [[[0.052781780981, -0.0, 0.0], [-0.0, 0.001988185152, 0.0], [0.0, 0.0,
                                                                                             0.001988185152]],
                                 [[-0.0, 0.001988185152, 0.0], [0.001988185152, -0.0, 0.0], [0.0, 0.0, -0.0]],
                                 [[0.0, 0.0, 0.001988185152], [0.0, 0.0, -0.0], [0.001988185152, -0.0, 0.0]]],
    'SAPT EXCH10 OCTUPOLE - B': [[[-0.228646359778, -0.0, 0.0], [-0.0, -0.024198503603, 0.0],
                                  [0.0, 0.0, -0.024198503603]],
                                 [[-0.0, -0.024198503603, 0.0], [-0.024198503603, -0.0, 0.0], [0.0, 0.0, -0.0]],
                                 [[0.0, 0.0, -0.024198503603], [0.0, 0.0, -0.0], [-0.024198503603, -0.0, 0.0]]],
    'SAPT MONOMER A ESP AT CENTER 1':
    -3.361382870621,
    'SAPT MONOMER B ESP AT CENTER 1':
    0.004133157124,
    'SAPT POL10,R ESP AT CENTER 1':
    0.000106655344,
    'SAPT POL10,R ESP AT CENTER 1 - A':
    0.000522242428,
    'SAPT POL10,R ESP AT CENTER 1 - B':
    -0.000415587084,
    'SAPT EXCH10(S^2) ESP AT CENTER 1':
    -0.000706402372,
    'SAPT EXCH10(S^2) ESP AT CENTER 1 - A':
    -0.002808647429,
    'SAPT EXCH10(S^2) ESP AT CENTER 1 - B':
    0.002102245057,
    'SAPT EXCH10 ESP AT CENTER 1':
    -0.000708497287,
    'SAPT EXCH10 ESP AT CENTER 1 - A':
    -0.002814116834,
    'SAPT EXCH10 ESP AT CENTER 1 - B':
    0.002105619547,
    'SAPT MONOMER A ESP AT CENTER 2':
    -5.917675e-06,
    'SAPT MONOMER B ESP AT CENTER 2':
    -1.080195502955,
    'SAPT POL10,R ESP AT CENTER 2':
    7.1595918e-05,
    'SAPT POL10,R ESP AT CENTER 2 - A':
    -0.000321601819,
    'SAPT POL10,R ESP AT CENTER 2 - B':
    0.000393197737,
    'SAPT EXCH10(S^2) ESP AT CENTER 2':
    -0.00160438739,
    'SAPT EXCH10(S^2) ESP AT CENTER 2 - A':
    0.000623063627,
    'SAPT EXCH10(S^2) ESP AT CENTER 2 - B':
    -0.002227451017,
    'SAPT EXCH10 ESP AT CENTER 2':
    -0.001607286981,
    'SAPT EXCH10 ESP AT CENTER 2 - A':
    0.000624070302,
    'SAPT EXCH10 ESP AT CENTER 2 - B':
    -0.002231357284,
    'SAPT MONOMER A ESP AT CENTER 3':
    -8.408828e-06,
    'SAPT MONOMER B ESP AT CENTER 3':
    -1.080346961277,
    'SAPT POL10,R ESP AT CENTER 3':
    0.000614436384,
    'SAPT POL10,R ESP AT CENTER 3 - A':
    -0.000163675551,
    'SAPT POL10,R ESP AT CENTER 3 - B':
    0.000778111936,
    'SAPT EXCH10(S^2) ESP AT CENTER 3':
    -0.003615947949,
    'SAPT EXCH10(S^2) ESP AT CENTER 3 - A':
    0.000326364675,
    'SAPT EXCH10(S^2) ESP AT CENTER 3 - B':
    -0.003942312624,
    'SAPT EXCH10 ESP AT CENTER 3':
    -0.003622308246,
    'SAPT EXCH10 ESP AT CENTER 3 - A':
    0.000326885089,
    'SAPT EXCH10 ESP AT CENTER 3 - B':
    -0.003949193335,
    'SAPT MONOMER A MULLIKEN CHARGES': [-1.514676e-06, -0.009460828631, 0.009462343307],
    'SAPT MONOMER B MULLIKEN CHARGES': [0.002991947833, 0.009281786833, -0.012273734666],
    'SAPT POL10,R MULLIKEN CHARGES': [0.00014463965, 0.003430380682, -0.003575020332],
    'SAPT POL10,R MULLIKEN CHARGES - A': [0.000741352294, -0.000164531442, -0.000576820852],
    'SAPT POL10,R MULLIKEN CHARGES - B': [-0.000596712644, 0.003594912125, -0.00299819948],
    'SAPT EXCH10(S^2) MULLIKEN CHARGES': [0.002230418664, 0.006644802607, -0.008875221271],
    'SAPT EXCH10(S^2) MULLIKEN CHARGES - A': [-0.001902943944, 0.005345452533, -0.003442508589],
    'SAPT EXCH10(S^2) MULLIKEN CHARGES - B': [0.004133362608, 0.001299350074, -0.005432712682],
    'SAPT EXCH10 MULLIKEN CHARGES': [0.002233628107, 0.006642777429, -0.008876405537],
    'SAPT EXCH10 MULLIKEN CHARGES - A': [-0.001907036078, 0.005353426777, -0.003446390699],
    'SAPT EXCH10 MULLIKEN CHARGES - B': [0.004140664185, 0.001289350652, -0.005430014838],
    'SAPT MONOMER A LOWDIN CHARGES': [0.031800498127, -0.029186902364, -0.002613595763],
    'SAPT MONOMER B LOWDIN CHARGES': [-0.011993091586, 0.005265171033, 0.006727920554],
    'SAPT POL10,R LOWDIN CHARGES': [-9.2277485e-05, -9.3014725e-05, 0.00018529221],
    'SAPT POL10,R LOWDIN CHARGES - A': [0.000520709988, -0.000460763459, -5.994653e-05],
    'SAPT POL10,R LOWDIN CHARGES - B': [-0.000612987473, 0.000367748734, 0.000245238739],
    'SAPT EXCH10(S^2) LOWDIN CHARGES': [0.001754269752, -0.000878385841, -0.00087588391],
    'SAPT EXCH10(S^2) LOWDIN CHARGES - A': [-0.001549558554, 0.000681714985, 0.000867843569],
    'SAPT EXCH10(S^2) LOWDIN CHARGES - B': [0.003303828305, -0.001560100826, -0.001743727479],
    'SAPT EXCH10 LOWDIN CHARGES': [0.001756878871, -0.000879567097, -0.000877311774],
    'SAPT EXCH10 LOWDIN CHARGES - A': [-0.001552499202, 0.000683429975, 0.000869069228],
    'SAPT EXCH10 LOWDIN CHARGES - B': [0.003309378073, -0.001562997072, -0.001746381001]
}
REF_HE_H2_SAPT0_ENERGY = 0.000638971671

# water dimer, SAPT(DFT)/cc-pVDZ, PBE0 with the GRAC shifts of test_saptdft.py, SCF_TYPE DF
REF_WATER_SAPTDFT = {
    'SAPT MONOMER A DIPOLE': [0.346057973577, 0.0, -0.659044801174],
    'SAPT MONOMER B DIPOLE': [0.453579823242, 0.0, 0.660148028415],
    'SAPT POL1,R DIPOLE': [0.150341387086, 0.0, 0.032247870966],
    'SAPT POL1,R DIPOLE - A': [0.044633843585, 0.0, -0.002501384728],
    'SAPT POL1,R DIPOLE - B': [0.105707543501, 0.0, 0.034749255694],
    'SAPT EXCH1(S^2) DIPOLE': [-0.026277786946, -0.0, -0.00977538475],
    'SAPT EXCH1(S^2) DIPOLE - A': [0.090713980043, 0.0, 0.001586441409],
    'SAPT EXCH1(S^2) DIPOLE - B': [-0.116991766989, -0.0, -0.011361826159],
    'SAPT EXCH1 DIPOLE': [-0.026727751611, -0.0, -0.009967586018],
    'SAPT EXCH1 DIPOLE - A': [0.092118808655, 0.0, 0.001603920953],
    'SAPT EXCH1 DIPOLE - B': [-0.118846560266, -0.0, -0.011571506971],
    'SAPT MONOMER A QUADRUPOLE': [[-5.264710148982, -0.0, 2.296178587363], [-0.0, -5.169764690645, -0.0],
                                  [2.296178587363, -0.0, -3.767674090673]],
    'SAPT MONOMER B QUADRUPOLE': [[-2.667477720181, 1e-12, 2.247806552181], [1e-12, -3.108287027381, 0.0],
                                  [2.247806552181, 0.0, -4.390481498409]],
    'SAPT POL1,R QUADRUPOLE': [[0.150761736158, 0.0, 0.089731956532], [0.0, 0.049575584081, 0.0],
                               [0.089731956532, 0.0, 0.00644757382]],
    'SAPT POL1,R QUADRUPOLE - A': [[-0.235742624465, -0.0, -0.016419452125], [-0.0, -0.016306855628, 0.0],
                                   [-0.016419452125, 0.0, -0.026570554395]],
    'SAPT POL1,R QUADRUPOLE - B': [[0.386504360623, 0.0, 0.106151408657], [0.0, 0.065882439709, 0.0],
                                   [0.106151408657, 0.0, 0.033018128215]],
    'SAPT EXCH1(S^2) QUADRUPOLE': [[-0.469197050418, -0.0, -0.039299764774], [-0.0, -0.00353690507, -0.0],
                                   [-0.039299764774, -0.0, -0.004116247823]],
    'SAPT EXCH1(S^2) QUADRUPOLE - A': [[-0.179885276495, -0.0, -0.014078214699], [-0.0, 0.004957381503, 0.0],
                                       [-0.014078214699, 0.0, -0.004793696134]],
    'SAPT EXCH1(S^2) QUADRUPOLE - B': [[-0.289311773922, -0.0, -0.025221550076], [-0.0, -0.008494286572, -0.0],
                                       [-0.025221550076, -0.0, 0.000677448311]],
    'SAPT EXCH1 QUADRUPOLE': [[-0.47739944249, -0.0, -0.039995888199], [-0.0, -0.003671095593, -0.0],
                              [-0.039995888199, -0.0, -0.004212737374]],
    'SAPT EXCH1 QUADRUPOLE - A': [[-0.18312790913, -0.0, -0.014314069336], [-0.0, 0.005006915246, 0.0],
                                  [-0.014314069336, 0.0, -0.004885386788]],
    'SAPT EXCH1 QUADRUPOLE - B': [[-0.29427153336, -0.0, -0.025681818863], [-0.0, -0.008678010839, -0.0],
                                  [-0.025681818863, -0.0, 0.000672649414]]
}
REF_WATER_SAPTDFT_ENERGY = -0.004399083969


def _check_variables(refs, wfn):
    """Every reference variable on the wavefunction and in the globals, and the parts of
    monomers A and B adding up to their terms."""

    for key, ref in refs.items():
        assert compare_values(ref, wfn.variable(key), 8, key)
        assert compare_values(ref, psi4.variable(key), 8, f"{key} (global)")

    for key in refs:
        if key.endswith(" - A"):
            total = key[:-len(" - A")]
            parts = np.asarray(wfn.variable(f"{total} - A")) + np.asarray(wfn.variable(f"{total} - B"))
            assert compare_values(np.asarray(wfn.variable(total)), parts, 10, f"{total} = A + B")


@pytest.mark.quick
def test_sapt0_properties_he_h2():
    """SAPT0 terms of every operator kind: multipoles to the octupole, the ESP at the
    nuclei and the population charges."""

    dimer = psi4.geometry(HE_H2)
    psi4.set_options({"basis": "aug-cc-pvdz", "scf_type": "direct", **CONVERGENCE})

    properties = ["dipole", "quadrupole", "multipole(3)", "esp_at_nuclei", "mulliken_charges", "lowdin_charges"]
    energy, wfn = psi4.properties("sapt0", properties=properties, molecule=dimer, return_wfn=True)

    _check_variables(REF_HE_H2_SAPT0, wfn)

    # the monomer values are those of OEProp on the monomers in the dimer-centered basis
    for label, (real, ghost) in (("A", (1, 2)), ("B", (2, 1))):
        _, monomer = psi4.energy("scf", molecule=dimer.extract_subsets(real, ghost), return_wfn=True)
        oeprop = psi4.core.OEProp(monomer)
        oeprop.set_title("MONOMER")
        for prop in properties:
            oeprop.add(prop.upper())
        oeprop.compute()

        for name in ("DIPOLE", "QUADRUPOLE", "OCTUPOLE"):
            assert compare_values(psi4.variable(f"MONOMER {name}"), wfn.variable(f"SAPT MONOMER {label} {name}"), 8,
                                  f"monomer {label} {name} vs OEProp")
        for name in ["ESP AT CENTER 1", "ESP AT CENTER 2", "ESP AT CENTER 3", "MULLIKEN CHARGES", "LOWDIN CHARGES"]:
            assert compare_values(monomer.variable(name), wfn.variable(f"SAPT MONOMER {label} {name}"), 8,
                                  f"monomer {label} {name} vs OEProp")

    # the energy is the SAPT(DFT) code's with HF monomers
    assert compare_values(REF_HE_H2_SAPT0_ENERGY, energy, 8, "CURRENT ENERGY")
    psi4.set_options({"sapt_dft_functional": "hf", "sapt_dft_do_dhf": False})
    assert compare_values(psi4.energy("sapt(dft)", molecule=dimer), energy, 8, "energy vs energy('sapt(dft)')")


@pytest.mark.saptdft
def test_saptdft_properties_water():
    """SAPT(DFT) terms, which drop the order index of the intramonomer correlation."""

    dimer = psi4.geometry(WATER_DIMER)
    psi4.set_options({
        "basis": "cc-pvdz",
        "sapt_dft_functional": "pbe0",
        "sapt_dft_grac_shift_a": 0.19807358,
        "sapt_dft_grac_shift_b": 0.19830016,
        **CONVERGENCE
    })

    energy, wfn = psi4.properties("sapt(dft)", properties=["dipole", "quadrupole"], molecule=dimer, return_wfn=True)

    _check_variables(REF_WATER_SAPTDFT, wfn)
    assert "SAPT POL10,R DIPOLE" not in wfn.variables()

    # the energy is that of SAPT(DFT) without the delta HF correction
    assert compare_values(REF_WATER_SAPTDFT_ENERGY, energy, 8, "CURRENT ENERGY")
    psi4.set_options({"sapt_dft_do_dhf": False})
    assert compare_values(psi4.energy("sapt(dft)", molecule=dimer), energy, 8, "energy vs energy('sapt(dft)')")


@pytest.mark.quick
@pytest.mark.parametrize("prop, message", [
    ("grid_esp", "one operator per grid point"),
    ("mayer_indices", "not linear in the one-particle density"),
    ("dipole_polarizabilities", "not linear in the one-particle density"),
    ("dipol", "Did you mean"),
])
def test_sapt_properties_rejected(prop, message):
    """Properties that are not linear in the density, or not yet available."""

    psi4.geometry(HE_H2)
    psi4.set_options({"basis": "cc-pvdz"})

    with pytest.raises(psi4.ValidationError, match=message):
        psi4.properties("sapt0", properties=["dipole", prop])
