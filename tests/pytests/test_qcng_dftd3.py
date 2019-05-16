import copy

import numpy as np
import pytest
import qcelemental as qcel
from qcelemental.testing import compare, compare_recursive, compare_values, tnm

import qcengine as qcng
from qcengine.programs import dftd3
from qcengine.testing import is_program_new_enough, using_dftd3, using_dftd3_321, using_psi4, using_qcdb, using_mp2d

pytestmark = [pytest.mark.quick]

@using_dftd3
@pytest.mark.parametrize("method", [
    "b3lyp-d3",
    "b3lyp-d3m",
    "b3lyp-d3bj",
    "b3lyp-d3mbj",
])
def test_dftd3_task(method):
    json_data = {"molecule": qcng.get_molecule("eneyne"), "driver": "energy", "model": {"method": method}}

    ret = qcng.compute(json_data, "dftd3", raise_error=True, return_dict=True)

    assert ret["driver"] == "energy"
    assert "provenance" in ret
    assert "normal termination of dftd3" in ret["stdout"]

    for key in ["cpu", "hostname", "username", "wall_time"]:
        assert key in ret["provenance"]

    assert ret["success"] is True


## Resources

ref = {}
dmm = ['dimer', 'mA', 'mB', 'mAgB', 'gAmB']
ref['eneyne'] = {}
ref['eneyne']['B3LYP-D2'] = dict(zip(dmm, [-0.00390110, -0.00165271, -0.00058118, -0.00165271, -0.00058118]))
ref['eneyne']['B3LYP-D3'] = dict(zip(dmm, [-0.00285088, -0.00084340, -0.00031923, -0.00084340, -0.00031923]))
ref['eneyne']['B3LYP-D3(BJ)'] = dict(zip(dmm, [-0.00784595, -0.00394347, -0.00226683, -0.00394347, -0.00226683]))
ref['eneyne']['PBE-D2'] = dict(zip(dmm, [-0.00278650, -0.00118051, -0.00041513, -0.00118051, -0.00041513]))
ref['eneyne']['PBE-D3'] = dict(zip(dmm, [-0.00175474, -0.00045421, -0.00016839, -0.00045421, -0.00016839]))
ref['eneyne']['PBE-D3(BJ)'] = dict(zip(dmm, [-0.00475937, -0.00235265, -0.00131239, -0.00235265, -0.00131239]))
ref['eneyne']['ATM'] = dict(
    zip(dmm, [-0.000000175571, 0.000000216003, -0.000000055859, 0.000000216003, -0.000000055859]))
ref['eneyne']['MP2-DMP2'] = dict(
    zip(dmm, [0.00632174635953, 0.00265335573161, 0.00344334929607, 0.00265335573161, 0.00344334929607]))
ref['ne'] = {}
ref['ne']['B3LYP-D3(BJ)'] = {'atom': 0.0}
ref['ne']['MP2-DMP2'] = {'atom': 0.0}
ref['ne']['ATM'] = {'atom': 0.0}

gref = {}
gref['eneyne'] = {}
gref['eneyne']['B3LYP-D2'] = dict(zip(dmm, [
    np.array([
  0.00000000000000E+00,  0.48816402308826E-03, -0.52644615172697E-03,
  0.00000000000000E+00, -0.48816402308826E-03, -0.52644615172697E-03,
 -0.73597441492032E-03, -0.91579236339614E-04, -0.84500341812746E-04,
  0.73597441492032E-03, -0.91579236339614E-04, -0.84500341812746E-04,
  0.73597441492032E-03,  0.91579236339614E-04, -0.84500341812746E-04,
 -0.73597441492032E-03,  0.91579236339614E-04, -0.84500341812746E-04,
  0.00000000000000E+00,  0.00000000000000E+00, -0.49418952404353E-03,
  0.00000000000000E+00,  0.00000000000000E+00,  0.10115907534998E-02,
  0.00000000000000E+00,  0.00000000000000E+00,  0.13962586025551E-02,
  0.00000000000000E+00,  0.00000000000000E+00, -0.52276616130647E-03]).reshape((-1, 3)),
    np.array([
  0.00000000000000E+00,  0.60147644572925E-03, -0.48929640608715E-06,
  0.00000000000000E+00, -0.60147644572925E-03, -0.48929640608715E-06,
 -0.76078100143016E-03, -0.58483420364762E-04,  0.24464820304358E-06,
  0.76078100143016E-03, -0.58483420364762E-04,  0.24464820304358E-06,
  0.76078100143016E-03,  0.58483420364762E-04,  0.24464820304358E-06,
 -0.76078100143016E-03,  0.58483420364762E-04,  0.24464820304358E-06]).reshape((-1, 3)),
    np.array([
  0.00000000000000E+00,  0.00000000000000E+00, -0.56705458935397E-03,
  0.00000000000000E+00,  0.00000000000000E+00,  0.56456483332009E-03,
  0.00000000000000E+00,  0.00000000000000E+00,  0.53090524837336E-03,
  0.00000000000000E+00,  0.00000000000000E+00, -0.52841549233948E-03]).reshape((-1, 3)),
    np.array([
  0.00000000000000E+00,  0.60147644572925E-03, -0.48929640608715E-06,
  0.00000000000000E+00, -0.60147644572925E-03, -0.48929640608715E-06,
 -0.76078100143016E-03, -0.58483420364762E-04,  0.24464820304358E-06,
  0.76078100143016E-03, -0.58483420364762E-04,  0.24464820304358E-06,
  0.76078100143016E-03,  0.58483420364762E-04,  0.24464820304358E-06,
 -0.76078100143016E-03,  0.58483420364762E-04,  0.24464820304358E-06,
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00]).reshape((-1, 3)),
    np.array([
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
  0.00000000000000E+00,  0.00000000000000E+00, -0.56705458935397E-03,
  0.00000000000000E+00,  0.00000000000000E+00,  0.56456483332009E-03,
  0.00000000000000E+00,  0.00000000000000E+00,  0.53090524837336E-03,
  0.00000000000000E+00,  0.00000000000000E+00, -0.52841549233948E-03]).reshape((-1, 3)),
]))  # yapf: disable
gref['eneyne']['B3LYP-D3'] = dict(zip(dmm, [
    np.array([
  0.67762635780344E-20,  0.19657186672293E-03, -0.23180716200687E-03,
  0.50821976835258E-20, -0.19657186672293E-03, -0.23180716200687E-03,
 -0.83754349667195E-04,  0.45844828386013E-04, -0.92969637976992E-04,
  0.83754349667195E-04,  0.45844828386013E-04, -0.92969637976992E-04,
  0.83754349667195E-04, -0.45844828386013E-04, -0.92969637976992E-04,
 -0.83754349667195E-04, -0.45844828386013E-04, -0.92969637976992E-04,
  0.00000000000000E+00,  0.00000000000000E+00, -0.11800508571549E-03,
 -0.33881317890172E-20, -0.50821976835258E-20,  0.62302832736499E-03,
  0.50821976835258E-20,  0.33881317890172E-20,  0.50037535445493E-03,
  0.00000000000000E+00, -0.52939559203394E-22, -0.16990572018272E-03]).reshape((-1, 3)),
    np.array([
  0.20328790734103E-19,  0.24171499732116E-03, -0.20480842481032E-06,
 -0.16940658945086E-20, -0.24171499732116E-03, -0.20480842481032E-06,
 -0.10776189540054E-03,  0.78926689997812E-04,  0.10240421240516E-06,
  0.10776189540054E-03,  0.78926689997812E-04,  0.10240421240516E-06,
  0.10776189540054E-03, -0.78926689997812E-04,  0.10240421240516E-06,
 -0.10776189540054E-03, -0.78926689997812E-04,  0.10240421240516E-06]).reshape((-1, 3)),
    np.array([
  0.00000000000000E+00,  0.00000000000000E+00, -0.21752286612122E-03,
  0.00000000000000E+00,  0.00000000000000E+00,  0.21634915516554E-03,
  0.00000000000000E+00,  0.00000000000000E+00,  0.17823532330490E-03,
  0.00000000000000E+00,  0.00000000000000E+00, -0.17706161234922E-03]).reshape((-1, 3)),
    np.array([
  0.20328790734103E-19,  0.24171499732116E-03, -0.20480842481032E-06,
 -0.16940658945086E-20, -0.24171499732116E-03, -0.20480842481032E-06,
 -0.10776189540054E-03,  0.78926689997812E-04,  0.10240421240516E-06,
  0.10776189540054E-03,  0.78926689997812E-04,  0.10240421240516E-06,
  0.10776189540054E-03, -0.78926689997812E-04,  0.10240421240516E-06,
 -0.10776189540054E-03, -0.78926689997812E-04,  0.10240421240516E-06,
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00]).reshape((-1, 3)),
    np.array([
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
  0.00000000000000E+00,  0.00000000000000E+00, -0.21752286612122E-03,
  0.00000000000000E+00,  0.00000000000000E+00,  0.21634915516554E-03,
  0.00000000000000E+00,  0.00000000000000E+00,  0.17823532330490E-03,
  0.00000000000000E+00,  0.00000000000000E+00, -0.17706161234922E-03]).reshape((-1, 3)),
]))  # yapf: disable
gref['eneyne']['B3LYP-D3(BJ)'] = dict(zip(dmm, [
    np.array([
  0.16940658945086E-20, -0.10896372137622E-03, -0.28496931787936E-03,
  0.33881317890172E-20,  0.10896372137622E-03, -0.28496931787936E-03,
  0.56547183189867E-04, -0.10791733716132E-03, -0.81750328898176E-04,
 -0.56547183189867E-04, -0.10791733716132E-03, -0.81750328898176E-04,
 -0.56547183189867E-04,  0.10791733716132E-03, -0.81750328898176E-04,
  0.56547183189867E-04,  0.10791733716132E-03, -0.81750328898176E-04,
  0.00000000000000E+00,  0.00000000000000E+00,  0.14698301085008E-03,
 -0.16940658945086E-20, -0.33881317890172E-20,  0.43655907696000E-03,
  0.00000000000000E+00,  0.33881317890172E-20,  0.23688438591518E-03,
  0.00000000000000E+00, -0.52939559203394E-22,  0.76513477626168E-04]).reshape((-1, 3)),
    np.array([
 -0.33881317890172E-20, -0.54157860939394E-04,  0.11299781801723E-07,
 -0.93173624197973E-20,  0.54157860939394E-04,  0.11299781801723E-07,
  0.35880725530239E-04, -0.79323052619042E-04, -0.56498909008614E-08,
 -0.35880725530239E-04, -0.79323052619042E-04, -0.56498909008614E-08,
 -0.35880725530239E-04,  0.79323052619042E-04, -0.56498909008614E-08,
  0.35880725530239E-04,  0.79323052619042E-04, -0.56498909008614E-08]).reshape((-1, 3)),
    np.array([
  0.00000000000000E+00,  0.00000000000000E+00,  0.45552310986933E-04,
  0.00000000000000E+00,  0.00000000000000E+00, -0.45561218665227E-04,
  0.00000000000000E+00,  0.00000000000000E+00, -0.69342175541743E-04,
  0.00000000000000E+00,  0.00000000000000E+00,  0.69351083220036E-04]).reshape((-1, 3)),
    np.array([
 -0.33881317890172E-20, -0.54157860939394E-04,  0.11299781801723E-07,
 -0.93173624197973E-20,  0.54157860939394E-04,  0.11299781801723E-07,
  0.35880725530239E-04, -0.79323052619042E-04, -0.56498909008614E-08,
 -0.35880725530239E-04, -0.79323052619042E-04, -0.56498909008614E-08,
 -0.35880725530239E-04,  0.79323052619042E-04, -0.56498909008614E-08,
  0.35880725530239E-04,  0.79323052619042E-04, -0.56498909008614E-08,
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00]).reshape((-1, 3)),
    np.array([
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
  0.00000000000000E+00,  0.00000000000000E+00,  0.45552310986933E-04,
  0.00000000000000E+00,  0.00000000000000E+00, -0.45561218665227E-04,
  0.00000000000000E+00,  0.00000000000000E+00, -0.69342175541743E-04,
  0.00000000000000E+00,  0.00000000000000E+00,  0.69351083220036E-04]).reshape((-1, 3)),
]))  # yapf: disable
gref['eneyne']['PBE-D2'] = dict(zip(dmm, [
    np.array([
  0.00000000000000E+00,  0.34868860375520E-03, -0.37603298259607E-03,
  0.00000000000000E+00, -0.34868860375520E-03, -0.37603298259607E-03,
 -0.52569603453084E-03, -0.65413743213220E-04, -0.60357389750118E-04,
  0.52569603453084E-03, -0.65413743213220E-04, -0.60357389750118E-04,
  0.52569603453084E-03,  0.65413743213220E-04, -0.60357389750118E-04,
 -0.52569603453084E-03,  0.65413743213220E-04, -0.60357389750118E-04,
  0.00000000000000E+00,  0.00000000000000E+00, -0.35299253320442E-03,
  0.00000000000000E+00,  0.00000000000000E+00,  0.72256485674234E-03,
  0.00000000000000E+00,  0.00000000000000E+00,  0.99732761854534E-03,
  0.00000000000000E+00,  0.00000000000000E+00, -0.37340441789063E-03]).reshape((-1, 3)),
    np.array([
  0.00000000000000E+00,  0.42962605217439E-03, -0.34949744879114E-06,
  0.00000000000000E+00, -0.42962605217439E-03, -0.34949744879114E-06,
 -0.54341502569968E-03, -0.41773873586195E-04,  0.17474872439557E-06,
  0.54341502569968E-03, -0.41773873586195E-04,  0.17474872439557E-06,
  0.54341502569968E-03,  0.41773873586195E-04,  0.17474872439557E-06,
 -0.54341502569968E-03,  0.41773873586195E-04,  0.17474872439557E-06]).reshape((-1, 3)),
    np.array([
  0.00000000000000E+00,  0.00000000000000E+00, -0.40503901078976E-03,
  0.00000000000000E+00,  0.00000000000000E+00,  0.40326061354193E-03,
  0.00000000000000E+00,  0.00000000000000E+00,  0.37921805177386E-03,
  0.00000000000000E+00,  0.00000000000000E+00, -0.37743965452603E-03]).reshape((-1, 3)),
    np.array([
  0.00000000000000E+00,  0.42962605217439E-03, -0.34949744879114E-06,
  0.00000000000000E+00, -0.42962605217439E-03, -0.34949744879114E-06,
 -0.54341502569968E-03, -0.41773873586195E-04,  0.17474872439557E-06,
  0.54341502569968E-03, -0.41773873586195E-04,  0.17474872439557E-06,
  0.54341502569968E-03,  0.41773873586195E-04,  0.17474872439557E-06,
 -0.54341502569968E-03,  0.41773873586195E-04,  0.17474872439557E-06,
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00]).reshape((-1, 3)),
    np.array([
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
  0.00000000000000E+00,  0.00000000000000E+00, -0.40503901078976E-03,
  0.00000000000000E+00,  0.00000000000000E+00,  0.40326061354193E-03,
  0.00000000000000E+00,  0.00000000000000E+00,  0.37921805177386E-03,
  0.00000000000000E+00,  0.00000000000000E+00, -0.37743965452603E-03]).reshape((-1, 3)),
]))  # yapf: disable
gref['eneyne']['PBE-D3'] = dict(zip(dmm, [
    np.array([
  0.33881317890172E-20,  0.97730853016389E-04, -0.71901324069440E-04,
  0.29646153153901E-20, -0.97730853016389E-04, -0.71901324069440E-04,
 -0.31222554291636E-04,  0.29545643062003E-04, -0.67132324795951E-04,
  0.31222554291636E-04,  0.29545643062003E-04, -0.67132324795951E-04,
  0.31222554291636E-04, -0.29545643062003E-04, -0.67132324795951E-04,
 -0.31222554291636E-04, -0.29545643062003E-04, -0.67132324795951E-04,
  0.00000000000000E+00,  0.00000000000000E+00, -0.20867204655394E-04,
  0.00000000000000E+00,  0.00000000000000E+00,  0.30651555323142E-03,
  0.16940658945086E-20,  0.16940658945086E-20,  0.18769576731255E-03,
  0.00000000000000E+00,  0.00000000000000E+00, -0.61012168565887E-04]).reshape((-1, 3)),
    np.array([
  0.33881317890172E-20,  0.11021182403760E-03, -0.93982803767757E-07,
 -0.42351647362715E-21, -0.11021182403760E-03, -0.93982803767758E-07,
 -0.48220259417857E-04,  0.52933097691669E-04,  0.46991401883879E-07,
  0.48220259417857E-04,  0.52933097691669E-04,  0.46991401883879E-07,
  0.48220259417857E-04, -0.52933097691669E-04,  0.46991401883879E-07,
 -0.48220259417857E-04, -0.52933097691669E-04,  0.46991401883879E-07]).reshape((-1, 3)),
    np.array([
  0.00000000000000E+00,  0.00000000000000E+00, -0.99901515312115E-04,
  0.00000000000000E+00,  0.00000000000000E+00,  0.99340886211351E-04,
  0.00000000000000E+00,  0.00000000000000E+00,  0.67878928346081E-04,
  0.00000000000000E+00,  0.00000000000000E+00, -0.67318299245317E-04]).reshape((-1, 3)),
    np.array([
  0.33881317890172E-20,  0.11021182403760E-03, -0.93982803767757E-07,
 -0.42351647362715E-21, -0.11021182403760E-03, -0.93982803767758E-07,
 -0.48220259417857E-04,  0.52933097691669E-04,  0.46991401883879E-07,
  0.48220259417857E-04,  0.52933097691669E-04,  0.46991401883879E-07,
  0.48220259417857E-04, -0.52933097691669E-04,  0.46991401883879E-07,
 -0.48220259417857E-04, -0.52933097691669E-04,  0.46991401883879E-07,
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00]).reshape((-1, 3)),
    np.array([
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
  0.00000000000000E+00,  0.00000000000000E+00, -0.99901515312115E-04,
  0.00000000000000E+00,  0.00000000000000E+00,  0.99340886211351E-04,
  0.00000000000000E+00,  0.00000000000000E+00,  0.67878928346081E-04,
  0.00000000000000E+00,  0.00000000000000E+00, -0.67318299245317E-04]).reshape((-1, 3)),
]))  # yapf: disable
gref['eneyne']['PBE-D3(BJ)'] = dict(zip(dmm, [
    np.array([
  0.00000000000000E+00, -0.61939589939064E-04, -0.16066534797355E-03,
  0.25410988417629E-20,  0.61939589939064E-04, -0.16066534797355E-03,
  0.35330272921363E-04, -0.65816270722009E-04, -0.53748175167354E-04,
 -0.35330272921363E-04, -0.65816270722009E-04, -0.53748175167354E-04,
 -0.35330272921363E-04,  0.65816270722009E-04, -0.53748175167354E-04,
  0.35330272921363E-04,  0.65816270722009E-04, -0.53748175167354E-04,
  0.00000000000000E+00,  0.00000000000000E+00,  0.98332430764322E-04,
 -0.25410988417629E-20, -0.33881317890172E-20,  0.25661736016373E-03,
  0.00000000000000E+00,  0.00000000000000E+00,  0.13371752089002E-03,
  0.00000000000000E+00, -0.26469779601697E-22,  0.47656084798449E-04]).reshape((-1, 3)),
    np.array([
 -0.25410988417629E-20, -0.31329250082804E-04,  0.61088639781542E-08,
 -0.50821976835258E-20,  0.31329250082804E-04,  0.61088639781542E-08,
  0.21959764459240E-04, -0.47293026603847E-04, -0.30544319890771E-08,
 -0.21959764459240E-04, -0.47293026603847E-04, -0.30544319890771E-08,
 -0.21959764459240E-04,  0.47293026603847E-04, -0.30544319890771E-08,
  0.21959764459240E-04,  0.47293026603847E-04, -0.30544319890771E-08]).reshape((-1, 3)),
    np.array([
  0.00000000000000E+00,  0.00000000000000E+00,  0.25685884880777E-04,
  0.00000000000000E+00,  0.00000000000000E+00, -0.25704336611069E-04,
  0.00000000000000E+00,  0.00000000000000E+00, -0.41528315631943E-04,
  0.00000000000000E+00,  0.00000000000000E+00,  0.41546767362234E-04]).reshape((-1, 3)),
    np.array([
 -0.25410988417629E-20, -0.31329250082804E-04,  0.61088639781542E-08,
 -0.50821976835258E-20,  0.31329250082804E-04,  0.61088639781542E-08,
  0.21959764459240E-04, -0.47293026603847E-04, -0.30544319890771E-08,
 -0.21959764459240E-04, -0.47293026603847E-04, -0.30544319890771E-08,
 -0.21959764459240E-04,  0.47293026603847E-04, -0.30544319890771E-08,
  0.21959764459240E-04,  0.47293026603847E-04, -0.30544319890771E-08,
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00]).reshape((-1, 3)),
    np.array([
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
  0.00000000000000E+00,  0.00000000000000E+00,  0.25685884880777E-04,
  0.00000000000000E+00,  0.00000000000000E+00, -0.25704336611069E-04,
  0.00000000000000E+00,  0.00000000000000E+00, -0.41528315631943E-04,
  0.00000000000000E+00,  0.00000000000000E+00,  0.41546767362234E-04]).reshape((-1, 3)),
]))  # yapf: disable
gref['eneyne']['ATM'] = dict(zip(dmm, [
    np.array([
  0.00000000000000E+00, -0.57988139838201E-06, -0.71628554331971E-06,
  0.00000000000000E+00,  0.57988139838201E-06, -0.71628554331971E-06,
  0.53149296386534E-06, -0.41638019417978E-06,  0.52694338024860E-06,
 -0.53149296386534E-06, -0.41638019417978E-06,  0.52694338024860E-06,
 -0.53149296386533E-06,  0.41638019417978E-06,  0.52694338024860E-06,
  0.53149296386533E-06,  0.41638019417978E-06,  0.52694338024858E-06,
  0.00000000000000E+00,  0.00000000000000E+00, -0.92557313363084E-06,
  0.00000000000000E+00,  0.00000000000000E+00,  0.31010265235900E-06,
  0.00000000000000E+00,  0.00000000000000E+00,  0.10194777599160E-05,
  0.00000000000000E+00,  0.00000000000000E+00, -0.10792097129990E-05]).reshape((-1, 3)),
    np.array([
  0.00000000000000E+00, -0.45154573778694E-07,  0.11131133827146E-09,
  0.00000000000000E+00,  0.45154573778707E-07,  0.11131133827146E-09,
  0.10133274017225E-06, -0.72175367263952E-07, -0.55655669135736E-10,
 -0.10133274017225E-06, -0.72175367263966E-07, -0.55655669135736E-10,
 -0.10133274017227E-06,  0.72175367263952E-07, -0.55655669135736E-10,
  0.10133274017227E-06,  0.72175367263966E-07, -0.55655669135736E-10]).reshape((-1, 3)),
    np.array([
  0.00000000000000E+00,  0.00000000000000E+00,  0.14079365564105E-07,
  0.00000000000000E+00,  0.00000000000000E+00, -0.14067311316192E-07,
  0.00000000000000E+00,  0.00000000000000E+00,  0.67277034390041E-07,
  0.00000000000000E+00,  0.00000000000000E+00, -0.67289088637954E-07]).reshape((-1, 3)),
    np.array([
  0.00000000000000E+00, -0.45154573778694E-07,  0.11131133827146E-09,
  0.00000000000000E+00,  0.45154573778707E-07,  0.11131133827146E-09,
  0.10133274017225E-06, -0.72175367263952E-07, -0.55655669135736E-10,
 -0.10133274017225E-06, -0.72175367263966E-07, -0.55655669135736E-10,
 -0.10133274017227E-06,  0.72175367263952E-07, -0.55655669135736E-10,
  0.10133274017227E-06,  0.72175367263966E-07, -0.55655669135736E-10,
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00]).reshape((-1, 3)),
    np.array([
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
  0.00000000000000E+00,  0.00000000000000E+00,  0.00000000000000E+00,
  0.00000000000000E+00,  0.00000000000000E+00,  0.14079365564105E-07,
  0.00000000000000E+00,  0.00000000000000E+00, -0.14067311316192E-07,
  0.00000000000000E+00,  0.00000000000000E+00,  0.67277034390041E-07,
  0.00000000000000E+00,  0.00000000000000E+00, -0.67289088637954E-07]).reshape((-1, 3)),
]))  # yapf: disable
gref['ne'] = {}
gref['ne']['B3LYP-D3(BJ)'] = {'atom': np.zeros(3).reshape((-1, 3))}
gref['ne']['MP2-DMP2'] = {'atom': np.zeros(3).reshape((-1, 3))}
gref['ne']['ATM'] = {'atom': np.zeros(3).reshape((-1, 3))}

seneyne = """
C   0.000000  -0.667578  -2.124659
C   0.000000   0.667578  -2.124659
H   0.923621  -1.232253  -2.126185
H  -0.923621  -1.232253  -2.126185
H  -0.923621   1.232253  -2.126185
H   0.923621   1.232253  -2.126185
--
C   0.000000   0.000000   2.900503
C   0.000000   0.000000   1.693240
H   0.000000   0.000000   0.627352
H   0.000000   0.000000   3.963929
"""

sne = """
Ne 0 0 0
"""


def eneyne_ne_qcdbmols():
    if not is_program_new_enough("psi4", "1.3rc2"):
        pytest.skip("Psi4 requires at least Psi4 v1.3rc2")
    from psi4.driver import qcdb

    eneyne = qcdb.Molecule(seneyne)
    ne = qcdb.Molecule(sne)
    mols = {
        'eneyne': {
            'dimer': eneyne,
            'mA': eneyne.extract_subsets(1),
            'mB': eneyne.extract_subsets(2),
            'mAgB': eneyne.extract_subsets(1, 2),
            'gAmB': eneyne.extract_subsets(2, 1),
        },
        'ne': {
            'atom': ne,
        }
    }
    return mols


def eneyne_ne_psi4mols():
    if not is_program_new_enough("psi4", "1.3rc2"):
        pytest.skip("Psi4 requires at least Psi4 v1.3rc2")
    import psi4

    eneyne = psi4.core.Molecule.from_string(seneyne)
    ne = psi4.core.Molecule.from_string(sne)
    mols = {
        'eneyne': {
            'dimer': eneyne,
            'mA': eneyne.extract_subsets(1),
            'mB': eneyne.extract_subsets(2),
            'mAgB': eneyne.extract_subsets(1, 2),
            'gAmB': eneyne.extract_subsets(2, 1),
        },
        'ne': {
            'atom': ne,
        }
    }
    return mols


def eneyne_ne_qcschemamols():

    eneyne = qcel.molparse.to_schema(qcel.molparse.from_string(seneyne)['qm'], dtype=2)
    mA = qcel.molparse.to_schema(qcel.molparse.from_string('\n'.join(seneyne.splitlines()[:7]))['qm'], dtype=2)
    mB = qcel.molparse.to_schema(qcel.molparse.from_string('\n'.join(seneyne.splitlines()[-4:]))['qm'], dtype=2)
    ne = qcel.molparse.to_schema(qcel.molparse.from_string(sne)['qm'], dtype=2)

    mAgB = qcel.molparse.from_string(seneyne)['qm']
    mAgB['real'] = [(iat < mAgB['fragment_separators'][0])
                    for iat in range(len(mAgB['elem']))]  # works b/c chgmult doesn't need refiguring
    mAgB = qcel.molparse.to_schema(mAgB, dtype=2)

    gAmB = qcel.molparse.from_string(seneyne)['qm']
    gAmB['real'] = [(iat >= gAmB['fragment_separators'][0]) for iat in range(len(gAmB['elem']))]
    gAmB = qcel.molparse.to_schema(gAmB, dtype=2)

    mols = {
        'eneyne': {
            'dimer': eneyne,
            'mA': mA,
            'mB': mB,
            'mAgB': mAgB,
            'gAmB': gAmB,
        },
        'ne': {
            'atom': ne,
        }
    }
    return mols


db3lypd3bj = {
    'dashlevel': 'd3bj',
    'dashparams': {
        's8': 1.9889,
        's6': 1.0,
        'a2': 4.4211,
        'a1': 0.3981
    },
    'dashparams_citation': '',
    'fctldash': 'b3lyp-d3(bj)'
}
db3lypd3bjcustom = copy.deepcopy(db3lypd3bj)
db3lypd3bjcustom['fctldash'] = ''
db3lypd3bjcustom['dashparams']['a2'] = 5.4211

dpbed3zero = {
    'dashlevel': 'd3zero',
    'dashparams': {
        's6': 1.0,
        's8': 0.722,
        'sr6': 1.217,
        'sr8': 1.0,
        'alpha6': 14.0
    },
    'dashparams_citation': '',
    'fctldash': 'pbe-d3'
}

atmgr = {
    'dashlevel': 'atmgr',
    'dashparams': {
        'alpha6': 14.0,
    },
    'dashparams_citation': '',
    'fctldash': 'atm(gr)',
}

chg = {
    'dashlevel': 'chg',
    'dashparams': {
        's6': 1.0,
    },
    'dashparams_citation': '',
    'fctldash': 'chg',
}

dmp2dmp2 = {
    'dashlevel': 'dmp2',
    'dashparams': {
        's8': 1.187,
        'a1': 0.944,
        'a2': 0.480,
        'rcut': 0.72,
        'w': 0.20,
    },
    'dashparams_citation': '',
    'fctldash': 'mp2-dmp2'
}


def _compute_key(pjrec):
    return pjrec['fctldash'].upper()


## Tests


@pytest.mark.parametrize("inp,expected", [
    (({'name_hint': 'b3lyp', 'level_hint': 'd3bj'}, 'B3LYP-D3(BJ)'), db3lypd3bj),
    (({'name_hint': 'b3LYP', 'level_hint': 'D3bj'}, 'B3LYP-D3(BJ)'), db3lypd3bj),
    (({'param_tweaks': {'s8': 1.9889, 's6': 1.0, 'a2': 4.4211, 'a1': 0.3981}, 'level_hint': 'd3bj'}, 'B3LYP-D3(BJ)'), db3lypd3bj),
    (({'name_hint': 'b3lyp', 'level_hint': 'd3bJ', 'param_tweaks': {'a2': 4.4211}}, 'B3LYP-D3(BJ)'), db3lypd3bj),
    (({'verbose': 3, 'name_hint': 'b3lyp', 'level_hint': 'd3bJ', 'param_tweaks': {'a2': 5.4211}}, ''), db3lypd3bjcustom),
    (({'name_hint': 'b3lyp-d3bj', 'param_tweaks': {'a2': 4.4211}}, 'B3LYP-D3(BJ)'), db3lypd3bj),
    (({'name_hint': 'pbe', 'level_hint': 'd3zero'}, 'PBE-D3'), dpbed3zero),
    (({'name_hint': 'pbe', 'level_hint': 'd3'}, 'PBE-D3'), dpbed3zero),
    (({'name_hint': 'pbe-d3'}, 'PBE-D3'), dpbed3zero),
    (({'name_hint': 'atm(gr)', 'level_hint': 'atmgr'}, 'ATM(GR)'), atmgr),
    (({'name_hint': 'atmgr'}, 'ATM(GR)'), atmgr),
    (({'name_hint': 'bp86-atmgr'}, 'ATM(GR)'), atmgr),
    (({'name_hint': 'asdf-chg'}, 'CHG'), chg),
    (({'name_hint': 'mp2-dmp2'}, 'MP2-DMP2'), dmp2dmp2),
    (({'name_hint': 'MP2', 'level_hint': 'dmp2'}, 'MP2-DMP2'), dmp2dmp2),
])  # yapf: disable
def test_dftd3__from_arrays(inp, expected):
    res = dftd3.from_arrays(**inp[0])
    assert compare_recursive(expected, res, atol=1.e-4)
    assert compare(inp[1], _compute_key(res), 'key')
    res = dftd3.from_arrays(name_hint=res['fctldash'], level_hint=res['dashlevel'], param_tweaks=res['dashparams'])
    assert compare_recursive(expected, res, tnm() + ' idempotent', atol=1.e-4)


@pytest.mark.parametrize("inp", [
    ({'name_hint': 'b3lyp', 'level_hint': 'd3bJ', 'param_tweaks': {'a3': 5.4211}}),
    ({'name_hint': 'fakeb3lyp', 'level_hint': 'd3bJ', 'param_tweaks': {'s6': 5.4211}}),
    ({'level_hint': 'd3bJ', 'param_tweaks': {'s6': 5.4211}}),
    ({'name_hint': 'b3lyp-d3bj', 'param_tweaks': {'a2': 4.4211, 'zzz': 0.0}}),
    ({'name_hint': 'asdf-d4'}),
    ({'name_hint': 'atm(gr)', 'level_hint': 'chg'}),
])  # yapf:disable
def test_dftd3__from_arrays__error(inp):
    with pytest.raises(ValueError):
        dftd3.from_arrays(**inp)


def test_dftd3__from_arrays__supplement():
    ans = {
        'dashlevel': 'chg',
        'dashparams': {
            's6': 4.05
        },
        'fctldash': 'asdf-d4',
        'dashparams_citation': '    mypaper\n'
    }
    supp = {'chg': {'definitions': {'asdf-d4': {'params': {'s6': 4.05}, 'citation': '    mypaper\n'}}}}

    res = dftd3.from_arrays(name_hint='asdf-d4', level_hint='chg', dashcoeff_supplement=supp)
    assert compare_recursive(ans, res, atol=1.e-4)
    with pytest.raises(ValueError) as e:
        dftd3.from_arrays(name_hint=res['fctldash'], level_hint=res['dashlevel'], param_tweaks=res['dashparams'])
    assert "Can't guess -D correction level" in str(e)
    res = dftd3.from_arrays(
        name_hint=res['fctldash'],
        level_hint=res['dashlevel'],
        param_tweaks=res['dashparams'],
        dashcoeff_supplement=supp)
    assert compare_recursive(ans, res, tnm() + ' idempotent', atol=1.e-4)


@using_dftd3
def test_3():
    sys = qcel.molparse.from_string(seneyne)['qm']

    resinp = {
        'schema_name': 'qcschema_input',
        'schema_version': 1,
        'molecule': qcel.molparse.to_schema(sys, dtype=2),
        'driver': 'energy',
        'model': {
            'method': 'b3lyp',
        },
        'keywords': {
            'level_hint': 'd3bj'
        },
    }
    res = qcng.compute(resinp, 'dftd3', raise_error=True)
    res = res.dict()

    #res = dftd3.run_dftd3_from_arrays(molrec=sys, name_hint='b3lyp', level_hint='d3bj')
    assert compare('B3LYP-D3(BJ)', _compute_key(res['extras']['info']), 'key')


@using_dftd3
@pytest.mark.parametrize(
    "subjects",
    [
        pytest.param(eneyne_ne_psi4mols, marks=using_psi4),
        pytest.param(eneyne_ne_qcdbmols,
                     marks=using_psi4),  # needs qcdb.Molecule, presently more common in psi4 than in qcdb
    ],
    ids=['qmol', 'pmol'])
@pytest.mark.parametrize(
    "inp", [
        ({'first': 'b3lyp', 'second': 'd', 'parent': 'eneyne', 'subject': 'dimer', 'lbl': 'B3LYP-D2'}),
        ({'first': 'b3lyp', 'second': 'd3bj', 'parent': 'eneyne', 'subject': 'mA', 'lbl': 'B3LYP-D3(BJ)'}),
        ({'first': 'pbe', 'second': 'd3zero', 'parent': 'eneyne', 'subject': 'mB', 'lbl': 'PBE-D3'}),
        ({'first': 'pbe', 'second': 'd3zero', 'parent': 'eneyne', 'subject': 'gAmB', 'lbl': 'PBE-D3'}),
        ({'first': 'pbe', 'second': 'd2', 'parent': 'eneyne', 'subject': 'mAgB', 'lbl': 'PBE-D2'}),
        ({'first': 'b3lyp', 'second': 'd3bj', 'parent': 'ne', 'subject': 'atom', 'lbl': 'B3LYP-D3(BJ)'}),
        #({'first': '', 'second': 'atmgr', 'parent': 'eneyne', 'subject': 'dimer', 'lbl': 'ATM'}),
        #({'first': 'b3lyp', 'second': 'atmgr', 'parent': 'eneyne', 'subject': 'mA', 'lbl': 'ATM'}),
        #({'first': 'pbe', 'second': 'atm(gr)', 'parent': 'eneyne', 'subject': 'mB', 'lbl': 'ATM'}),
        #({'first': '', 'second': 'ATMgr', 'parent': 'eneyne', 'subject': 'mAgB', 'lbl': 'ATM'}),
        # below two xfail until dftd3 that's only 2-body is out of psi4 proper
        pytest.param({'first': 'atmgr', 'second': 'atmgr', 'parent': 'eneyne', 'subject': 'gAmB', 'lbl': 'ATM'}, marks=[using_dftd3_321, pytest.mark.xfail]),
        pytest.param({'first': 'pbe-atmgr', 'second': None, 'parent': 'ne', 'subject': 'atom', 'lbl': 'ATM'}, marks=[using_dftd3_321, pytest.mark.xfail]),
    ])  # yapf: disable
def test_molecule__run_dftd3__23body(inp, subjects):
    subject = subjects()[inp['parent']][inp['subject']]
    expected = ref[inp['parent']][inp['lbl']][inp['subject']]
    gexpected = gref[inp['parent']][inp['lbl']][inp['subject']]

    E, G = subject.run_dftd3(inp['first'], inp['second'])
    assert compare_values(expected, E, atol=1.e-7)
    assert compare_values(gexpected, G, atol=1.e-7)


@using_qcdb
def test_qcdb__energy_d3():
    eneyne = qcdb.set_molecule(seneyne)
    eneyne.update_geometry()

    E, jrec = qcdb.energy('d3-b3lyp-d2', return_wfn=True)
    assert compare_values(ref['eneyne']['B3LYP-D2']['dimer'], E, 7, 'P: Ethene-Ethyne -D2')
    assert compare_values(ref['eneyne']['B3LYP-D2']['dimer'], jrec['qcvars']['DISPERSION CORRECTION ENERGY'].data, 7,
                          tnm())
    assert compare_values(ref['eneyne']['B3LYP-D2']['dimer'],
                          jrec['qcvars']['B3LYP-D2 DISPERSION CORRECTION ENERGY'].data, 7, tnm())

    mA = eneyne.extract_subsets(1)

    E, jrec = qcdb.energy('d3-b3lyp-d3bj', return_wfn=True, molecule=mA)
    assert compare_values(ref['eneyne']['B3LYP-D3(BJ)']['mA'], E, 7, tnm())
    assert compare_values(ref['eneyne']['B3LYP-D3(BJ)']['mA'], jrec['qcvars']['DISPERSION CORRECTION ENERGY'].data, 7,
                          tnm())
    assert compare_values(ref['eneyne']['B3LYP-D3(BJ)']['mA'],
                          jrec['qcvars']['B3LYP-D3(BJ) DISPERSION CORRECTION ENERGY'].data, 7, tnm())


@using_mp2d
@pytest.mark.parametrize(
    "subjects",
    [
        pytest.param(eneyne_ne_psi4mols, marks=using_psi4),
        pytest.param(eneyne_ne_qcdbmols,
                     marks=using_psi4),  # needs qcdb.Molecule, presently more common in psi4 than in qcdb
        pytest.param(eneyne_ne_qcschemamols),
    ],
    ids=['qmol', 'pmol', 'qcmol'])
@pytest.mark.parametrize("inp", [
    ({'parent': 'eneyne', 'name': 'mp2d-mp2-dmp2', 'subject': 'dimer', 'lbl': 'MP2-DMP2'}),
    ({'parent': 'eneyne', 'name': 'mp2d-mp2-dmp2', 'subject': 'mA', 'lbl': 'MP2-DMP2'}),
    ({'parent': 'eneyne', 'name': 'mp2d-mp2-dmp2', 'subject': 'mB', 'lbl': 'MP2-DMP2'}),
    ({'parent': 'eneyne', 'name': 'mp2d-mp2-dmp2', 'subject': 'gAmB', 'lbl': 'MP2-DMP2'}),
    ({'parent': 'eneyne', 'name': 'mp2d-mp2-dmp2', 'subject': 'mAgB', 'lbl': 'MP2-DMP2'}),
    ({'parent': 'ne', 'name': 'mp2d-mp2-dmp2', 'subject': 'atom', 'lbl': 'MP2-DMP2'}),
])  # yapf: disable
def test_mp2d__run_mp2d__2body(inp, subjects, request):
    subject = subjects()[inp['parent']][inp['subject']]
    expected = ref[inp['parent']][inp['lbl']][inp['subject']]
    #gexpected = gref[inp['parent']][inp['lbl']][inp['subject']].ravel()

    if 'qcmol' in request.node.name:
        mol = subject
    else:
        mol = subject.to_schema(dtype=2)

    resinp = {
        'schema_name': 'qcschema_input',
        'schema_version': 1,
        'molecule': mol,
        'driver': 'energy', #gradient',
        'model': {
            'method': inp['name']
        },
        'keywords': {},
    }
    jrec = qcng.compute(resinp, 'mp2d', raise_error=True)
    jrec = jrec.dict()

    #assert len(jrec['extras']['qcvars']) == 8

    assert compare_values(expected, jrec['extras']['qcvars']['CURRENT ENERGY'], atol=1.e-7)
    assert compare_values(expected, jrec['extras']['qcvars']['DISPERSION CORRECTION ENERGY'], atol=1.e-7)
    assert compare_values(expected, jrec['extras']['qcvars'][inp['lbl'] + ' DISPERSION CORRECTION ENERGY'], atol=1.e-7)


@using_dftd3
@pytest.mark.parametrize(
    "subjects",
    [
        pytest.param(eneyne_ne_psi4mols, marks=using_psi4),
        pytest.param(eneyne_ne_qcdbmols,
                     marks=using_psi4),  # needs qcdb.Molecule, presently more common in psi4 than in qcdb
        pytest.param(eneyne_ne_qcschemamols),
    ],
    ids=['qmol', 'pmol', 'qcmol'])
@pytest.mark.parametrize("inp", [
    ({'parent': 'eneyne', 'name': 'd3-b3lyp-d', 'subject': 'dimer', 'lbl': 'B3LYP-D2'}),
    ({'parent': 'eneyne', 'name': 'd3-b3lyp-d3bj', 'subject': 'mA', 'lbl': 'B3LYP-D3(BJ)'}),
    ({'parent': 'eneyne', 'name': 'd3-PBE-D3zero', 'subject': 'mB', 'lbl': 'PBE-D3'}),
    ({'parent': 'eneyne', 'name': 'd3-PBE-D3zero', 'subject': 'gAmB', 'lbl': 'PBE-D3'}),
    ({'parent': 'eneyne', 'name': 'd3-PBE-D2', 'subject': 'mAgB', 'lbl': 'PBE-D2'}),
    ({'parent': 'ne', 'name': 'd3-b3lyp-d3bj', 'subject': 'atom', 'lbl': 'B3LYP-D3(BJ)'}),
])  # yapf: disable
def test_dftd3__run_dftd3__2body(inp, subjects, request):
    subject = subjects()[inp['parent']][inp['subject']]
    expected = ref[inp['parent']][inp['lbl']][inp['subject']]
    gexpected = gref[inp['parent']][inp['lbl']][inp['subject']].ravel()

    if 'qcmol' in request.node.name:
        mol = subject
    else:
        mol = subject.to_schema(dtype=2)

    resinp = {
        'schema_name': 'qcschema_input',
        'schema_version': 1,
        'molecule': mol,
        'driver': 'gradient',
        'model': {
            'method': inp['name']
        },
        'keywords': {},
    }
    jrec = qcng.compute(resinp, 'dftd3', raise_error=True)
    jrec = jrec.dict()

    assert len(jrec['extras']['qcvars']) == 8

    assert compare_values(expected, jrec['extras']['qcvars']['CURRENT ENERGY'], atol=1.e-7)
    assert compare_values(expected, jrec['extras']['qcvars']['DISPERSION CORRECTION ENERGY'], atol=1.e-7)
    assert compare_values(expected, jrec['extras']['qcvars']['2-BODY DISPERSION CORRECTION ENERGY'], atol=1.e-7)
    assert compare_values(expected, jrec['extras']['qcvars'][inp['lbl'] + ' DISPERSION CORRECTION ENERGY'], atol=1.e-7)

    assert compare_values(gexpected, jrec['extras']['qcvars']['CURRENT GRADIENT'], atol=1.e-7)
    assert compare_values(gexpected, jrec['extras']['qcvars']['DISPERSION CORRECTION GRADIENT'], atol=1.e-7)
    assert compare_values(gexpected, jrec['extras']['qcvars']['2-BODY DISPERSION CORRECTION GRADIENT'], atol=1.e-7)
    assert compare_values(
        gexpected, jrec['extras']['qcvars'][inp['lbl'] + ' DISPERSION CORRECTION GRADIENT'], atol=1.e-7)


@using_dftd3_321
@pytest.mark.parametrize(
    "subjects",
    [
        pytest.param(eneyne_ne_psi4mols, marks=using_psi4),
        pytest.param(eneyne_ne_qcdbmols,
                     marks=using_psi4),  # needs qcdb.Molecule, presently more common in psi4 than in qcdb
        pytest.param(eneyne_ne_qcschemamols),
    ],
    ids=['qmol', 'pmol', 'qcmol'])
@pytest.mark.parametrize("inp", [
    ({'parent': 'eneyne', 'name': 'd3-atmgr', 'subject': 'dimer', 'lbl': 'ATM'}),
    ({'parent': 'eneyne', 'name': 'd3-b3lyp-atmgr', 'subject': 'mA', 'lbl': 'ATM'}),
    ({'parent': 'eneyne', 'name': 'd3-pbe-atm(gr)', 'subject': 'mB', 'lbl': 'ATM'}),
    ({'parent': 'eneyne', 'name': 'd3-ATMgr', 'subject': 'mAgB', 'lbl': 'ATM'}),
    ({'parent': 'eneyne', 'name': 'd3-atmgr', 'subject': 'gAmB', 'lbl': 'ATM'}),
    ({'parent': 'ne', 'name': 'd3-atmgr', 'subject': 'atom', 'lbl': 'ATM'}),
])  # yapf: disable
def test_dftd3__run_dftd3__3body(inp, subjects, request):
    subject = subjects()[inp['parent']][inp['subject']]
    expected = ref[inp['parent']][inp['lbl']][inp['subject']]
    gexpected = gref[inp['parent']][inp['lbl']][inp['subject']].ravel()

    if 'qcmol' in request.node.name:
        mol = subject
    else:
        mol = subject.to_schema(dtype=2)

    resinp = {
        'schema_name': 'qcschema_input',
        'schema_version': 1,
        'molecule': mol,
        'driver': 'gradient',
        'model': {
            'method': inp['name']
        },
        'keywords': {},
    }
    jrec = qcng.compute(resinp, 'dftd3', raise_error=True)
    jrec = jrec.dict()

    assert len(jrec['extras']['qcvars']) == 8

    assert compare_values(expected, jrec['extras']['qcvars']['CURRENT ENERGY'], atol=1.e-7)
    assert compare_values(expected, jrec['extras']['qcvars']['DISPERSION CORRECTION ENERGY'], atol=1.e-7)
    assert compare_values(expected, jrec['extras']['qcvars']['3-BODY DISPERSION CORRECTION ENERGY'], atol=1.e-7)
    assert compare_values(
        expected, jrec['extras']['qcvars']['AXILROD-TELLER-MUTO 3-BODY DISPERSION CORRECTION ENERGY'], atol=1.e-7)

    assert compare_values(gexpected, jrec['extras']['qcvars']['CURRENT GRADIENT'], atol=1.e-7)
    assert compare_values(gexpected, jrec['extras']['qcvars']['DISPERSION CORRECTION GRADIENT'], atol=1.e-7)
    assert compare_values(gexpected, jrec['extras']['qcvars']['3-BODY DISPERSION CORRECTION GRADIENT'], atol=1.e-7)
    assert compare_values(
        gexpected, jrec['extras']['qcvars']['AXILROD-TELLER-MUTO 3-BODY DISPERSION CORRECTION GRADIENT'], atol=1.e-7)
