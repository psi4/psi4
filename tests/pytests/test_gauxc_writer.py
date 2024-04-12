import pytest

import os 
import sys

from utils import compare
from addons import uusing

import psi4 

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.long]

@pytest.fixture
def mols():
    return {
        "benzene" : psi4.geometry(""" 
                0 1
                H   1.0131690738465817   -1.5835021154162925   -1.7978401310386911
                H   -1.0131690738465817   1.5835021154162925   -5.012160365319985
                H   -1.8164621323959855   -0.7357020536229176   -4.931121359413318
                H   1.8164621323959855   0.7357020536229176   -1.8788791369453572
                H   -0.7729940563410085   -2.3568841717855844   -3.3212372420741265
                H   0.7729940563410085   2.3568841717855844   -3.4887632542845504
                C   0.5690300414747392   -0.9024360657756142   -2.502675182411813
                C   -0.5690300414747392   0.9024360657756142   -4.307325313946864
                C   -1.0176030741697608   -0.4210740306907108   -4.263060310720532
                C   1.0176030741697608   0.4210740306907108   -2.546940185638147
                C   -0.44857303269502175   -1.3122060956424129   -3.358011244754464
                C   0.44857303269502175   1.3122060956424129   -3.4519892516042137
                symmetry c1
                no_reorient
                no_com
            """),
    }

@uusing("h5py")
@uusing("gauxc")
@pytest.mark.parametrize( 
    "inp",
    [
        pytest.param("benzene", id="benzene"),
    ],
)
@pytest.mark.parametrize(
    "basis",
    [
        pytest.param("bse:cc-pVDZ", id="cc-pVDZ"),
        #pytest.param("bse:cc-pVTZ", id="cc-pVTZ"),
    ],
)
#@pytest.mark.parametrize(
#    "snlink_force_cartesian",
#    [
#        pytest.param(True, id="cartesian"),
#        pytest.param(False, id="spherical"),
#    ],
#)
def test_gauxc_writer(inp, basis, mols, request):
    """Writes GauXC HDF5 data files using data from Psi4 calculations,
    then compares the resulting files to reference GauXC data files.
    All data (molecule, basis set, density, exchange) should match to
    within tolerance."""

    #== have to import some modules here ==#
    sys.path.append(os.path.join(os.path.dirname(__file__), "test_gauxc_writer/")) 
    import gauxc_writer as gxcw 

    #== initial set up ==#
    psi4.core.be_quiet()
    
    test_id = request.node.callspec.id
    
    psi4.set_options({
        "scf_type": "dfdirj+snlink",
        "basis": basis,
        "guess": "core",
        "df_scf_guess": False,
        "e_convergence": 1.0e-10,
        "incfock": False,
        "maxiter": 1,
        "fail_on_maxiter": False,
        "ints_tolerance": 1e-10,
        "snlink_ints_tolerance": 1e-10,
        "snlink_basis_tolerance": 2.22e-16,
        "snlink_use_debug_grid": True,
        "snlink_grid_batch_size": 512,
        #"snlink_force_cartesian": snlink_force_cartesian,
        "save_jk": True,
    })
    
    #== execute calculation==#
    E, wfn = psi4.energy("scf", molecule=mols[inp], return_wfn=True)

    #== write results ==#
    file_written = gxcw.write_results(test_id, wfn)
    assert compare(file_written, True, f'{test_id} HDF5 file written')

    #== cross-check against reference results ==#
    mol_same, basis_same, D_same, K_same = gxcw.validate_results(test_id, 14, write_output=True)

    assert compare(mol_same, True, f'{test_id} same molecule as reference')
    assert compare(basis_same, True, f'{test_id} same basis as reference')
    assert compare(D_same, True, f'{test_id} D matrices same within E-14 RMS tolerance')
    assert compare(K_same, True, f'{test_id} K matrices same within E-14 RMS tolerance')
