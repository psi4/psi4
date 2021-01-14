import pytest
import psi4


pytestmark = [pytest.mark.quick]

def test_dft_grid_threaded_raise():
    dimer = psi4.geometry("""
      1 1
      K -4.067042 -1.894214 0.002270
    """)
    
    psi4.set_options({
        "dft_grid_name": "SG1",
        "dft_vv10_radial_points": 50,
        "dft_vv10_spherical_points": 194,
        "dft_nuclear_scheme": "treutler",
        "dft_radial_scheme": "EM",
        "basis": "def2-TZVPPD",
    })
    
    with pytest.raises(RuntimeError) as e:
        ene = psi4.energy("wB97M-V")

    assert "There is no SG-1 grid defined for the requested atomic number" in str(e.value)
