import pytest
import psi4
import itertools
from utils import compare, compare_integers, compare_values
from addons import using

pytestmark = [pytest.mark.psi, pytest.mark.api]

@pytest.mark.parametrize(
    "scf_type,scf_subtype", [
        # avoid redundant combinations of type and subtype
        pytest.param("PK", "INCORE"),
        pytest.param("PK", "OUT_OF_CORE"),
        pytest.param("PK", "YOSHIMINE_OUT_OF_CORE"),
        pytest.param("PK", "REORDER_OUT_OF_CORE"),

        pytest.param("DIRECT", "AUTO"),

        pytest.param("OUT_OF_CORE", "AUTO"),

        pytest.param("DISK_DF", "INCORE"),
        pytest.param("DISK_DF", "OUT_OF_CORE"),

        pytest.param("MEM_DF", "INCORE"),
        pytest.param("MEM_DF", "OUT_OF_CORE"),

        pytest.param("DFDIRJ+LINK", "AUTO"),

        pytest.param("DFDIRJ+COSX", "AUTO"),

        pytest.param("DFDIRJ+SNLINK", "AUTO", marks=using('gauxc')),
    ]
)
@pytest.mark.parametrize("screening", [ "SCHWARZ", "DENSITY", "CSAM", "NONE" ])
def test_comprehensive_jk_screening(scf_type, scf_subtype, screening):
    """Checks the energy values computed by different JK methods using different
    screening types. The differences in energies should be insignificant.""" 

    #== define reference energies ==#
    Eref = {  
        "Nuclear"       :   30.7884922572,
        "Singlet": {
            "Canonical" : -149.58723684929720,
            "DF"        : -149.58715054487624,
            "Composite": {
              "DFDIRJ+COSX"    : -149.58722317236171,
              "DFDIRJ+LINK"    : -149.58726772171027,
              "DFDIRJ+SNLINK"  : -149.58726759282922, 
              "DFDIRJ+COSX_OOO": -149.58726961908508,  # like scf5, OOO finds lower root
            }
        }
    }
    
    #== define molecule ==#
    singlet_o2 = """ 
        0 1
        O
        O 1 1.1
        units    angstrom
    """
    
    psi4.geometry(singlet_o2)
    
    #=== define options ==#
    psi4.set_options({
        "scf_type": scf_type,
        "scf_subtype": scf_subtype,
        "screening": screening, 
        "df_scf_guess": False,
        "basis": "cc-pvtz",
        "df_basis_scf": "cc-pvtz-jkfit",
        "print": 2,
    })

    #== certain combinations of SCF_TYPE and SCREENING should throw an exception by design ==#
    should_throw = False
    #== for now, this is non-integral-direct methods and DFDirJ+COSX with SCREENING=DENSITY... ==#
    should_throw = should_throw or (scf_type not in [ "DIRECT", "DFDIRJ+LINK" ] and screening == "DENSITY")

    E = 0.0 
    
    #== if expected, test if current option combo throws exception ==# 
    if should_throw:
        with pytest.raises(Exception) as e_info:
            E = psi4.energy('scf')

        # we keep this line just for printout purposes; should always pass if done correctly
        assert compare(type(e_info), pytest.ExceptionInfo, f'Singlet {scf_type}({scf_subtype})+{screening}  throws exception')

    #== otherwise, test if current option combo gives right answer ==#
    else: 
        E = psi4.energy('scf')
  
        E_ref = 0.0
        if scf_type in [ "DIRECT", "PK", "OUT_OF_CORE" ]:
            E_ref = Eref["Singlet"]["Canonical"] 
        elif scf_type in [ "MEM_DF", "DISK_DF" ]:
            E_ref = Eref["Singlet"]["DF"]
        elif scf_type in Eref["Singlet"]["Composite"].keys(): 
            if psi4.core.get_global_option("orbital_optimizer_package") != "INTERNAL" and scf_type == "DFDIRJ+COSX":  # KP-DIFF-ANS
                E_ref = Eref["Singlet"]["Composite"][f"{scf_type}_OOO"]
            else:
                E_ref = Eref["Singlet"]["Composite"][scf_type]
        else:
            raise Exception("Invalid JK method used!") 

        assert compare_values(E_ref, E, 6, f'Singlet {scf_type}({scf_subtype})+{screening}  RHF energy')
