import pytest
import psi4
import itertools
from utils import compare, compare_integers, compare_values
from addons import using

pytestmark = [pytest.mark.psi, pytest.mark.api] 

@pytest.mark.parametrize(
    "scf_type", [ 
        pytest.param("PK"), 
        pytest.param("DIRECT"), 
        pytest.param("OUT_OF_CORE"), 
        pytest.param("DISK_DF"), 
        pytest.param("MEM_DF"), 
        pytest.param("DFDIRJ+LINK"), 
        pytest.param("DFDIRJ+COSX"),
        pytest.param("DFDIRJ+SNLINK", marks=using('gauxc')),
    ]
)
@pytest.mark.parametrize("scf_subtype", [ "AUTO", "INCORE", "OUT_OF_CORE", "YOSHIMINE_OUT_OF_CORE", "REORDER_OUT_OF_CORE" ])
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
 
    #== skip redundant option combinations based on type/subtype combination ==#   
    if scf_type not in [ "PK", "DISK_DF", "MEM_DF"] and scf_subtype != "AUTO":
        pytest.skip(f'Singlet {scf_type}({scf_subtype})+{screening}  skipped: redundant test') 
    elif scf_type in [ "PK", "DISK_DF", "MEM_DF"] and scf_subtype == "AUTO":
        pytest.skip(f'Singlet {scf_type}({scf_subtype})+{screening}  skipped: redundant test') 
    elif scf_type in [ "DISK_DF", "MEM_DF"] and scf_subtype in [ "YOSHIMINE_OUT_OF_CORE", "REORDER_OUT_OF_CORE" ]:
        pytest.skip(f'Singlet {scf_type}({scf_subtype})+{screening}  skipped: redundant test') 

    #== certain combinations of SCF_TYPE and SCREENING should throw an exception by design ==#
    should_throw = False
    #== specifically, non-integral-direct methods and DFDirJ+COSX/SNLINK with SCREENING = DENSITY... ==# 
    should_throw = should_throw or (scf_type not in [ "DIRECT", "DFDIRJ+LINK" ] and screening == "DENSITY")
    #== ... Composite methods with SCREENING=NONE... ==#
    should_throw = should_throw or (scf_type in Eref["Singlet"]["Composite"].keys() and screening == "NONE")
    #== .. DISK_DF, DIRECT, or PK with SCREENING=NONE ==#
    should_throw = should_throw or (scf_type == "PK" and screening == "NONE")
    should_throw = should_throw or (scf_type == "DISK_DF" and screening == "NONE")
    should_throw = should_throw or (scf_type == "DIRECT" and screening == "NONE")
    #== .. and DFDIRJ+LINK with SCREENING=SCHWARZ or CSAM... ==#
    should_throw = should_throw or (scf_type == "DFDIRJ+LINK" and screening in [ "SCHWARZ", "CSAM" ])
 
    E = 0.0 
    
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
            E_ref = Eref["Singlet"]["Composite"][scf_type]
        else:
            raise Exception("Invalid JK method used!") 

        assert compare_values(E_ref, E, 6, f'Singlet {scf_type}({scf_subtype})+{screening}  RHF energy')
