from pathlib import Path

import numpy as np
import pytest

import psi4
from addons import *
from psi4.driver.p4util.solvers import davidson_solver, hamiltonian_solver
from psi4.driver.procrouting.response.scf_products import (TDRSCFEngine,
                                                           TDUSCFEngine)
from utils import *

## marks
# reference type
UHF = pytest.mark.unrestricted
RHF_singlet = pytest.mark.restricted_singlet
RHF_triplet = pytest.mark.restricted_triplet
# functional types
hf = pytest.mark.hf
lda = pytest.mark.lda
gga = pytest.mark.gga
hyb_gga = pytest.mark.hyb_gga
hyb_gga_lrc = pytest.mark.hyb_gga_lrc
# response type
RPA = pytest.mark.RPA
TDA = pytest.mark.TDA


@pytest.fixture
def tddft_systems():
    psi4.core.clean()

    # Canonical unrestricted system
    ch2 = psi4.geometry("""0 3
    C           0.000000    0.000000    0.159693
    H          -0.000000    0.895527   -0.479080
    H          -0.000000   -0.895527   -0.479080
    no_reorient
    no_com
    """)

    # Canonical restricted system
    h2o = psi4.geometry("""0 1
    O           0.000000    0.000000    0.135446
    H          -0.000000    0.866812   -0.541782
    H          -0.000000   -0.866812   -0.541782
    no_reorient
    no_com
    """)

    return {'UHF': ch2, 'RHF': h2o}


@pytest.fixture
def wfn_factory(tddft_systems):
    def _build_wfn(ref, func, basis, nosym):
        if ref.startswith('RHF'):
            mol = tddft_systems['RHF']
        else:
            mol = tddft_systems['UHF']
            psi4.set_options({'reference': 'UHF'})
        if nosym:
            mol.reset_point_group('c1')
        psi4.set_options({'scf_type': 'pk', 'e_convergence': 8, 'd_convergence': 8, 'save_jk': True})
        e, wfn = psi4.energy("{}/{}".format(func, basis), return_wfn=True, molecule=mol)
        return wfn

    return _build_wfn


@pytest.fixture
def solver_funcs():
    return {'TDA': davidson_solver, 'RPA': hamiltonian_solver}


@pytest.fixture
def engines():
    return {
        'RHF-1': lambda w, p: TDRSCFEngine(w, ptype=p.lower(), triplet=False),
        'RHF-3': lambda w, p: TDRSCFEngine(w, ptype=p.lower(), triplet=True),
        'UHF': lambda w, p: TDUSCFEngine(w, ptype=p.lower())
    }


@pytest.fixture
def expected():
    """Fixture holding expected values"""
    return {
        "UHF-SVWN-RPA-cc-pvdz": [{
            "e": 0.18569065156695178,
            "sym": "B2"
        }, {
            "e": 0.2371470529055981,
            "sym": "A2"
        }, {
            "e": 0.2606702942260291,
            "sym": "B1"
        }, {
            "e": 0.304699657491157,
            "sym": "A2"
        }, {
            "e": 0.3264956806549558,
            "sym": "A1"
        }, {
            "e": 0.37615871505650317,
            "sym": "B2"
        }, {
            "e": 0.39731897494137147,
            "sym": "A1"
        }, {
            "e": 0.4074250386237938,
            "sym": "B2"
        }, {
            "e": 0.44300573265041665,
            "sym": "B2"
        }, {
            "e": 0.4471437063618376,
            "sym": "A1"
        }],
        "UHF-SVWN-TDA-cc-pvdz": [{
            "e": 0.18742889452032843,
            "sym": "B2"
        }, {
            "e": 0.2379445132034474,
            "sym": "A2"
        }, {
            "e": 0.2614052806756598,
            "sym": "B1"
        }, {
            "e": 0.3049128035615499,
            "sym": "A2"
        }, {
            "e": 0.3271020444759012,
            "sym": "A1"
        }, {
            "e": 0.37727956939219004,
            "sym": "B2"
        }, {
            "e": 0.4011813287341809,
            "sym": "A1"
        }, {
            "e": 0.40830334743110247,
            "sym": "B2"
        }, {
            "e": 0.44456757885588194,
            "sym": "B2"
        }, {
            "e": 0.44881580053474746,
            "sym": "A1"
        }],
        "RHF-1-SVWN-RPA-cc-pvdz": [{
            "e": 0.22569596402035155,
            "sym": "B1"
        }, {
            "e": 0.2900256530242798,
            "sym": "A2"
        }, {
            "e": 0.32071133729636214,
            "sym": "A1"
        }, {
            "e": 0.38457063497252675,
            "sym": "B2"
        }, {
            "e": 0.4422266470138081,
            "sym": "B2"
        }, {
            "e": 0.5495126190664035,
            "sym": "A1"
        }, {
            "e": 0.6878003195644223,
            "sym": "A2"
        }, {
            "e": 0.7199890511259994,
            "sym": "B1"
        }, {
            "e": 0.798062986738022,
            "sym": "B2"
        }, {
            "e": 0.804067826031505,
            "sym": "A1"
        }],
        "RHF-1-SVWN-TDA-cc-pvdz": [{
            "e": 0.2272761848870576,
            "sym": "B1"
        }, {
            "e": 0.290485019555299,
            "sym": "A2"
        }, {
            "e": 0.32546669962547287,
            "sym": "A1"
        }, {
            "e": 0.38741503253259757,
            "sym": "B2"
        }, {
            "e": 0.44712165676834864,
            "sym": "B2"
        }, {
            "e": 0.5652927781399749,
            "sym": "A1"
        }, {
            "e": 0.6879399669898522,
            "sym": "A2"
        }, {
            "e": 0.721392875244794,
            "sym": "B1"
        }, {
            "e": 0.8013924753548491,
            "sym": "B2"
        }, {
            "e": 0.8054973746760368,
            "sym": "A1"
        }],
        "RHF-3-SVWN-RPA-cc-pvdz": [{
            "e": 0.19919602757891613,
            "sym": "B1"
        }, {
            "e": 0.2697069526242387,
            "sym": "A2"
        }, {
            "e": 0.2731025900215326,
            "sym": "A1"
        }, {
            "e": 0.3404714479946839,
            "sym": "B2"
        }, {
            "e": 0.39101279120353993,
            "sym": "B2"
        }, {
            "e": 0.4437113196420621,
            "sym": "A1"
        }, {
            "e": 0.6798587909761624,
            "sym": "A2"
        }, {
            "e": 0.6987920419186496,
            "sym": "B1"
        }, {
            "e": 0.7436776643975973,
            "sym": "B2"
        }, {
            "e": 0.7614459618174199,
            "sym": "A1"
        }],
        "RHF-3-SVWN-TDA-cc-pvdz": [{
            "e": 0.19988323990932086,
            "sym": "B1"
        }, {
            "e": 0.2701736690197542,
            "sym": "A2"
        }, {
            "e": 0.2748922820263834,
            "sym": "A1"
        }, {
            "e": 0.3421398672353456,
            "sym": "B2"
        }, {
            "e": 0.39202707250403035,
            "sym": "B2"
        }, {
            "e": 0.44628377221576965,
            "sym": "A1"
        }, {
            "e": 0.6799727138758552,
            "sym": "A2"
        }, {
            "e": 0.6990125378535389,
            "sym": "B1"
        }, {
            "e": 0.7447874939365398,
            "sym": "B2"
        }, {
            "e": 0.7633679513832042,
            "sym": "A1"
        }],
        "UHF-HF-RPA-cc-pvdz": [{
            "e": 0.2445704160468683,
            "sym": "B2"
        }, {
            "e": 0.2878574429978692,
            "sym": "A2"
        }, {
            "e": 0.3179110389232691,
            "sym": "B1"
        }, {
            "e": 0.3547301851175197,
            "sym": "A2"
        }, {
            "e": 0.3879221731828428,
            "sym": "A1"
        }, {
            "e": 0.4038089052916107,
            "sym": "B2"
        }, {
            "e": 0.43529939972603865,
            "sym": "A1"
        }, {
            "e": 0.4508039388809985,
            "sym": "B2"
        }, {
            "e": 0.4834961361605727,
            "sym": "B2"
        }, {
            "e": 0.515831865012076,
            "sym": "A1"
        }],
        "UHF-HF-TDA-cc-pvdz": [{
            "e": 0.24880026306449304,
            "sym": "B2"
        }, {
            "e": 0.28968755925744966,
            "sym": "A2"
        }, {
            "e": 0.32054964027744337,
            "sym": "B1"
        }, {
            "e": 0.35741288565867185,
            "sym": "A2"
        }, {
            "e": 0.39502214228627547,
            "sym": "A1"
        }, {
            "e": 0.41144173957102564,
            "sym": "B2"
        }, {
            "e": 0.4445528791268893,
            "sym": "A1"
        }, {
            "e": 0.4535932124573471,
            "sym": "B2"
        }, {
            "e": 0.48482278670215617,
            "sym": "B2"
        }, {
            "e": 0.5203446818128086,
            "sym": "A1"
        }],
        "RHF-1-HF-RPA-cc-pvdz": [{
            "e": 0.27878771020942616,
            "sym": "B1"
        }, {
            "e": 0.333654448674359,
            "sym": "A2"
        }, {
            "e": 0.3755413264388134,
            "sym": "A1"
        }, {
            "e": 0.42481114308980833,
            "sym": "B2"
        }, {
            "e": 0.45601499280888025,
            "sym": "B2"
        }, {
            "e": 0.5776368755615227,
            "sym": "A1"
        }, {
            "e": 0.8037995559773897,
            "sym": "A2"
        }, {
            "e": 0.839747743228828,
            "sym": "B1"
        }, {
            "e": 0.9021958669217016,
            "sym": "B2"
        }, {
            "e": 0.9229188098690396,
            "sym": "A1"
        }],
        "RHF-1-HF-TDA-cc-pvdz": [{
            "e": 0.282253171319435,
            "sym": "B1"
        }, {
            "e": 0.33726690707429396,
            "sym": "A2"
        }, {
            "e": 0.3798850963561309,
            "sym": "A1"
        }, {
            "e": 0.4300699711369161,
            "sym": "B2"
        }, {
            "e": 0.45887776503019195,
            "sym": "B2"
        }, {
            "e": 0.5918368137683881,
            "sym": "A1"
        }, {
            "e": 0.8045639418850057,
            "sym": "A2"
        }, {
            "e": 0.8421474739868723,
            "sym": "B1"
        }, {
            "e": 0.9055290304707769,
            "sym": "B2"
        }, {
            "e": 0.9250906948276983,
            "sym": "A1"
        }],
        "RHF-3-HF-RPA-cc-pvdz": [{
            "e": 0.2357358789223071,
            "sym": "B1"
        }, {
            "e": 0.26804588324807327,
            "sym": "A1"
        }, {
            "e": 0.3013517942130891,
            "sym": "A2"
        }, {
            "e": 0.30942929529453067,
            "sym": "B2"
        }, {
            "e": 0.4023793566470789,
            "sym": "B2"
        }, {
            "e": 0.42341834376775817,
            "sym": "A1"
        }, {
            "e": 0.7930871284740221,
            "sym": "A2"
        }, {
            "e": 0.8015468225092716,
            "sym": "B2"
        }, {
            "e": 0.8053062781991327,
            "sym": "A1"
        }, {
            "e": 0.8117190349721608,
            "sym": "B1"
        }],
        "RHF-3-HF-TDA-cc-pvdz": [{
            "e": 0.2428836221449658,
            "sym": "B1"
        }, {
            "e": 0.2954976271417805,
            "sym": "A1"
        }, {
            "e": 0.30883763120257796,
            "sym": "A2"
        }, {
            "e": 0.33788429569198375,
            "sym": "B2"
        }, {
            "e": 0.4084209452630434,
            "sym": "B2"
        }, {
            "e": 0.4426382394256013,
            "sym": "A1"
        }, {
            "e": 0.7945056523218094,
            "sym": "A2"
        }, {
            "e": 0.808852587818601,
            "sym": "B2"
        }, {
            "e": 0.8142547382233867,
            "sym": "B1"
        }, {
            "e": 0.8170954608512094,
            "sym": "A1"
        }],
        "UHF-HCTH93-RPA-cc-pvdz": [{
            "e": 0.1856355275832295,
            "sym": "B2"
        }, {
            "e": 0.2451767798678136,
            "sym": "A2"
        }, {
            "e": 0.26726679761146477,
            "sym": "B1"
        }, {
            "e": 0.3089515541022707,
            "sym": "A2"
        }, {
            "e": 0.3376895242828316,
            "sym": "A1"
        }, {
            "e": 0.38421049161220766,
            "sym": "B2"
        }, {
            "e": 0.41623385122261786,
            "sym": "B2"
        }, {
            "e": 0.41855640840345093,
            "sym": "A1"
        }, {
            "e": 0.45103913454488037,
            "sym": "B2"
        }, {
            "e": 0.4544715212646558,
            "sym": "A1"
        }],
        "UHF-HCTH93-TDA-cc-pvdz": [{
            "e": 0.18766409018421026,
            "sym": "B2"
        }, {
            "e": 0.24640420723869694,
            "sym": "A2"
        }, {
            "e": 0.26811570696078824,
            "sym": "B1"
        }, {
            "e": 0.3091904246984007,
            "sym": "A2"
        }, {
            "e": 0.33847595978393646,
            "sym": "A1"
        }, {
            "e": 0.3854121944573539,
            "sym": "B2"
        }, {
            "e": 0.41763767534141255,
            "sym": "B2"
        }, {
            "e": 0.4238887350955218,
            "sym": "A1"
        }, {
            "e": 0.4524539834604195,
            "sym": "B2"
        }, {
            "e": 0.4571064476865819,
            "sym": "A1"
        }],
        "RHF-1-HCTH93-RPA-cc-pvdz": [{
            "e": 0.2278274247242806,
            "sym": "B1"
        }, {
            "e": 0.2894192892033345,
            "sym": "A2"
        }, {
            "e": 0.3289064362097446,
            "sym": "A1"
        }, {
            "e": 0.3884917876813066,
            "sym": "B2"
        }, {
            "e": 0.4481984119170576,
            "sym": "B2"
        }, {
            "e": 0.5517690274667697,
            "sym": "A1"
        }, {
            "e": 0.6738466018181832,
            "sym": "A2"
        }, {
            "e": 0.707226011428162,
            "sym": "B1"
        }, {
            "e": 0.7881149451422703,
            "sym": "B2"
        }, {
            "e": 0.7978351409386365,
            "sym": "A1"
        }],
        "RHF-1-HCTH93-TDA-cc-pvdz": [{
            "e": 0.22904382729841943,
            "sym": "B1"
        }, {
            "e": 0.28974635817342015,
            "sym": "A2"
        }, {
            "e": 0.3335662503004033,
            "sym": "A1"
        }, {
            "e": 0.39102014106803623,
            "sym": "B2"
        }, {
            "e": 0.4529464243816721,
            "sym": "B2"
        }, {
            "e": 0.5679313794941492,
            "sym": "A1"
        }, {
            "e": 0.6739678745823723,
            "sym": "A2"
        }, {
            "e": 0.7085967611567234,
            "sym": "B1"
        }, {
            "e": 0.7913415356561492,
            "sym": "B2"
        }, {
            "e": 0.7994116868730944,
            "sym": "A1"
        }],
        "RHF-3-HCTH93-RPA-cc-pvdz": [{
            "e": 0.2006807002071702,
            "sym": "B1"
        }, {
            "e": 0.26766736522651347,
            "sym": "A2"
        }, {
            "e": 0.2785010654940702,
            "sym": "A1"
        }, {
            "e": 0.3420810683193752,
            "sym": "B2"
        }, {
            "e": 0.39419895746268907,
            "sym": "B2"
        }, {
            "e": 0.4394557480987003,
            "sym": "A1"
        }, {
            "e": 0.6662358171322572,
            "sym": "A2"
        }, {
            "e": 0.6870212339278138,
            "sym": "B1"
        }, {
            "e": 0.732255974970336,
            "sym": "B2"
        }, {
            "e": 0.7517257660210537,
            "sym": "A1"
        }],
        "RHF-3-HCTH93-TDA-cc-pvdz": [{
            "e": 0.20200735074875364,
            "sym": "B1"
        }, {
            "e": 0.26861182281428897,
            "sym": "A2"
        }, {
            "e": 0.281727656007949,
            "sym": "A1"
        }, {
            "e": 0.3449548653374312,
            "sym": "B2"
        }, {
            "e": 0.3960474483835103,
            "sym": "B2"
        }, {
            "e": 0.44406043820563673,
            "sym": "A1"
        }, {
            "e": 0.6665408365088539,
            "sym": "A2"
        }, {
            "e": 0.6875430743070515,
            "sym": "B1"
        }, {
            "e": 0.7339868680592164,
            "sym": "B2"
        }, {
            "e": 0.7549670562639251,
            "sym": "A1"
        }],
        "UHF-LRC-wPBE-RPA-cc-pvdz": [{
            "e": 0.20857445467620409,
            "sym": "B2"
        }, {
            "e": 0.25750617756036887,
            "sym": "A2"
        }, {
            "e": 0.29567402388969183,
            "sym": "B1"
        }, {
            "e": 0.3353963665599838,
            "sym": "A2"
        }, {
            "e": 0.3676659466310203,
            "sym": "A1"
        }, {
            "e": 0.4059771153180213,
            "sym": "B2"
        }, {
            "e": 0.42314272384914664,
            "sym": "A1"
        }, {
            "e": 0.4309152055539914,
            "sym": "B2"
        }, {
            "e": 0.47108588995855805,
            "sym": "B2"
        }, {
            "e": 0.47426103142096276,
            "sym": "A1"
        }],
        "UHF-LRC-wPBE-TDA-cc-pvdz": [{
            "e": 0.2111138328596782,
            "sym": "B2"
        }, {
            "e": 0.2588953019501709,
            "sym": "A2"
        }, {
            "e": 0.296850002209101,
            "sym": "B1"
        }, {
            "e": 0.33591453200697347,
            "sym": "A2"
        }, {
            "e": 0.36902567156283717,
            "sym": "A1"
        }, {
            "e": 0.40870758997839934,
            "sym": "B2"
        }, {
            "e": 0.4285081249314508,
            "sym": "A1"
        }, {
            "e": 0.4325872997269013,
            "sym": "B2"
        }, {
            "e": 0.4722361437522301,
            "sym": "B2"
        }, {
            "e": 0.4778293906339199,
            "sym": "A1"
        }],
        "RHF-1-LRC-wPBE-RPA-cc-pvdz": [{
            "e": 0.24823799843052555,
            "sym": "B1"
        }, {
            "e": 0.3133577978678068,
            "sym": "A2"
        }, {
            "e": 0.34407655653012253,
            "sym": "A1"
        }, {
            "e": 0.4111771444491586,
            "sym": "B2"
        }, {
            "e": 0.45620241435353615,
            "sym": "B2"
        }, {
            "e": 0.5713821408751654,
            "sym": "A1"
        }, {
            "e": 0.7121320459794471,
            "sym": "A2"
        }, {
            "e": 0.7467646074860463,
            "sym": "B1"
        }, {
            "e": 0.8237985372718415,
            "sym": "B2"
        }, {
            "e": 0.8347057361843613,
            "sym": "A1"
        }],
        "RHF-1-LRC-wPBE-TDA-cc-pvdz": [{
            "e": 0.25019673731879144,
            "sym": "B1"
        }, {
            "e": 0.3137951148053371,
            "sym": "A2"
        }, {
            "e": 0.3480087340356468,
            "sym": "A1"
        }, {
            "e": 0.4139259937707775,
            "sym": "B2"
        }, {
            "e": 0.4595355779026114,
            "sym": "B2"
        }, {
            "e": 0.5866661840952361,
            "sym": "A1"
        }, {
            "e": 0.7122790432693733,
            "sym": "A2"
        }, {
            "e": 0.7482823545045337,
            "sym": "B1"
        }, {
            "e": 0.8268119817153274,
            "sym": "B2"
        }, {
            "e": 0.8363484308992859,
            "sym": "A1"
        }],
        "RHF-3-LRC-wPBE-RPA-cc-pvdz": [{
            "e": 0.21493576239775788,
            "sym": "B1"
        }, {
            "e": 0.28128298920592243,
            "sym": "A1"
        }, {
            "e": 0.2873576522121203,
            "sym": "A2"
        }, {
            "e": 0.34026197685653914,
            "sym": "B2"
        }, {
            "e": 0.40093878320580273,
            "sym": "B2"
        }, {
            "e": 0.44324092831429845,
            "sym": "A1"
        }, {
            "e": 0.7026102965244811,
            "sym": "A2"
        }, {
            "e": 0.7221756358136507,
            "sym": "B1"
        }, {
            "e": 0.7560327866158897,
            "sym": "B2"
        }, {
            "e": 0.7725589569358363,
            "sym": "A1"
        }],
        "RHF-3-LRC-wPBE-TDA-cc-pvdz": [{
            "e": 0.2164314598227564,
            "sym": "B1"
        }, {
            "e": 0.28700485871629755,
            "sym": "A1"
        }, {
            "e": 0.28869165261820007,
            "sym": "A2"
        }, {
            "e": 0.34761184135284623,
            "sym": "B2"
        }, {
            "e": 0.40273215014290165,
            "sym": "B2"
        }, {
            "e": 0.44947361340716696,
            "sym": "A1"
        }, {
            "e": 0.7029814646815448,
            "sym": "A2"
        }, {
            "e": 0.722752600176611,
            "sym": "B1"
        }, {
            "e": 0.7587595863440195,
            "sym": "B2"
        }, {
            "e": 0.7765352336283384,
            "sym": "A1"
        }],
        "UHF-PBE0-RPA-cc-pvdz": [{
            "e": 0.2057925309643518,
            "sym": "B2"
        }, {
            "e": 0.25788469558192867,
            "sym": "A2"
        }, {
            "e": 0.27774770438319873,
            "sym": "B1"
        }, {
            "e": 0.3209796073504773,
            "sym": "A2"
        }, {
            "e": 0.34905608972637064,
            "sym": "A1"
        }, {
            "e": 0.3920417722330229,
            "sym": "B2"
        }, {
            "e": 0.4166160441764259,
            "sym": "B2"
        }, {
            "e": 0.4181815653141393,
            "sym": "A1"
        }, {
            "e": 0.4582640513447503,
            "sym": "B2"
        }, {
            "e": 0.46098717614063206,
            "sym": "A1"
        }],
        "UHF-PBE0-TDA-cc-pvdz": [{
            "e": 0.20835395874131485,
            "sym": "B2"
        }, {
            "e": 0.2593583434134383,
            "sym": "A2"
        }, {
            "e": 0.2788171096674114,
            "sym": "B1"
        }, {
            "e": 0.32146837333948175,
            "sym": "A2"
        }, {
            "e": 0.3502541176392687,
            "sym": "A1"
        }, {
            "e": 0.3942357067851706,
            "sym": "B2"
        }, {
            "e": 0.41852333401321756,
            "sym": "B2"
        }, {
            "e": 0.42350654214171385,
            "sym": "A1"
        }, {
            "e": 0.45957232722509295,
            "sym": "B2"
        }, {
            "e": 0.46426154077373694,
            "sym": "A1"
        }],
        "RHF-1-PBE0-RPA-cc-pvdz": [{
            "e": 0.24152757214539713,
            "sym": "B1"
        }, {
            "e": 0.30322600965964747,
            "sym": "A2"
        }, {
            "e": 0.3397842356642792,
            "sym": "A1"
        }, {
            "e": 0.3989212454015664,
            "sym": "B2"
        }, {
            "e": 0.450576093081613,
            "sym": "B2"
        }, {
            "e": 0.5607321872200163,
            "sym": "A1"
        }, {
            "e": 0.7111912633239198,
            "sym": "A2"
        }, {
            "e": 0.7451917364838366,
            "sym": "B1"
        }, {
            "e": 0.8210717375437115,
            "sym": "B2"
        }, {
            "e": 0.8330520166726922,
            "sym": "A1"
        }],
        "RHF-1-PBE0-TDA-cc-pvdz": [{
            "e": 0.24295344585768072,
            "sym": "B1"
        }, {
            "e": 0.3037111007164037,
            "sym": "A2"
        }, {
            "e": 0.3439369091046927,
            "sym": "A1"
        }, {
            "e": 0.40143489905930346,
            "sym": "B2"
        }, {
            "e": 0.4543833228907001,
            "sym": "B2"
        }, {
            "e": 0.5759280320661313,
            "sym": "A1"
        }, {
            "e": 0.711338260613846,
            "sym": "A2"
        }, {
            "e": 0.7466102603316238,
            "sym": "B1"
        }, {
            "e": 0.8240815070549493,
            "sym": "B2"
        }, {
            "e": 0.8345991631491648,
            "sym": "A1"
        }],
        "RHF-3-PBE0-RPA-cc-pvdz": [{
            "e": 0.20871042716938573,
            "sym": "B1"
        }, {
            "e": 0.2750098798583243,
            "sym": "A1"
        }, {
            "e": 0.2767775222696862,
            "sym": "A2"
        }, {
            "e": 0.33116284461011086,
            "sym": "B2"
        }, {
            "e": 0.3922953425581455,
            "sym": "B2"
        }, {
            "e": 0.43251747601418633,
            "sym": "A1"
        }, {
            "e": 0.7010631500480085,
            "sym": "A2"
        }, {
            "e": 0.7203308198250776,
            "sym": "B1"
        }, {
            "e": 0.752831920627748,
            "sym": "B2"
        }, {
            "e": 0.7687590769912455,
            "sym": "A1"
        }],
        "RHF-3-PBE0-TDA-cc-pvdz": [{
            "e": 0.21080513855083327,
            "sym": "B1"
        }, {
            "e": 0.27866276251298894,
            "sym": "A2"
        }, {
            "e": 0.28137118757987817,
            "sym": "A1"
        }, {
            "e": 0.3390676238758892,
            "sym": "B2"
        }, {
            "e": 0.3944745773813006,
            "sym": "B2"
        }, {
            "e": 0.4398416159847564,
            "sym": "A1"
        }, {
            "e": 0.7015372163080204,
            "sym": "A2"
        }, {
            "e": 0.7211099054616862,
            "sym": "B1"
        }, {
            "e": 0.7557167424425485,
            "sym": "B2"
        }, {
            "e": 0.7733747918949263,
            "sym": "A1"
        }],
        "UHF-wB97X-RPA-cc-pvdz": [{
            "e": 0.18313657365448505,
            "sym": "B2"
        }, {
            "e": 0.24031484450350646,
            "sym": "A2"
        }, {
            "e": 0.2859023790418515,
            "sym": "B1"
        }, {
            "e": 0.3266610526061227,
            "sym": "A2"
        }, {
            "e": 0.36160598335381505,
            "sym": "A1"
        }, {
            "e": 0.4077778321196165,
            "sym": "B2"
        }, {
            "e": 0.4338588262847624,
            "sym": "B2"
        }, {
            "e": 0.4340278731681775,
            "sym": "A1"
        }, {
            "e": 0.4666281971415478,
            "sym": "B2"
        }, {
            "e": 0.48119562857322856,
            "sym": "A1"
        }],
        "UHF-wB97X-TDA-cc-pvdz": [{
            "e": 0.1852129103746918,
            "sym": "B2"
        }, {
            "e": 0.24126665195577823,
            "sym": "A2"
        }, {
            "e": 0.28685786142637143,
            "sym": "B1"
        }, {
            "e": 0.3270946946114049,
            "sym": "A2"
        }, {
            "e": 0.3624254932451533,
            "sym": "A1"
        }, {
            "e": 0.4091853311706593,
            "sym": "B2"
        }, {
            "e": 0.43516710216510507,
            "sym": "B2"
        }, {
            "e": 0.4401539852258495,
            "sym": "A1"
        }, {
            "e": 0.4679474978186349,
            "sym": "B2"
        }, {
            "e": 0.484293596458422,
            "sym": "A1"
        }],
        "RHF-1-wB97X-RPA-cc-pvdz": [{
            "e": 0.2472126923332907,
            "sym": "B1"
        }, {
            "e": 0.3131997757811362,
            "sym": "A2"
        }, {
            "e": 0.34444772468718604,
            "sym": "A1"
        }, {
            "e": 0.41038703401580556,
            "sym": "B2"
        }, {
            "e": 0.45278105243050515,
            "sym": "B2"
        }, {
            "e": 0.5674462884373929,
            "sym": "A1"
        }, {
            "e": 0.7190997175219462,
            "sym": "A2"
        }, {
            "e": 0.7518360139884983,
            "sym": "B1"
        }, {
            "e": 0.8295240817144648,
            "sym": "B2"
        }, {
            "e": 0.840034387944184,
            "sym": "A1"
        }],
        "RHF-1-wB97X-TDA-cc-pvdz": [{
            "e": 0.24888846143844873,
            "sym": "B1"
        }, {
            "e": 0.31362606792192205,
            "sym": "A2"
        }, {
            "e": 0.34848280029565865,
            "sym": "A1"
        }, {
            "e": 0.41301093564098723,
            "sym": "B2"
        }, {
            "e": 0.4563016375242363,
            "sym": "B2"
        }, {
            "e": 0.5822268659394665,
            "sym": "A1"
        }, {
            "e": 0.7192503897441206,
            "sym": "A2"
        }, {
            "e": 0.7532949620910153,
            "sym": "B1"
        }, {
            "e": 0.8325669256159359,
            "sym": "B2"
        }, {
            "e": 0.8418130551522903,
            "sym": "A1"
        }],
        "RHF-3-wB97X-RPA-cc-pvdz": [{
            "e": 0.21928320724732356,
            "sym": "B1"
        }, {
            "e": 0.29090028689934033,
            "sym": "A2"
        }, {
            "e": 0.29346538960855156,
            "sym": "A1"
        }, {
            "e": 0.3566558496155522,
            "sym": "B2"
        }, {
            "e": 0.4047901122018677,
            "sym": "B2"
        }, {
            "e": 0.45354911327036923,
            "sym": "A1"
        }, {
            "e": 0.7135101455725047,
            "sym": "A2"
        }, {
            "e": 0.7332151322871041,
            "sym": "B1"
        }, {
            "e": 0.7665614675068497,
            "sym": "B2"
        }, {
            "e": 0.7836168280705303,
            "sym": "A1"
        }],
        "RHF-3-wB97X-TDA-cc-pvdz": [{
            "e": 0.22053635914394396,
            "sym": "B1"
        }, {
            "e": 0.29207259028650134,
            "sym": "A2"
        }, {
            "e": 0.2973130436723683,
            "sym": "A1"
        }, {
            "e": 0.36116866641628476,
            "sym": "B2"
        }, {
            "e": 0.4065136554262517,
            "sym": "B2"
        }, {
            "e": 0.4588630653011993,
            "sym": "A1"
        }, {
            "e": 0.7138849886618163,
            "sym": "A2"
        }, {
            "e": 0.7338067963790569,
            "sym": "B1"
        }, {
            "e": 0.7689722230616385,
            "sym": "B2"
        }, {
            "e": 0.7875931047630326,
            "sym": "A1"
        }]
    }


@pytest.mark.tdscf
@pytest.mark.parametrize("ref,func,ptype,basis", [
    pytest.param(   'UHF',      'SVWN',  'RPA',  'cc-pvdz', marks=[lda, UHF, RPA]), # G09 rev E.01
    pytest.param(   'UHF',      'SVWN',  'TDA',  'cc-pvdz', marks=[lda, UHF, TDA]), # G09 rev E.01
    pytest.param( 'RHF-1',      'SVWN',  'RPA',  'cc-pvdz', marks=[lda, RHF_singlet, RPA]), # G09 rev E.01
    pytest.param( 'RHF-1',      'SVWN',  'TDA',  'cc-pvdz', marks=[lda, RHF_singlet, TDA]), # G09 rev E.01
    pytest.param( 'RHF-3',      'SVWN',  'RPA',  'cc-pvdz', marks=[lda, RHF_triplet, RPA]), # G09 rev E.01
    pytest.param( 'RHF-3',      'SVWN',  'TDA',  'cc-pvdz', marks=[lda, RHF_triplet, TDA]), # G09 rev E.01
    pytest.param(   'UHF',        'HF',  'RPA',  'cc-pvdz', marks=[hf, UHF, RPA, pytest.mark.quick]), # G09 rev E.01
    pytest.param(   'UHF',        'HF',  'TDA',  'cc-pvdz', marks=[hf, UHF, TDA, pytest.mark.quick]), # G09 rev E.01
    pytest.param( 'RHF-1',        'HF',  'RPA',  'cc-pvdz', marks=[hf, RHF_singlet, RPA, pytest.mark.quick]), # G09 rev E.01
    pytest.param( 'RHF-1',        'HF',  'TDA',  'cc-pvdz', marks=[hf, RHF_singlet, TDA, pytest.mark.quick]), # G09 rev E.01
    pytest.param( 'RHF-3',        'HF',  'RPA',  'cc-pvdz', marks=[hf, RHF_triplet, RPA, pytest.mark.quick]), # G09 rev E.01
    pytest.param( 'RHF-3',        'HF',  'TDA',  'cc-pvdz', marks=[hf, RHF_triplet, TDA, pytest.mark.quick]), # G09 rev E.01
    pytest.param(   'UHF',    'HCTH93',  'RPA',  'cc-pvdz', marks=[gga, UHF, RPA]), # G09 rev E.01
    pytest.param(   'UHF',    'HCTH93',  'TDA',  'cc-pvdz', marks=[gga, UHF, TDA]), # G09 rev E.01
    pytest.param( 'RHF-1',    'HCTH93',  'RPA',  'cc-pvdz', marks=[gga, RHF_singlet, RPA]), # G09 rev E.01
    pytest.param( 'RHF-1',    'HCTH93',  'TDA',  'cc-pvdz', marks=[gga, RHF_singlet, TDA]), # G09 rev E.01
    pytest.param( 'RHF-3',    'HCTH93',  'RPA',  'cc-pvdz', marks=[gga, RHF_triplet, RPA]), # G09 rev E.01
    pytest.param( 'RHF-3',    'HCTH93',  'TDA',  'cc-pvdz', marks=[gga, RHF_triplet, TDA]), # G09 rev E.01
    pytest.param(   'UHF',      'PBE0',  'RPA',  'cc-pvdz', marks=[hyb_gga, UHF, RPA]), # G09 rev E.01
    pytest.param(   'UHF',      'PBE0',  'TDA',  'cc-pvdz', marks=[hyb_gga, UHF, TDA]), # G09 rev E.01
    pytest.param( 'RHF-1',      'PBE0',  'RPA',  'cc-pvdz', marks=[hyb_gga, RHF_singlet, RPA]), # G09 rev E.01
    pytest.param( 'RHF-1',      'PBE0',  'TDA',  'cc-pvdz', marks=[hyb_gga, RHF_singlet, TDA]), # G09 rev E.01
    pytest.param( 'RHF-3',      'PBE0',  'RPA',  'cc-pvdz', marks=[hyb_gga, RHF_triplet, RPA]), # G09 rev E.01
    pytest.param( 'RHF-3',      'PBE0',  'TDA',  'cc-pvdz', marks=[hyb_gga, RHF_triplet, TDA]), # G09 rev E.01
    pytest.param(   'UHF',     'wB97X',  'RPA',  'cc-pvdz', marks=[hyb_gga_lrc, UHF, RPA]), # G09 rev E.01
    pytest.param(   'UHF',     'wB97X',  'TDA',  'cc-pvdz', marks=[hyb_gga_lrc, UHF, TDA]), # G09 rev E.01
    pytest.param( 'RHF-1',     'wB97X',  'RPA',  'cc-pvdz', marks=[hyb_gga_lrc, RHF_singlet, RPA]), # G09 rev E.01
    pytest.param( 'RHF-1',     'wB97X',  'TDA',  'cc-pvdz', marks=[hyb_gga_lrc, RHF_singlet, TDA]), # G09 rev E.01
    pytest.param( 'RHF-3',     'wB97X',  'RPA',  'cc-pvdz', marks=[hyb_gga_lrc, RHF_triplet, RPA]), # G09 rev E.01
    pytest.param( 'RHF-3',     'wB97X',  'TDA',  'cc-pvdz', marks=[hyb_gga_lrc, RHF_triplet, TDA]), # G09 rev E.01
]) # yapf: disable
def test_tdscf(ref, func, ptype, basis, expected, wfn_factory, solver_funcs, engines):
    if (ref == 'RHF-1') or (func == "HF"):
        # RHF-singlet everything works and TDHF/CIS works for RHF-triplet, UHF
        pass
    elif (ref == 'RHF-3'):
        pytest.xfail("RKS Vx kernel only Spin Adapted for Singlet")
    elif (ref == 'UHF' and func != 'SVWN'):
        pytest.xfail("UKS Vx kernel bug for non-lda")

    ### setup
    # Look up expected in fixture (easier to read)
    exp_lookup = "{}-{}-{}-{}".format(ref, func, ptype, basis)
    exp_low_per_sym = {'A1': 100, 'A2': 100, 'B1': 100, "B2": 100}
    for x in expected[exp_lookup]:
        exp_low_per_sym[x['sym']] = min(x['e'], exp_low_per_sym[x['sym']])

    exp_energy_sorted = np.array([x['e'] for x in expected[exp_lookup]])
    exp_energy_sorted.sort()

    # get wfn (don't use symmetry b/c too slow)
    wfn = wfn_factory(ref, func, basis, nosym=True)
    # select solver function (TDA->davidson/RPA->hamiltonian)
    solver = solver_funcs[ptype]

    # build engine
    engine = engines[ref](wfn, ptype)

    # skipping the entrypoint, just call the solver
    out = solver(
        engine=engine,
        guess=engine.generate_guess(16),
        max_vecs_per_root=10,
        nroot=4,
        verbose=1,
        maxiter=30,
        e_tol=1.0e-7,
        r_tol=1.0e-5,
        schmidt_tol=1.0e-12)

    test_vals = out[0]
    stats = out[-1]
    assert stats[-1]['done'], "Solver did not converge"

    for i, my_v in enumerate(test_vals):
        ref_v = exp_energy_sorted[i]
        assert compare_values(ref_v, my_v, 4, "{}-{}-{}-ROOT-{}".format(ref, func, ptype, i + 1))
