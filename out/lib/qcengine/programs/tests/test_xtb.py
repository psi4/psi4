"""
Testing for the extended tight binding (xtb) program harness.

Most of the tests use mindless molecules for the diversity of element species
to test as much different interactions as possible.
"""

import numpy as np
import pytest
import qcelemental as qcel
from qcelemental.testing import compare_recursive

import qcengine as qcng
from qcengine.testing import using


@using("xtb")
def test_xtb_task_gfn1xtb_m01():

    thr = 1.0e-7

    return_result = np.array(
        [
            [-0.000358030800412827, +0.000500962437171796, -0.009256666294614777],
            [+0.001228927299439432, +0.000322925070197139, +0.002327951513695912],
            [-0.000336383869223903, +0.023413253286994400, -0.003233956461896982],
            [+0.004257485938529293, -0.003712584842768384, +0.000112655448628981],
            [+0.001349484431325273, +0.001591703446467543, -0.001444564232184169],
            [+0.001902096016850091, -0.010806282150986885, -0.002612898358134064],
            [-0.000694337971868170, -0.003558134599236037, +0.002566433298563669],
            [-0.003819207888143676, -0.000529277869208852, +0.004041924390207515],
            [-0.000170900544431631, +0.009573597947311165, +0.008463731042553293],
            [+0.004036505931939967, +0.003025020761338415, -0.004557936832112239],
            [+0.004790375712217767, -0.006895073623195738, +0.006111995172986280],
            [-0.001527116100844721, -0.005269034163367844, -0.008340510436175520],
            [-0.005441057023457164, +0.003102606807119279, +0.000817922187552893],
            [-0.003996677758010005, -0.021839390583571504, +0.010361632086933262],
            [-0.003219908482914395, +0.017669626267348096, -0.011084538678920158],
            [+0.001998745109004805, -0.006589918191612544, +0.005726826152916107],
        ]
    )

    atomic_input = qcel.models.AtomicInput(
        molecule=qcng.get_molecule("mindless-01"),
        model={"method": "GFN1-xTB"},
        driver="gradient",
    )

    atomic_result = qcng.compute(atomic_input, "xtb")

    assert atomic_result.success
    assert pytest.approx(atomic_result.return_result, abs=thr) == return_result


@using("xtb")
def test_xtb_task_gfn1xtb_m02():

    thr = 1.0e-8

    return_result = np.array(
        [
            [+0.0030701065605734, -0.0002868799003678, -0.0027410974319138],
            [-0.0088629069947518, +0.0107644916504071, -0.0022851679479647],
            [-0.0022557927039338, +0.0097846218878783, +0.0082851159007886],
            [-0.0077313130066463, +0.0051272865317215, -0.0113382914323220],
            [-0.0056136016465032, -0.0142776861097860, -0.0087142303224132],
            [-0.0002886024284553, -0.0097634247790545, +0.0076436761747340],
            [+0.0018264919145715, -0.0006763806034670, +0.0014564482472876],
            [+0.0076140201508601, -0.0046761007465963, +0.0001723200093288],
            [-0.0149881572808628, -0.0032937272168366, +0.0287922167009266],
            [+0.0086150890774961, +0.0015081700859065, +0.0027847818912378],
            [+0.0117238707826711, -0.0012461659655145, -0.0006573178440729],
            [-0.0070967442352033, -0.0038359045518911, -0.0046752059918041],
            [+0.0103475500130377, +0.0070999031854269, +0.0075036436034086],
            [+0.0048130271334382, -0.0072967572754967, +0.0021987765873840],
            [-0.0062043239425508, -0.0083018326109514, -0.0061290889642350],
            [+0.0050312866062591, +0.0193703864186214, -0.0222965791803703],
        ]
    )

    atomic_input = qcel.models.AtomicInput(
        molecule=qcng.get_molecule("mindless-02"),
        model={"method": "GFN1-xTB"},
        driver="gradient",
        keywords={
            "accuracy": 0.1,
            "electronic_temperature": 500.0,
        },
    )

    atomic_result = qcng.compute(atomic_input, "xtb")

    assert atomic_result.success
    assert pytest.approx(atomic_result.return_result, abs=thr) == return_result


@using("xtb")
def test_xtb_task_gfn1xtb_m03():

    thr = 1.0e-8

    return_result = np.array(
        [
            [+0.010871358951482300, +0.004305794599339098, -0.019573836237906372],
            [+0.009678292282897322, +0.000426352757690409, +0.002515056015587392],
            [+0.001466416808346703, +0.001883246716371389, -0.000305228304878110],
            [-0.008120495156064160, +0.002564316083979744, -0.005119006489452139],
            [-0.000358357509908501, -0.001936784273655625, -0.002535739405639419],
            [+0.015070726015984316, +0.000278829270926101, +0.009656901918435126],
            [-0.015351695279260991, -0.001257824054680808, +0.018629896846315233],
            [-0.003076325534402210, +0.005496614140839041, +0.005850232528002299],
            [+0.001968655423412034, -0.000723978488125643, -0.000559968512005294],
            [+0.002806198234445627, +0.002329417172370834, +0.001133651685314686],
            [-0.000176086574912473, +0.002193520870220370, -0.001163714333146813],
            [-0.006971737078166070, +0.005435531350132502, +0.001219924082264835],
            [-0.013489242089583017, -0.016619202727944704, -0.012805207719079817],
            [-0.003224316220896839, -0.002988773757177052, +0.004070779421829900],
            [+0.002236894208466977, -0.002020664316295901, +0.002364564453769414],
            [+0.006669713518158995, +0.000633604656010196, -0.003378305949410931],
        ]
    )

    atomic_input = qcel.models.AtomicInput(
        molecule=qcng.get_molecule("mindless-03"),
        model={"method": "GFN1-xTB"},
        driver="gradient",
        keywords={
            "solvent": "chcl3",
        },
    )

    atomic_result = qcng.compute(atomic_input, "xtb")

    assert atomic_result.success
    assert pytest.approx(atomic_result.return_result, abs=thr) == return_result


@using("xtb")
def test_xtb_task_gfn1xtb_m04():

    thr = 1.0e-6

    return_result = {
        "dipole": np.array([-0.6691734918153486, +0.2578255474279201, -0.1528081650463131]),
        "mulliken_charges": np.array(
            [
                -0.0370404774665983,
                +0.0015354242523749,
                -0.0660191018009798,
                -0.3208773810551390,
                +0.5826454505938811,
                -0.0279475114762505,
                +0.0559054711521047,
                +0.1201622158540841,
                -0.0587401615749275,
                -0.0832435708476980,
                -0.2022582682225247,
                +0.3853632328414913,
                +0.1003824016037064,
                -0.5199127882822523,
                +0.0476813673425397,
                +0.0223636970861876,
            ]
        ),
    }

    atomic_input = qcel.models.AtomicInput(
        molecule=qcng.get_molecule("mindless-04"),
        model={"method": "GFN1-xTB"},
        driver="properties",
    )

    atomic_result = qcng.compute(atomic_input, "xtb")

    assert atomic_result.success
    assert pytest.approx(atomic_result.return_result["dipole"], abs=thr) == return_result["dipole"]
    assert pytest.approx(atomic_result.return_result["mulliken_charges"], abs=thr) == return_result["mulliken_charges"]
    assert "mayer_indices" in atomic_result.return_result


@using("xtb")
def test_xtb_task_gfn1xtb_m05():

    thr = 1.0e-8

    return_result = -29.038403257613453

    atomic_input = qcel.models.AtomicInput(
        molecule=qcng.get_molecule("mindless-05"),
        model={"method": "GFN1-xTB"},
        driver="energy",
    )

    atomic_result = qcng.compute(atomic_input, "xtb")

    assert atomic_result.success
    assert pytest.approx(atomic_result.return_result, abs=thr) == return_result
    assert "xtb" in atomic_result.extras


@using("xtb")
def test_xtb_task_gfn2xtb_m01():

    thr = 1.0e-7

    return_result = np.array(
        [
            [+1.7603980711827444e-03, -1.0892843650556535e-03, -1.4000925937447548e-02],
            [+7.2828551455461962e-03, -3.7216708563619287e-04, -1.7598863191610557e-03],
            [-4.3636138516365761e-04, +3.0781059500885385e-02, -2.4027777032994326e-03],
            [-4.9615553784293671e-03, +2.4184700438294312e-03, -1.8273143544031206e-03],
            [-3.1941473089998765e-03, -2.0177319546503510e-03, +6.0230398095918502e-03],
            [+3.5807344500847016e-03, -2.3730921963802166e-03, -2.0271824139507818e-03],
            [+1.0323930350176554e-03, +3.3301666520180575e-04, +2.6197885276710883e-03],
            [+6.4140580960895879e-03, -1.0424376542844720e-02, +1.7654830175669090e-02],
            [+2.1841001857333085e-03, +8.5100245114072097e-03, -6.7175792414255900e-03],
            [+8.1719403581595913e-03, +7.3797085798300455e-03, -7.5901731542113542e-03],
            [+1.1557950311395539e-03, -4.1907703068314668e-04, +1.7047910821890132e-03],
            [-4.5861988817335895e-03, -2.1269884595520816e-02, -1.4002736618139961e-02],
            [+1.1542155825459586e-04, +8.6721397831818186e-03, +3.9953426868203148e-03],
            [-1.8552117264088834e-03, -2.0185103100139296e-02, +1.3339454522070686e-02],
            [-2.9643338287258375e-03, -4.2193920947953519e-03, +6.6460843884674770e-06],
            [-1.3699887421746685e-02, +4.2756898813700889e-03, +4.9846828536382450e-03],
        ]
    )

    atomic_input = qcel.models.AtomicInput(
        molecule=qcng.get_molecule("mindless-01"),
        model={"method": "GFN2-xTB"},
        driver="gradient",
    )

    atomic_result = qcng.compute(atomic_input, "xtb")

    assert atomic_result.success
    assert pytest.approx(atomic_result.return_result, abs=thr) == return_result


@using("xtb")
def test_xtb_task_gfn2xtb_m02():

    thr = 1.0e-8

    return_result = np.array(
        [
            [+3.8555105053060808e-03, -2.2109578613126895e-03, -9.1601564191558012e-03],
            [+1.0970312247719392e-03, +1.3996306191282932e-02, -5.3472776842492294e-03],
            [-8.0585598914535359e-03, -2.3489450628686069e-03, +1.9653181376190359e-02],
            [-4.4938919015412053e-03, -4.5670953184344815e-04, -1.5746266791920045e-02],
            [+4.6136585412341665e-03, -8.3294935486403383e-03, -6.7825684743429791e-03],
            [-4.0644714892632226e-05, -7.6816568320453270e-03, +7.7872915896373372e-03],
            [+7.8804497759029323e-04, -4.9765353392257116e-03, +5.9586247342918366e-03],
            [+5.3161448621741650e-04, +5.1936947267340041e-03, +8.0564751525841038e-03],
            [-9.2290119389112895e-03, -1.5299670210190750e-02, -2.1407890260109627e-02],
            [+6.9445890241134674e-03, -3.7434087434538437e-04, +5.8421593307521448e-03],
            [+3.5614525125536099e-03, +7.6527534931603086e-03, +4.8040792443537744e-03],
            [-1.1121288887907211e-02, -6.4733203973410449e-03, +7.1241182667347342e-04],
            [+5.0828298212722894e-03, +1.2142928028482842e-02, +9.1892848637567797e-04],
            [+1.5094869632378746e-03, -2.3321110942833823e-03, +4.8244261495168301e-03],
            [+9.2444414237502198e-03, +1.2118975270529935e-03, -1.2700207457819403e-03],
            [-4.2852621453414374e-03, +1.0286160785383667e-02, +1.1566024851840597e-03],
        ]
    )

    atomic_input = qcel.models.AtomicInput(
        molecule=qcng.get_molecule("mindless-02"),
        model={"method": "GFN2-xTB"},
        driver="gradient",
        keywords={
            "accuracy": 0.1,
            "electronic_temperature": 500.0,
        },
    )

    atomic_result = qcng.compute(atomic_input, "xtb")

    assert atomic_result.success
    assert pytest.approx(atomic_result.return_result, abs=thr) == return_result


@using("xtb")
def test_xtb_task_gfn2xtb_m03():

    thr = 1.0e-8

    return_result = np.array(
        [
            [+0.0070908037236785, -0.0040563410767159, -0.0079330935995817],
            [+0.0064835107069477, -0.0016161473196322, +0.0038446333826145],
            [+0.0015491578122356, +0.0017538999500591, +0.0023145617382237],
            [-0.0072614197876866, -0.0045215532254949, -0.0067907981445658],
            [-0.0057299589336137, +0.0153458988830468, -0.0096315407254071],
            [-0.0062312740620834, +0.0050171730261916, -0.0018412044311888],
            [-0.0052673804035766, -0.0144351034054319, +0.0059596174012332],
            [-0.0008246247806272, -0.0004517884963982, +0.0019491982249825],
            [+0.0028511549191436, -0.0014262045374241, -0.0002516625700107],
            [+0.0055212287468860, -0.0022357394178147, +0.0016154868861744],
            [+0.0075014451430198, +0.0015965713047808, +0.0018304027741802],
            [-0.0086660517089299, +0.0117934694403132, +0.0058233201373099],
            [-0.0046887504989642, -0.0043855041726185, -0.0061676141160093],
            [-0.0022854923487687, -0.0025230789065119, +0.0037380292185573],
            [+0.0059302325276057, -0.0013873942542295, +0.0045321840851373],
            [+0.0040274189447333, +0.0015318422078802, +0.0010084797383503],
        ]
    )

    atomic_input = qcel.models.AtomicInput(
        molecule=qcng.get_molecule("mindless-03"),
        model={"method": "GFN2-xTB"},
        driver="gradient",
        keywords={
            "solvent": "chcl3",
        },
    )

    atomic_result = qcng.compute(atomic_input, "xtb")

    assert atomic_result.success
    assert pytest.approx(atomic_result.return_result, abs=thr) == return_result


@using("xtb")
def test_xtb_task_gfn2xtb_m04():

    thr = 1.0e-6

    return_result = {
        "dipole": np.array([-0.6566343439892861, +0.3553407721120402, -0.0768293825933504]),
        "mulliken_charges": np.array(
            [
                -0.0899891425563753,
                +0.0942167152764677,
                -0.1493904808718346,
                -0.2991140763033323,
                +0.4855285445087461,
                -0.0683158235225622,
                +0.0149989651434585,
                +0.2793618023805854,
                -0.1240704542128441,
                -0.0936742758503764,
                -0.2190620093200218,
                +0.2145448240515186,
                +0.3061608537234609,
                -0.3861071843835874,
                -0.0015140477313429,
                +0.0364257896680373,
            ]
        ),
    }

    atomic_input = qcel.models.AtomicInput(
        molecule=qcng.get_molecule("mindless-04"),
        model={"method": "GFN2-xTB"},
        driver="properties",
    )

    atomic_result = qcng.compute(atomic_input, "xtb")

    assert atomic_result.success
    assert pytest.approx(atomic_result.return_result["dipole"], abs=thr) == return_result["dipole"]
    assert pytest.approx(atomic_result.return_result["mulliken_charges"], abs=thr) == return_result["mulliken_charges"]
    assert "mayer_indices" in atomic_result.return_result


@using("xtb")
def test_xtb_task_gfn2xtb_m05():

    thr = 1.0e-8

    return_result = -27.73598761779656

    atomic_input = qcel.models.AtomicInput(
        molecule=qcng.get_molecule("mindless-05"),
        model={"method": "GFN2-xTB"},
        driver="energy",
    )

    atomic_result = qcng.compute(atomic_input, "xtb")

    assert atomic_result.success
    assert pytest.approx(atomic_result.return_result, abs=thr) == return_result
    assert "xtb" in atomic_result.extras


@using("xtb")
def test_xtb_task_unknown_method():

    atomic_input = qcel.models.AtomicInput(
        molecule=qcng.get_molecule("water"),
        model={"method": "GFN-xTB"},
        driver="energy",
    )
    error = qcel.models.ComputeError(error_type="input_error", error_message="Invalid method GFN-xTB provided in model")

    atomic_result = qcng.compute(atomic_input, "xtb")

    assert not atomic_result.success
    assert atomic_result.error == error


@using("xtb")
def test_xtb_task_unsupported_driver():

    atomic_input = qcel.models.AtomicInput(
        molecule=qcng.get_molecule("water"),
        model={"method": "GFN2-xTB"},
        driver="hessian",
    )
    error = qcel.models.ComputeError(
        error_type="input_error", error_message="Calculation succeeded but invalid driver request provided"
    )

    atomic_result = qcng.compute(atomic_input, "xtb")

    assert not atomic_result.success
    assert atomic_result.error == error


@using("xtb")
def test_xtb_task_cold_fusion():

    atomic_input = qcel.models.AtomicInput(
        molecule={
            "symbols": ["Li", "Li", "Li", "Li"],
            "geometry": [
                [-1.58746019997201, +1.58746019997201, +1.58746019997201],
                [-1.58746019997201, +1.58746019997201, +1.58746019997201],
                [-1.58746019997201, -1.58746019997201, -1.58746019997201],
                [+1.58746019997201, +1.58746019997201, -1.58746019997201],
            ],
            "validated": True,  # Force a nuclear fusion input, to make xtb fail
        },
        model={"method": "GFN2-xTB"},
        driver="energy",
    )
    error = qcel.models.ComputeError(
        error_type="runtime_error",
        error_message="Setup of molecular structure failed:\n-1- xtb_api_newMolecule: Could not generate molecular structure",
    )

    atomic_result = qcng.compute(atomic_input, "xtb")

    assert not atomic_result.success
    assert atomic_result.error == error
