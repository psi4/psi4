import pytest

import qcdb

co_dominant = (59, 27, 'Co', 58.933195048, True, '')
co_dominant_mine = (59, 27, 'Co', 58.933195048, True, '_mine')
co_dominant_shortmass = (59, 27, 'Co', 58.933, True, '')
co60 = (60, 27, 'Co', 59.933817059, True, '')
co60ghost = (60, 27, 'Co', 59.933817059, False, '')
co_unspecified = (-1, 27, 'Co', 60.6, True, '')


def test_reconcile_nucleus_1():
    assert co_dominant == qcdb.molparse.reconcile_nucleus(E='co')


def test_reconcile_nucleus_2():
    assert co_dominant == qcdb.molparse.reconcile_nucleus(Z=27)


def test_reconcile_nucleus_3():
    assert co_dominant == qcdb.molparse.reconcile_nucleus(A=59, Z=27)


def test_reconcile_nucleus_4():
    assert co_dominant == qcdb.molparse.reconcile_nucleus(E='cO', mass=58.933195048)


def test_reconcile_nucleus_5():
    assert co_dominant == qcdb.molparse.reconcile_nucleus(A=59, Z=27, E='CO')


def test_reconcile_nucleus_6():
    assert co_dominant == qcdb.molparse.reconcile_nucleus(A=59, E='cO', mass=58.933195048)


def test_reconcile_nucleus_7():
    assert co_dominant == qcdb.molparse.reconcile_nucleus(label='co')


def test_reconcile_nucleus_8():
    assert co_dominant == qcdb.molparse.reconcile_nucleus(label='59co')


def test_reconcile_nucleus_9():
    assert co_dominant == qcdb.molparse.reconcile_nucleus(label='co@58.933195048')


def test_reconcile_nucleus_10():
    assert co_dominant == qcdb.molparse.reconcile_nucleus(
        A=59, Z=27, E='cO', mass=58.933195048, label='co@58.933195048')


def test_reconcile_nucleus_11():
    assert co_dominant == qcdb.molparse.reconcile_nucleus(
        A=59, Z=27, E='cO', mass=58.933195048, label='27@58.933195048')


def test_reconcile_nucleus_12():
    assert co_dominant == qcdb.molparse.reconcile_nucleus(label='27')


def test_reconcile_nucleus_13():
    assert co_dominant_mine == qcdb.molparse.reconcile_nucleus(label='co_miNe')


def test_reconcile_nucleus_14():
    assert co_dominant_mine == qcdb.molparse.reconcile_nucleus(label='co_mIne@58.933195048')


def test_reconcile_nucleus_15():
    assert co_dominant_shortmass == qcdb.molparse.reconcile_nucleus(E='cO', mass=58.933)


def test_reconcile_nucleus_16():
    assert co_dominant_shortmass == qcdb.molparse.reconcile_nucleus(label='cO@58.933')


def test_reconcile_nucleus_17():
    with pytest.raises(AssertionError):
        assert co_dominant_shortmass == qcdb.molparse.reconcile_nucleus(E='cO', mass=58.933, mtol=1.e-4)


def test_reconcile_nucleus_18():
    with pytest.raises(AssertionError):
        assert co_dominant_shortmass == qcdb.molparse.reconcile_nucleus(label='27@58.933', mtol=1.e-4)


def test_reconcile_nucleus_19():
    assert co60 == qcdb.molparse.reconcile_nucleus(E='Co', A=60)


def test_reconcile_nucleus_20():
    assert co60 == qcdb.molparse.reconcile_nucleus(Z=27, A=60, real=True)


def test_reconcile_nucleus_21():
    assert co60 == qcdb.molparse.reconcile_nucleus(E='Co', A=60)


def test_reconcile_nucleus_22():
    assert co60 == qcdb.molparse.reconcile_nucleus(Z=27, mass=59.933817059)


def test_reconcile_nucleus_23():
    assert co60 == qcdb.molparse.reconcile_nucleus(A=60, Z=27, mass=59.933817059)


def test_reconcile_nucleus_24():
    assert co60 == qcdb.molparse.reconcile_nucleus(label='60Co')


def test_reconcile_nucleus_25():
    assert co60 == qcdb.molparse.reconcile_nucleus(label='27', mass=59.933817059)


def test_reconcile_nucleus_26():
    assert co60 == qcdb.molparse.reconcile_nucleus(label='Co', mass=59.933817059)


def test_reconcile_nucleus_27():
    assert co60 == qcdb.molparse.reconcile_nucleus(A=60, label='Co')


def test_reconcile_nucleus_28():
    assert co60ghost == qcdb.molparse.reconcile_nucleus(E='Co', A=60, real=False)


def test_reconcile_nucleus_29():
    assert co60ghost == qcdb.molparse.reconcile_nucleus(A=60, Z=27, mass=59.933817059, real=0)


def test_reconcile_nucleus_30():
    assert co60ghost == qcdb.molparse.reconcile_nucleus(label='@60Co')


def test_reconcile_nucleus_31():
    assert co60ghost == qcdb.molparse.reconcile_nucleus(label='Gh(27)', mass=59.933817059)


def test_reconcile_nucleus_32():
    assert co60ghost == qcdb.molparse.reconcile_nucleus(label='@Co', mass=59.933817059)


def test_reconcile_nucleus_33():
    assert co60ghost == qcdb.molparse.reconcile_nucleus(A=60, label='Gh(Co)')


def test_reconcile_nucleus_34():
    assert co_unspecified == qcdb.molparse.reconcile_nucleus(mass=60.6, Z=27)


def test_reconcile_nucleus_35():
    assert co_unspecified == qcdb.molparse.reconcile_nucleus(mass=60.6, E='Co')


def test_reconcile_nucleus_36():
    assert co_unspecified == qcdb.molparse.reconcile_nucleus(mass=60.6, label='27')


def test_reconcile_nucleus_37():
    assert co_unspecified == qcdb.molparse.reconcile_nucleus(label='Co@60.6')


def test_reconcile_nucleus_38():
    with pytest.raises(qcdb.ValidationError):
        assert co_unspecified == qcdb.molparse.reconcile_nucleus(mass=60.6, Z=27, A=61)


def test_reconcile_nucleus_39():
    with pytest.raises(qcdb.ValidationError):
        qcdb.molparse.reconcile_nucleus(A=80, Z=27)


def test_reconcile_nucleus_40():
    with pytest.raises(qcdb.ValidationError):
        qcdb.molparse.reconcile_nucleus(Z=27, mass=200)


def test_reconcile_nucleus_41():
    qcdb.molparse.reconcile_nucleus(Z=27, mass=200, nonphysical=True)


def test_reconcile_nucleus_42():
    with pytest.raises(qcdb.ValidationError):
        qcdb.molparse.reconcile_nucleus(Z=27, mass=-200, nonphysical=True)


def test_reconcile_nucleus_43():
    with pytest.raises(qcdb.ValidationError):
        qcdb.molparse.reconcile_nucleus(Z=-27, mass=200, nonphysical=True)


def test_reconcile_nucleus_44():
    with pytest.raises(qcdb.ValidationError):
        qcdb.molparse.reconcile_nucleus(Z=1, label='he')


def test_reconcile_nucleus_45():
    with pytest.raises(qcdb.ValidationError):
        qcdb.molparse.reconcile_nucleus(A=4, label='3he')


def test_reconcile_nucleus_46():
    with pytest.raises(qcdb.ValidationError):
        qcdb.molparse.reconcile_nucleus(label='@U', real=True)


def test_reconcile_nucleus_47():
    with pytest.raises(qcdb.ValidationError):
        qcdb.molparse.reconcile_nucleus(label='U', real=False)
