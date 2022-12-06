from addons import *

@ctest_labeler("misc")
def test_psithon2():
    ctest_runner(__file__, [
        "psiaux1/myccpvdzri.gbs",
        "psiaux1/S33.py",
        "psiaux1/myplugin1/doc.rst",
        "psiaux1/myplugin1/__init__.py",
        "psiaux1/myplugin1/inputalt.dat",
        "psiaux1/myplugin1/input.dat",
        "psiaux1/myplugin1/pymodule.py",
        "psiaux2/mysto3g.gbs",
        "psiaux2/S44.py",
    ])

