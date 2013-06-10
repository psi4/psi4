from grendel.gmath.tensor import Tensor
from grendel.gmath.vector import Vector
from grendel.util.strings import classname

def test_classname():
    assert classname(object) == "object"
    assert classname(Tensor) == "Tensor"
    assert classname("Torsion") == "Torsion"
    assert classname(Vector(1,2,3)) == "Vector"