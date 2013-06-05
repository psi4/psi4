
import sys

def py_version_at_least(*nums):
    return sys.version_info >= nums

def py_version_before(*nums):
    return sys.version_info < nums

# No way to do Abstract class methods in Python 2; Python 3.2 and later should have @abstractclassmethod
if py_version_before(3, 2):
    def abstractclassmethod(f):
        cm = classmethod(f)
        setattr(f, '__isabstractmethod__', True)
        return cm
else:
    from abc import abstractclassmethod
