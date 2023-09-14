try:
    from psi4.driver import constants
except ModuleNotFoundError:
    import qcelemental as qcel
    constants = qcel.PhysicalConstantsContext("CODATA2014")
