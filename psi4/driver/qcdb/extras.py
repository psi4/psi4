try:
    from psi4.driver import constants
except ModuleNotFoundError:
    import qcelemental
    constants = qcelemental.PhysicalConstantsContext("CODATA2018")
