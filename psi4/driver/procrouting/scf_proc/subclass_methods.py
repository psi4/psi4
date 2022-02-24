import math

from psi4 import core
from ..diis import DIIS, StoragePolicy, RemovalPolicy

def diis_engine_helper(self):
    engines = set()
    if core.get_option('SCF', 'DIIS'):
        engines.add('diis')
    restricted_open = self.same_a_b_orbs() and not self.same_a_b_dens()
    if not restricted_open:
        aediis = core.get_option('SCF', 'SCF_INITIAL_ACCELERATOR')
        if aediis != "NONE":
            engines.add(aediis.lower())
    return engines

def _RHF_orbital_gradient(self, save_fock: bool, max_diis_vectors: int) -> float:
    gradient = self.form_FDSmSDF(self.Fa(), self.Da())

    if save_fock:
        if not self.initialized_diis_manager_:
            storage_policy = StoragePolicy.InCore if self.scf_type() == "DIRECT" else StoragePolicy.OnDisk
            self.diis_manager_ = DIIS(max_diis_vectors, "HF DIIS vector", RemovalPolicy.LargestError, storage_policy, engines=diis_engine_helper(self))
            self.initialized_diis_manager_ = True

        entry = {"target": [self.Fa()]}
        if core.get_option('SCF', 'DIIS'):
            entry["error"] = [gradient]
        aediis = core.get_option('SCF', 'SCF_INITIAL_ACCELERATOR')
        if aediis != "NONE":
            entry["densities"] = [self.Da()]
            if aediis == "EDIIS":
                entry["energy"] = [self.compute_E()]
        self.diis_manager_.add_entry(entry)

    if self.options().get_bool("DIIS_RMS_ERROR"):
        return gradient.rms()
    else:
        return gradient.absmax()

def _UHF_orbital_gradient(self, save_fock: bool, max_diis_vectors: int) -> float:
    gradient_a = self.form_FDSmSDF(self.Fa(), self.Da())
    gradient_b = self.form_FDSmSDF(self.Fb(), self.Db())

    if save_fock:
        if not self.initialized_diis_manager_:
            self.diis_manager_ = DIIS(max_diis_vectors, "HF DIIS vector", RemovalPolicy.LargestError,
                                                          StoragePolicy.OnDisk, False, engines=diis_engine_helper(self))
            self.initialized_diis_manager_ = True

        entry = {"target": [self.Fa(), self.Fb()]}
        if core.get_option('SCF', 'DIIS'):
            entry["error"] = [gradient_a, gradient_b]
        aediis = core.get_option('SCF', 'SCF_INITIAL_ACCELERATOR')
        if aediis != "NONE":
            entry["densities"] = [self.Da(), self.Db()]
            if aediis == "EDIIS":
                entry["energy"] = [self.compute_E()]
        self.diis_manager_.add_entry(entry)

    if self.options().get_bool("DIIS_RMS_ERROR"):
        return math.sqrt(0.5 * (gradient_a.rms() ** 2 + gradient_b.rms() ** 2))
    else:
        return max(gradient_a.absmax(), gradient_b.absmax())

def _ROHF_orbital_gradient(self, save_fock: bool, max_diis_vectors: int) -> float:
    # Only the inact-act, inact-vir, and act-vir rotations are non-redundant
    dim_zero = core.Dimension(self.nirrep(), "Zero Dim")
    noccpi = self.doccpi() + self.soccpi()
    row_slice = core.Slice(dim_zero, noccpi)
    col_slice = core.Slice(self.doccpi(), self.nmopi())
    MOgradient = self.moFeff().get_block(row_slice, col_slice)

    # Zero the active-active block
    for h in range(MOgradient.nirrep()):
        socc = self.soccpi()[h]
        docc = self.doccpi()[h]

        MOgradient.nph[h][docc:docc+socc, 0:socc] = 0

    # Grab inact-act and act-vir orbs
    # Ct is (nmo x nmo), not the (nso x nmo) you would expect
    row_slice = core.Slice(dim_zero, self.nmopi())
    col_slice = core.Slice(dim_zero, noccpi)
    Cia = self.Ct().get_block(row_slice, col_slice)
    col_slice = core.Slice(self.doccpi(), self.nmopi())
    Cav = self.Ct().get_block(row_slice, col_slice)

    # Back transform MOgradient
    gradient = core.triplet(Cia, MOgradient, Cav, False, False, True)

    if save_fock:
        if not self.initialized_diis_manager_:
            self.diis_manager_ = DIIS(max_diis_vectors, "HF DIIS vector", RemovalPolicy.LargestError, StoragePolicy.OnDisk, engines=diis_engine_helper(self))
            self.diis_manager_.set_error_vector_size(gradient)
            self.diis_manager_.set_vector_size(self.soFeff())
            self.initialized_diis_manager_ = True

        self.diis_manager_.add_entry({"error": [gradient], "target": [self.soFeff()]})

    if self.options().get_bool("DIIS_RMS_ERROR"):
        return gradient.rms()
    else:
        return gradient.absmax()

core.RHF.compute_orbital_gradient = _RHF_orbital_gradient
core.UHF.compute_orbital_gradient = core.CUHF.compute_orbital_gradient = _UHF_orbital_gradient
core.ROHF.compute_orbital_gradient = _ROHF_orbital_gradient

def _RHF_diis(self, Dnorm):
    return self.diis_manager_.extrapolate(self.Fa(), Dnorm=Dnorm)

def _UHF_diis(self, Dnorm):
    return self.diis_manager_.extrapolate(self.Fa(), self.Fb(), Dnorm=Dnorm)

def _ROHF_diis(self, Dnorm):
    return self.diis_manager_.extrapolate(self.soFeff(), Dnorm=Dnorm)

core.RHF.diis = _RHF_diis
core.UHF.diis = core.CUHF.diis = _UHF_diis
core.ROHF.diis = _ROHF_diis
