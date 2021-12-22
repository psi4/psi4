from enum import Enum
from psi4 import core
from psi4.driver import psifiles as psif

import numpy as np

class RemovalPolicy(Enum):
    LargestError = 1
    OldestAdded = 2

class StoragePolicy(Enum):
    InCore = 1
    OnDisk = 2

def axpy(y, alpha, x):
    if isinstance(y, (core.Matrix, core.Vector)):
        y.axpy(alpha, x)
    elif isinstance(y, (core.dpdbuf4, core.dpdfile2)):
        y.axpy_matrix(x, alpha)
    else:
        raise TypeError("Unrecognized object type for DIIS.")

def template_helper(template, *args):
    if template:
        raise Exception("Template already set.")

    for arg in args:
        if isinstance(arg, core.Vector):
            template.append([arg.dimpi()])
        elif isinstance(arg, (core.Matrix, core.dpdfile2, core.dpdbuf4)):
            template.append([arg.rowdim(), arg.coldim()])
        else:
            raise TypeError("Unrecognized object type for DIIS.")

class DIIS:

    def __init__(self, max_vecs: int, name: str, removal_policy = RemovalPolicy.LargestError, storage_policy = StoragePolicy.OnDisk):
        if not isinstance(removal_policy, RemovalPolicy):
            raise TypeError(f"removal_policy must be a RemovalPolicy, not a {type(removal_policy)}")
        if not isinstance(storage_policy, StoragePolicy):
            raise TypeError(f"stroage_policy must be a StoragePolicy, not a {type(storage_policy)}")
        self.max_vecs = max_vecs
        self.name = name
        self.removal_policy = removal_policy
        self.storage_policy = storage_policy
        self.R_template = []
        self.T_template = []
        self.reset_subspace()
        self.opened_libdiis = False
        if self.storage_policy == StoragePolicy.OnDisk:
            psio = core.IO.shared_object()
            if not psio.open_check(psif.PSIF_LIBDIIS):
                psio.open(psif.PSIF_LIBDIIS, 1) # 1 = PSIO_OPEN_OLD
                self.opened_libdiis = True

    def __del__(self):
        if self.opened_libdiis:
            psio = core.IO.shared_object()
            if psio.open_check(psif.PSIF_LIBDIIS):
                psio.close(psif.PSIF_LIBDIIS, 1) # 1 = KEEP

    def reset_subspace(self):
        self.stored_vectors = [] # list[tuple[R entry, T entry]]
        self.last_added = -1
        self.cached_dot_products = dict()

    def copier(self, x, new_name):
        if isinstance(x, (core.Matrix, core.Vector)):
            copy = x.clone()
        elif isinstance(x, (core.dpdbuf4, core.dpdfile2)):
            copy = core.Matrix(x)
        else:
            raise TypeError("Unrecognized object type for DIIS.")

        copy.name = new_name

        if self.storage_policy == StoragePolicy.OnDisk:
            psio = core.IO.shared_object()
            if isinstance(x, core.Vector):
                copy.save(psio, psif.PSIF_LIBDIIS)
            else:
                copy.save(psio, psif.PSIF_LIBDIIS, core.SaveType.SubBlocks)
            copy = None

        return copy

    def get_name(self, name, entry_num, item_num):
        return f"{self.name}: {name} Entry {entry_num}, Item {item_num}"

    def get_dot_product(self, i, j):
        key = frozenset([i, j])
        try:
            return self.cached_dot_products[key]
        except KeyError:
            if self.storage_policy == StoragePolicy.InCore:
                Ri = self.stored_vectors[i][0]
                Rj = self.stored_vectors[j][0]
                dot_product = sum(Rix.vector_dot(Rjx) for Rix, Rjx in zip(Ri, Rj))
            elif self.storage_policy == StoragePolicy.OnDisk:
                dot_product = 0
                psio = core.IO.shared_object()
                for x, entry_dims in enumerate(self.R_template):
                    if len(entry_dims) == 2:
                        Rix = core.Matrix(self.get_name("R", i, x), *entry_dims)
                        Rjx = core.Matrix(self.get_name("R", j, x), *entry_dims)
                        Rix.load(psio, psif.PSIF_LIBDIIS, core.SaveType.SubBlocks)
                        Rjx.load(psio, psif.PSIF_LIBDIIS, core.SaveType.SubBlocks)
                    elif len(entry_dims) == 1:
                        Rix = core.Vector(self.get_name("R", i, x), *entry_dims)
                        Rjx = core.Vector(self.get_name("R", j, x), *entry_dims)
                        Rix.load(psio, psif.PSIF_LIBDIIS)
                        Rjx.load(psio, psif.PSIF_LIBDIIS)
                    else:
                        raise Exception("R_template may only have 1 or 2 dimensions. This is a bug: contact developers.")
                    dot_product += Rix.vector_dot(Rjx)
            else:
                raise Exception(f"StoragePolicy {self.storage_policy} not recognized. This is a bug: contact developers.")

            self.cached_dot_products[key] = dot_product
            return dot_product


    def set_error_vector_size(self, *args):
        template_helper(self.R_template, *args)

    def set_vector_size(self, *args):
        template_helper(self.T_template, *args)

    def build_entry(self, entry, target_index):
        if len(self.R_template) + len(self.T_template) != len(entry):
            raise Exception(f"Cannot build {len(self.R_template)} residuals and {len(self.T_template)} amplitudes from {len(entries)} items.")
        R = entry[:len(self.R_template)]
        T = entry[len(self.R_template):]

        return [
                [self.copier(Ri, self.get_name("R", target_index, i)) for i, Ri in enumerate(R)],
                [self.copier(Ti, self.get_name("T", target_index, i)) for i, Ti in enumerate(T)]
        ]

    def add_entry(self, *args):
        if self.max_vecs == 0:
            return False

        self.last_added += 1
        self.last_added % self.max_vecs

        if len(self.stored_vectors) >= self.max_vecs:
            if self.removal_policy == RemovalPolicy.OldestAdded:
                target_index = self.last_added
            elif self.removal_policy == RemovalPolicy.LargestError:
                target_index = np.argmax([self.get_dot_product(i, i) for i in range(len(self.stored_vectors))])
            else:
                raise Exception(f"RemovalPolicy {self.removal_policy} not recognized. This is a bug: contact developers.")
            # Purge imminently-outdated values from cache.
            self.cached_dot_products = {key: val for key, val in self.cached_dot_products.items() if target_index not in key}
            # Set the new entry.
            self.stored_vectors[target_index] = self.build_entry(args, target_index)
        else:
            self.stored_vectors.append(self.build_entry(args, self.last_added))

        return True

    def extrapolate(self, *args):

        dim = len(self.stored_vectors) + 1
        B = np.zeros((dim, dim))
        for i in range(len(self.stored_vectors)):
            for j in range(len(self.stored_vectors)):
                B[i, j] = self.get_dot_product(i, j)
        B[-1, :-1] = B[:-1, -1] = -1

        rhs = np.zeros((dim))
        rhs[-1] = -1

        # Trick to improve numerical conditioning.
        # Instead of solving B c = r, we solve D B D^-1 D c = D r, using
        # D r = r. D is the diagonals ^ -1/2 matrix.
        # This improves the conditioning of the problem.
        diagonals = B.diagonal().copy()
        diagonals[-1] = 1
        if np.all(diagonals > 0):
            diagonals = diagonals ** (- 0.5)
            B = np.einsum("i,ij,j -> ij", diagonals, B, diagonals)
            coeffs = np.linalg.lstsq(B, rhs, rcond=None)[0][:-1] * diagonals[:-1]
        else:
            coeffs = np.linalg.lstsq(B, rhs, rcond=None)[0][:-1]

        for j, Tj in enumerate(args):
            Tj.zero()
            if self.storage_policy == StoragePolicy.InCore:
                for ci, (_, Ti) in zip(coeffs, self.stored_vectors):
                    axpy(Tj, ci, Ti[j])
            elif self.storage_policy == StoragePolicy.OnDisk:
                for i, ci in enumerate(coeffs):
                    psio = core.IO.shared_object()
                    if isinstance(Tj, core.Vector):
                        Tij = core.Vector(self.get_name("T", i, j), *self.T_template[j])
                        Tij.load(psio, psif.PSIF_LIBDIIS)
                    elif isinstance(Tj, (core.Matrix, core.dpdfile2, core.dpdbuf4)):
                        Tij = core.Matrix(self.get_name("T", i, j), *self.T_template[j])
                        Tij.load(psio, psif.PSIF_LIBDIIS, core.SaveType.SubBlocks)
                    else:
                        raise TypeError("Unrecognized object type for DIIS.")
                    axpy(Tj, ci, Tij)
            else:
                raise Exception(f"StoragePolicy {self.storage_policy} not recognized. This is a bug: contact developers.")

        return True

    def delete_diis_file(self):
        psio = core.IO.shared_object()
        if not psio.open_check(psif.PSIF_LIBDIIS):
            psio.open(psif.PSIF_LIBDIIS, 1) # 1 = PSIO_OPEN_OLD
        psio.close(psif.PSIF_LIBDIIS, 0) # 0 = DELETE

