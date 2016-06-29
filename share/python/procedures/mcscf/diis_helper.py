import psi4
import numpy as np

class DIIS_helper(object):
    """
    Simple DIIS class.
    """

    def __init__(self, max_vec=6):
        self.error = []
        self.vector = []
        self.max_vec = max_vec

    def add(self, matrix, error):
        #if len(self.error) > 1:
        #    if self.error[-1].shape[0] != error.size:
        #        raise Exception("Error vector size does not match previous vector.")
        #    if self.vector[-1].shape != matrix.shape:
        #        raise Exception("Vector shape does not match previous vector.")

        self.error.append(error.clone())
        self.vector.append(matrix.clone())

    def extrapolate(self):
        # Limit size of DIIS vector
        diis_count = len(self.vector)

        if diis_count == 0:
            raise Exception("DIIS: No previous vectors.")
        if diis_count == 1:
            return self.vector[0]

        if diis_count > self.max_vec:
            # Remove oldest vector
            del self.vector[0]
            del self.error[0]
            diis_count -= 1

        # Build error matrix B
        B = np.empty((diis_count + 1, diis_count + 1))
        B[-1, :] = 1
        B[:, -1] = 1
        B[-1, -1] = 0
        for num1, e1 in enumerate(self.error):
            B[num1, num1] = e1.vector_dot(e1)
            for num2, e2 in enumerate(self.error):
                if num2 >= num1: continue
                val = e1.vector_dot(e2)
                B[num1, num2] = B[num2, num1] = val

        # Build residual vector
        resid = np.zeros(diis_count + 1)
        resid[-1] = 1

        # Solve pulay equations
        ci = np.linalg.solve(B, resid)

        # Yea, yea this is unstable make it stable
        iszero = np.any(np.diag(B)[:-1] <= 0.0)
        if iszero:
            S = np.ones((diis_count + 1))
        else:
            S = np.ones((diis_count + 1))
            S[:-1] = np.diag(B)[:-1]
            S = S ** -0.5
            S[-1] = 1

        # Then we gotta do a custom inverse
        B *= S[:, None] * S

        eigvals, eigvecs = np.linalg.eigh(B)
        maxval = np.max(np.abs(eigvals[[0, -1]])) * 1.e-12

        # If the relative is too small, zero it out
        eigvals[(np.abs(eigvals) < maxval)] = 0

        # Make sure we dont invert actual zeros!
        eigvals[np.abs(eigvals) > 1.e-16] = eigvals[np.abs(eigvals) > 1.e-16] ** -1

        invB = np.dot(eigvecs * eigvals, eigvecs.T)
        ci = np.dot(invB, resid) * S

        # combination of previous fock matrices
        V = psi4.Matrix("DIIS result", self.vector[0].rowdim(), self.vector[1].coldim())
        for num, c in enumerate(ci[:-1]):
            V.axpy(c, self.vector[num])

        return V
