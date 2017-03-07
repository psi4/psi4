import time
import pytest
import numpy as np
import multiprocessing
import psi4


# Test below is fine on its own but erratic through pytest. Most likely
#   to succeed as first test collected, so here it lies.
@pytest.mark.xfail(True, reason='threading treatment suspect', run=True)
def test_threaded_blas():
    threads = multiprocessing.cpu_count()
    threads = int(threads / 2)

    times = {}

    size = [200, 500, 2000, 5000]
    threads = [1, threads]

    for th in threads:
        psi4.set_num_threads(th)

        for sz in size:
            nruns = max(1, int(1.e10 / (sz ** 3)))

            a = psi4.core.Matrix(sz, sz)
            b = psi4.core.Matrix(sz, sz)
            c = psi4.core.Matrix(sz, sz)

            tp4 = time.time()
            for n in range(nruns):
                c.gemm(False, False, 1.0, a, b, 0.0)

            retp4 = (time.time() - tp4) / nruns

            tnp = time.time()
            for n in range(nruns):
                np.dot(a, b, out=np.asarray(c))

            retnp = (time.time() - tnp) / nruns
            print("Time for threads %2d, size %5d: Psi4: %12.6f  NumPy: %12.6f" % (th, sz, retp4, retnp))
            if sz == 5000:
                times["p4-n{}".format(th)] = retp4
                times["np-n{}".format(th)] = retnp
                assert psi4.get_num_threads() == th

    rat1 = times["np-n" + str(threads[-1])] / times["p4-n" + str(threads[-1])]
    rat2 = times["p4-n" + str(threads[0])] / times["p4-n" + str(threads[-1])]
    print("  NumPy@n%d : Psi4@n%d ratio (want ~1): %.2f" % (threads[-1], threads[-1], rat1))
    print("   Psi4@n%d : Psi4@n%d ratio (want ~%d): %.2f" % (threads[0], threads[-1], threads[-1], rat2))
    assert pytest.approx(rat1, 0.2) == 1.0
    assert pytest.approx(rat2, 0.8) == threads[-1]
