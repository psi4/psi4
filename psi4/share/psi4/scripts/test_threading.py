#!/usr/bin/env python

import time

# good
#import numpy as np
#import psi4

# bad
import psi4
import numpy as np


def test_threaded_blas(args):
    threads = int(args.nthread)

    times = {}

    size = [200, 500, 2000, 4000]
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
            if sz == 4000:
                times["p4-n{}".format(th)] = retp4
                times["np-n{}".format(th)] = retnp
                assert psi4.get_num_threads() == th

    rat1 = times["np-n" + str(threads[-1])] / times["p4-n" + str(threads[-1])]
    rat2 = times["p4-n" + str(threads[0])] / times["p4-n" + str(threads[-1])]
    print("  NumPy@n%d : Psi4@n%d ratio (want ~1): %.2f" % (threads[-1], threads[-1], rat1))
    print("   Psi4@n%d : Psi4@n%d ratio (want ~%d): %.2f" % (threads[0], threads[-1], threads[-1], rat2))


def test_psithon(args):
    import subprocess

    inputdat = """
molecule dimer {
0 1
C   0.000000  -0.667578  -2.124659
C   0.000000   0.667578  -2.124659
H   0.923621  -1.232253  -2.126185
H  -0.923621  -1.232253  -2.126185
H  -0.923621   1.232253  -2.126185
H   0.923621   1.232253  -2.126185
--
0 1
C   0.000000   0.000000   2.900503
C   0.000000   0.000000   1.693240
H   0.000000   0.000000   0.627352
H   0.000000   0.000000   3.963929
units angstrom
}

set {
    BASIS jun-cc-pVQZ
    SCF_TYPE DF
    FREEZE_CORE True
}

energy('sapt0')

compare_values(85.189064196429, dimer.nuclear_repulsion_energy(), 9,        "Nuclear Repulsion Energy")
compare_values(-0.00343131388,  psi4.get_variable("SAPT ELST ENERGY"), 6,   "SAPT0 elst")
compare_values( 0.00368418851,  psi4.get_variable("SAPT EXCH ENERGY"), 6,   "SAPT0 exch")
compare_values(-0.00094288048,  psi4.get_variable("SAPT IND ENERGY"), 6,    "SAPT0 ind")
compare_values(-0.00231901067,  psi4.get_variable("SAPT DISP ENERGY"), 6,   "SAPT0 disp")
compare_values(-0.00300901652,  psi4.get_variable("SAPT0 TOTAL ENERGY"), 6, "SAPT0")
"""

    tfn = '_thread_test_input_psi4_yo'
    with open(tfn + '.in', 'w') as fp:
        fp.write(inputdat)

    times = {}
    threads = [1, int(args.nthread)]

    for nt in threads:
        t0 = time.time()
        cmd = """psi4 -i {tfn}.in -o {tfn}_n{nt}.out -n{nt}""".format(tfn=tfn, nt=nt)
        print('Running {} ...'.format(cmd))
        subprocess.call(cmd.split())

        t1 = time.time()
        times[nt] = t1 - t0
        print("Time for threads %2d: Psi4: %12.6f" % (nt, times[nt]))

    rat1 = times[threads[0]] / times[threads[-1]]
    print("   Psi4@n%d : Psi4@n%d ratio (want ~%d): %.2f" % (threads[0], threads[-1], threads[-1], rat1))


def print_math_ldd():
    import sys
    import subprocess

    import psi4
    p4core = psi4.__file__[:-11] + 'core.so'

    if sys.platform.startswith('linux'):
        cmd = """ldd -v {} | grep -e ':' -e 'mkl' -e 'openblas' -e 'iomp5' -e 'gomp'""".format(p4core)
        print('Running {} ...'.format(cmd))
        subprocess.call(cmd, shell=True)
    else:
        print('Not available w/o `ldd`')


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description="Psi4 threading tester. `psi4` and `import psi4` expected in *PATHs")
    parser.add_argument("-n", "--nthread", default=4,
                        help="Number of threads to use. Psi4 disregards OMP_NUM_THREADS/MKL_NUM_THREADS.")

    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument('--psithon', action='store_true',
                        help="""Test Psi4 (PsiAPI) vs NumPy in threaded matrix multiply.""")
    group.add_argument('--psiapi', action='store_true',
                        help="""Test Psi4 (PSIthon) in threaded SAPT calc.""")
    group.add_argument('--ldd', action='store_true',
                        help="""Run ldd to examine BLAS and OMP linking of psi4/core.so.""")

    args, unknown = parser.parse_known_args()

    if args.psiapi:
        test_threaded_blas(args)
    elif args.psithon:
        test_psithon(args)
    elif args.ldd:
        print_math_ldd()
    else:
        test_threaded_blas(args)
        test_psithon(args)
        print_math_ldd()
