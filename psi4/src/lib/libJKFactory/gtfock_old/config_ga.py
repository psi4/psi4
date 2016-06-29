import os
import sys

# global arrays
ga_url = "http://hpc.pnl.gov/globalarrays/download/ga-5-3.tgz"
ga_dir = "ga-5-3"
ga_tgz = ga_dir + ".tgz"

#blas
blas_inc = "/opt/apps/intel/13/composer_xe_2013_sp1.1.106/mkl/include" 
blas_libdir = "/opt/apps/intel/13/composer_xe_2013_sp1.1.106/mkl/lib/intel64"
blas_libs = "-lmkl_core -lmkl_intel_lp64 -lmkl_intel_thread -liomp5"
lapack_inc = blas_inc
lapack_libdir = blas_libdir
lapack_libs = blas_libs
scalapack_inc = blas_inc
scalapack_libdir = blas_libdir
scalapack_libs = "-lmkl_scalapack_lp64"

#mpi
mpi_inc = "/opt/apps/intel14/mvapich2/2.0b/include"
mpi_libdir = "/opt/apps/intel14/mvapich2/2.0b/lib" 
mpi_libs = ""

#openib
openib_device = "mlx4_0"
openib_inc = "/opt/ofed/include"
openib_libdir = "/opt/ofed/lib64"
openib_libs = "-libumad -libverbs"

# initialize variales
str = os.popen("pwd").read()
pwdir = str.strip("\n")
install_dir = pwdir + "/GAlib-armci"
config_cmd = "./configure"
config_cmd += " --prefix=" + install_dir
config_cmd += " --with-blas=\"" + "-I" + blas_inc + " -L" + blas_libdir + " " + blas_libs + "\""
config_cmd += " --with-lapack=\"" + "-I" + lapack_inc + " -L" + lapack_libdir + " " + lapack_libs + "\""
config_cmd += " --with-scalapack=\"" + "-I" + scalapack_inc + " -L" + scalapack_libdir + " " + scalapack_libs + "\""
config_cmd += " --with-mpi=\"" + "-I" + mpi_inc + " -L" + mpi_libdir + " " + mpi_libs + "\""

with_openib = False
download_ga = False
# parse arguments
for str in sys.argv:
    if str == sys.argv[0]:
        print "Config global arrays"
    elif str == "download":
        download_ga = True
    elif str == "openib":
        with_openib = True
    else:
        print "Usage: " + sys.argv[0] + " <openib> <download>"
        sys.exit(0)       

if download_ga:
    os.system("wget " + ga_url + " -O " + ga_tgz)
    os.system("tar xzf " + ga_tgz)
    os.system ("rm -f " + ga_tgz)

if with_openib:
    config_cmd += " --with-openib=\"" + "-I" + openib_inc + " -L" + openib_libdir + " " + openib_libs + "\""

print config_cmd

os.chdir(ga_dir)
os.system("export ARMCI_OPENIB_DEVICE=\"" + openib_device + "\"")
os.system(config_cmd)
os.system("make CC=icc FC=ifort MPICC=mpicc")
os.system("make install")
