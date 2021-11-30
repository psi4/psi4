# FROM mcr.microsoft.com/vscode/devcontainers/base:ubuntu-20.04
FROM mcr.microsoft.com/vscode/devcontainers/miniconda:3

ARG USERNAME=vscode
ARG USER_UID=1000
ARG USER_GID=${USER_UID}

# Avoid warnings by switching to noninteractive
ENV DEBIAN_FRONTEND=noninteractive
ENV APT_KEY_DONT_WARN_ON_DANGEROUS_USAGE=1

# install OS tools
RUN apt-get update -y \
    && apt-get install -y -o=Dpkg::Use-Pty=0 \
    build-essential \
    gdb \
    pkg-config \
    ca-certificates \
    gnupg \
    git \
    ninja-build \
    zsh \
    libarchive13 \
    liblapack-dev \
    #
    # [Optional] Update UID/GID if needed
    && if [ "${USER_GID}" != "1000"] || [ "${USER_UID}" != "1000" ]; then \
    groupmod --gid ${USER_GID} ${USERNAME} \
    && usermod --uid ${USER_UID} --gid ${USER_GID} ${USERNAME} \
    && chown -R ${USER_UID}:${USER_GID} /home/${USERNAME}; \
    fi

# install latest CMake
# RUN wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null | gpg --dearmor - | tee /etc/apt/trusted.gpg.d/kitware.gpg >/dev/null \
# && echo "deb https://apt.kitware.com/ubuntu/ focal main" >> /etc/apt/sources.list \
# && apt-get update \
# && apt-get -y install cmake

# ### Import Intel(R) apt repo public keys
# # add apt repo public key
# ARG url=https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2023.PUB
# ADD $url /
# RUN file=$(basename "$url") && \
#     apt-key add "$file" && \
#     rm "$file"

# # configure the repository
# ARG repo=https://apt.repos.intel.com/oneapi
# RUN echo "deb $repo all main" > /etc/apt/sources.list.d/oneAPI.list

# # disable cert check
# ARG disable_cert_check=
# RUN if [ "$disable_cert_check" ] ; then echo "Acquire::https::Verify-Peer \"false\";\nAcquire::https::Verify-Host \"false\";" > /etc/apt/apt.conf.d/99-disable-cert-check ; fi

# # install Intel(R) oneAPI Base Toolkit
# RUN apt-get update -y \
#     && apt-get install -y --no-install-recommends -o=Dpkg::Use-Pty=0 \
#     intel-basekit-getting-started \
#     intel-oneapi-advisor \
#     intel-oneapi-ccl-devel \
#     intel-oneapi-common-licensing \
#     intel-oneapi-common-vars \
#     intel-oneapi-compiler-dpcpp-cpp \
#     # intel-oneapi-dal-devel \
#     intel-oneapi-dev-utilities \
#     # intel-oneapi-dnnl-devel \
#     intel-oneapi-dpcpp-debugger \
#     intel-oneapi-ipp-devel \
#     # intel-oneapi-ippcp-devel \
#     intel-oneapi-libdpstd-devel \
#     intel-oneapi-mkl-devel \
#     # intel-oneapi-onevpl-devel \
#     intel-oneapi-python \
#     intel-oneapi-tbb-devel \
#     intel-oneapi-vtune \
#     --

# # install Intel(R) oneAPI HPC Toolkit
# RUN apt-get update -y \
#     && apt-get install -y --no-install-recommends -o=Dpkg::Use-Pty=0 \
#     intel-hpckit-getting-started \
#     intel-oneapi-clck \
#     # intel-oneapi-common-licensing \
#     # intel-oneapi-common-vars \
#     intel-oneapi-compiler-dpcpp-cpp-and-cpp-classic \
#     intel-oneapi-compiler-fortran \
#     # intel-oneapi-dev-utilities \
#     intel-oneapi-inspector \
#     intel-oneapi-itac \
#     # intel-oneapi-mpi-devel \
#     --

# # setvars.sh environment variables
# ENV ACL_BOARD_VENDOR_PATH='/opt/Intel/OpenCLFPGA/oneAPI/Boards'
# ENV ADVISOR_2021_DIR='/opt/intel/oneapi/advisor/2021.2.0'
# ENV APM='/opt/intel/oneapi/advisor/2021.2.0/perfmodels'
# ENV CCL_CONFIGURATION='cpu_gpu_dpcpp'
# ENV CCL_ROOT='/opt/intel/oneapi/ccl/2021.2.0'
# ENV CLASSPATH='/opt/intel/oneapi/mpi/2021.2.0//lib/mpi.jar:/opt/intel/oneapi/dal/2021.2.0/lib/onedal.jar'
# ENV CLCK_ROOT='/opt/intel/oneapi/clck/2021.2.0'
# ENV CMAKE_PREFIX_PATH='/opt/intel/oneapi/vpl:/opt/intel/oneapi/tbb/2021.2.0/env/..:/opt/intel/oneapi/dal/2021.2.0'
# ENV CONDA_DEFAULT_ENV='base'
# ENV CONDA_EXE='/opt/intel/oneapi/intelpython/latest/bin/conda'
# ENV CONDA_PREFIX='/opt/intel/oneapi/intelpython/latest'
# ENV CONDA_PROMPT_MODIFIER='(base) '
# ENV CONDA_PYTHON_EXE='/opt/intel/oneapi/intelpython/latest/bin/python'
# ENV CONDA_SHLVL='1'
# ENV CPATH='/opt/intel/oneapi/vpl/2021.2.2/include:/opt/intel/oneapi/tbb/2021.2.0/env/../include:/opt/intel/oneapi/mpi/2021.2.0//include:/opt/intel/oneapi/mkl/latest/include:/opt/intel/oneapi/ippcp/2021.2.0/include:/opt/intel/oneapi/ipp/2021.2.0/include:/opt/intel/oneapi/dpl/2021.2.0/linux/include:/opt/intel/oneapi/dnnl/2021.2.0/cpu_dpcpp_gpu_dpcpp/lib:/opt/intel/oneapi/dev-utilities/2021.2.0/include:/opt/intel/oneapi/dal/2021.2.0/include:/opt/intel/oneapi/compiler/2021.2.0/linux/include:/opt/intel/oneapi/ccl/2021.2.0/include/cpu_gpu_dpcpp'
# ENV CPLUS_INCLUDE_PATH='/opt/intel/oneapi/clck/2021.2.0/include'
# ENV DAALROOT='/opt/intel/oneapi/dal/2021.2.0'
# ENV DALROOT='/opt/intel/oneapi/dal/2021.2.0'
# ENV DAL_MAJOR_BINARY='1'
# ENV DAL_MINOR_BINARY='1'
# ENV DNNLROOT='/opt/intel/oneapi/dnnl/2021.2.0/cpu_dpcpp_gpu_dpcpp'
# ENV FI_PROVIDER_PATH='/opt/intel/oneapi/mpi/2021.2.0//libfabric/lib/prov:/usr/lib64/libfabric'
# ENV INFOPATH='/opt/intel/oneapi/debugger/10.1.1/documentation/info/'
# ENV INSPECTOR_2021_DIR='/opt/intel/oneapi/inspector/2021.2.0'
# ENV INTELFPGAOCLSDKROOT='/opt/intel/oneapi/compiler/2021.2.0/linux/lib/oclfpga'
# ENV INTEL_LICENSE_FILE='/opt/intel/licenses:/root/intel/licenses:/opt/intel/oneapi/clck/2021.2.0/licensing:/opt/intel/licenses:/root/intel/licenses:/Users/Shared/Library/Application Support/Intel/Licenses'
# ENV INTEL_PYTHONHOME='/opt/intel/oneapi/debugger/10.1.1/dep'
# ENV IPPCP_TARGET_ARCH='intel64'
# ENV IPPCRYPTOROOT='/opt/intel/oneapi/ippcp/2021.2.0'
# ENV IPPROOT='/opt/intel/oneapi/ipp/2021.2.0'
# ENV IPP_TARGET_ARCH='intel64'
# ENV I_MPI_ROOT='/opt/intel/oneapi/mpi/2021.2.0'
# ENV LD_LIBRARY_PATH='/opt/intel/oneapi/vpl/2021.2.2/lib:/opt/intel/oneapi/tbb/2021.2.0/env/../lib/intel64/gcc4.8:/opt/intel/oneapi/mpi/2021.2.0//libfabric/lib:/opt/intel/oneapi/mpi/2021.2.0//lib/release:/opt/intel/oneapi/mpi/2021.2.0//lib:/opt/intel/oneapi/mkl/latest/lib/intel64:/opt/intel/oneapi/itac/2021.2.0/slib:/opt/intel/oneapi/ippcp/2021.2.0/lib/intel64:/opt/intel/oneapi/ipp/2021.2.0/lib/intel64:/opt/intel/oneapi/dnnl/2021.2.0/cpu_dpcpp_gpu_dpcpp/lib:/opt/intel/oneapi/debugger/10.1.1/dep/lib:/opt/intel/oneapi/debugger/10.1.1/libipt/intel64/lib:/opt/intel/oneapi/debugger/10.1.1/gdb/intel64/lib:/opt/intel/oneapi/dal/2021.2.0/lib/intel64:/opt/intel/oneapi/compiler/2021.2.0/linux/lib:/opt/intel/oneapi/compiler/2021.2.0/linux/lib/x64:/opt/intel/oneapi/compiler/2021.2.0/linux/lib/emu:/opt/intel/oneapi/compiler/2021.2.0/linux/lib/oclfpga/host/linux64/lib:/opt/intel/oneapi/compiler/2021.2.0/linux/lib/oclfpga/linux64/lib:/opt/intel/oneapi/compiler/2021.2.0/linux/compiler/lib/intel64_lin:/opt/intel/oneapi/compiler/2021.2.0/linux/compiler/lib:/opt/intel/oneapi/ccl/2021.2.0/lib/cpu_gpu_dpcpp'
# ENV LIBRARY_PATH='/opt/intel/oneapi/vpl/2021.2.2/lib:/opt/intel/oneapi/tbb/2021.2.0/env/../lib/intel64/gcc4.8:/opt/intel/oneapi/mpi/2021.2.0//libfabric/lib:/opt/intel/oneapi/mpi/2021.2.0//lib/release:/opt/intel/oneapi/mpi/2021.2.0//lib:/opt/intel/oneapi/mkl/latest/lib/intel64:/opt/intel/oneapi/ippcp/2021.2.0/lib/intel64:/opt/intel/oneapi/ipp/2021.2.0/lib/intel64:/opt/intel/oneapi/dnnl/2021.2.0/cpu_dpcpp_gpu_dpcpp/lib:/opt/intel/oneapi/dal/2021.2.0/lib/intel64:/opt/intel/oneapi/compiler/2021.2.0/linux/compiler/lib/intel64_lin:/opt/intel/oneapi/compiler/2021.2.0/linux/lib:/opt/intel/oneapi/clck/2021.2.0/lib/intel64:/opt/intel/oneapi/ccl/2021.2.0/lib/cpu_gpu_dpcpp'
# ENV MANPATH='/opt/intel/oneapi/mpi/2021.2.0/man:/opt/intel/oneapi/itac/2021.2.0/man:/opt/intel/oneapi/debugger/10.1.1/documentation/man:/opt/intel/oneapi/clck/2021.2.0/man::/opt/intel/oneapi/compiler/2021.2.0/documentation/en/man/common:::'
# ENV MKLROOT='/opt/intel/oneapi/mkl/latest'
# ENV NLSPATH='/opt/intel/oneapi/mkl/latest/lib/intel64/locale/%l_%t/%N'
# ENV OCL_ICD_FILENAMES='libintelocl_emu.so:libalteracl.so:/opt/intel/oneapi/compiler/2021.2.0/linux/lib/x64/libintelocl.so'
# ENV ONEAPI_ROOT='/opt/intel/oneapi'
# ENV PATH='/opt/intel/oneapi/vtune/2021.2.0/bin64:/opt/intel/oneapi/vpl/2021.2.2/bin:/opt/intel/oneapi/mpi/2021.2.0//libfabric/bin:/opt/intel/oneapi/mpi/2021.2.0//bin:/opt/intel/oneapi/mkl/latest/bin/intel64:/opt/intel/oneapi/itac/2021.2.0/bin:/opt/intel/oneapi/itac/2021.2.0/bin:/opt/intel/oneapi/intelpython/latest/bin:/opt/intel/oneapi/intelpython/latest/condabin:/opt/intel/oneapi/inspector/2021.2.0/bin64:/opt/intel/oneapi/dev-utilities/2021.2.0/bin:/opt/intel/oneapi/debugger/10.1.1/gdb/intel64/bin:/opt/intel/oneapi/compiler/2021.2.0/linux/lib/oclfpga/llvm/aocl-bin:/opt/intel/oneapi/compiler/2021.2.0/linux/lib/oclfpga/bin:/opt/intel/oneapi/compiler/2021.2.0/linux/bin/intel64:/opt/intel/oneapi/compiler/2021.2.0/linux/bin:/opt/intel/oneapi/compiler/2021.2.0/linux/ioc/bin:/opt/intel/oneapi/clck/2021.2.0/bin/intel64:/opt/intel/oneapi/advisor/2021.2.0/bin64:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin'
# ENV PKG_CONFIG_PATH='/opt/intel/oneapi/vtune/2021.2.0/include/pkgconfig/lib64:/opt/intel/oneapi/vpl/2021.2.2/lib/pkgconfig:/opt/intel/oneapi/mkl/latest/tools/pkgconfig:/opt/intel/oneapi/inspector/2021.2.0/include/pkgconfig/lib64:/opt/intel/oneapi/advisor/2021.2.0/include/pkgconfig/lib64:'
# ENV PYTHONPATH='/opt/intel/oneapi/advisor/2021.2.0/pythonapi'
# ENV SETVARS_COMPLETED='1'
# ENV SETVARS_VARS_PATH='/opt/intel/oneapi/vtune/latest/env/vars.sh'
# ENV TBBROOT='/opt/intel/oneapi/tbb/2021.2.0/env/..'
# ENV VPL_BIN='/opt/intel/oneapi/vpl/2021.2.2/bin'
# ENV VPL_INCLUDE='/opt/intel/oneapi/vpl/2021.2.2/include'
# ENV VPL_LIB='/opt/intel/oneapi/vpl/2021.2.2/lib'
# ENV VPL_ROOT='/opt/intel/oneapi/vpl/2021.2.2'
# ENV VTUNE_PROFILER_2021_DIR='/opt/intel/oneapi/vtune/2021.2.0'
# ENV VT_ADD_LIBS='-ldwarf -lelf -lvtunwind -lm -lpthread'
# ENV VT_LIB_DIR='/opt/intel/oneapi/itac/2021.2.0/lib'
# ENV VT_MPI='impi4'
# ENV VT_ROOT='/opt/intel/oneapi/itac/2021.2.0'
# ENV VT_SLIB_DIR='/opt/intel/oneapi/itac/2021.2.0/slib'
# ENV _CE_CONDA=''
# ENV _CE_M=''

# Psi4 dependencies (apt)
# RUN apt-get install -y --no-install-recommends \
# libeigen3-dev libhdf5-dev

# Install Python 3.9
RUN conda install -y python=3.9 \
    && pip install --no-cache-dir pipx \
    && pipx uninstall pipx \
    && pipx reinstall-all

# Psi4 dependencies (conda)
RUN conda install -y -c psi4/label/dev hdf5 pybind11-headers ambit chemps2 dkh gau2grid gdma libint2 pcmsolver simint libxc \
    && conda install -y -c psi4/label/dev cmake eigen mpfr dataclasses msgpack-python networkx pytest pytest-xdist qcelemental qcengine