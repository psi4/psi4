# C compiler
CC= clang

# Fortran compiler
FC= gfortran46

# install prefix
PREFIX= /home/ilya/libefp

# fragment library path
FRAGLIB= $(PREFIX)/libefp/fraglib

# additional link libraries
MYLIBS= -lalapack -lcblas -lf77blas -latlas -lgfortran -lmpi

# additional linker flags
MYLDFLAGS= -L/usr/local/lib -L/usr/local/lib/gcc46 -L/usr/local/mpi/openmpi/lib

# additional C flags
MYCFLAGS= -DWITH_MPI -std=c99 -g -Weverything -Wno-unreachable-code -Wno-c11-extensions -Wno-vla -Wno-sign-conversion -Wno-padded -Wno-format-nonliteral -Wno-disabled-macro-expansion -Wno-cast-align -I/usr/local/include -I/usr/local/mpi/openmpi/include

# additional Fortran flags
MYFFLAGS= -g

# ar command
AR= ar rcu

# ranlib command
RANLIB= ranlib

# rm command
RM= rm -f

# tags program
CTAGS= exctags

# tags file
TAGS= tags
