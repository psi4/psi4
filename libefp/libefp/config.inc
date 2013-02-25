# C compiler
CC= gcc

# Fortran compiler
FC= gfortran

# lapack libraries
LAPACK= -lalapack -lcblas -lf77blas -latlas

# install prefix
PREFIX= /usr/local

# fragment library path
FRAGLIB= $(PREFIX)/share/libefp/fraglib

# additional C flags
MYCFLAGS= -std=c99 -g -O2 -fopenmp -Wall -Wextra

# additional Fortran flags
MYFFLAGS= -O2 -fopenmp

# additional linker flags
MYLDFLAGS=

# ar command
AR= ar rcu

# ranlib command
RANLIB= ranlib

# rm command
RM= rm -f

# tags program
CTAGS= ctags

# tags file
TAGS= tags
