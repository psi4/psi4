include ../make.in

EXEC = scf

SRC = $(wildcard *.c)
OBJS := $(addsuffix .o, $(basename $(SRC)))

CFLAGS = -O3 -Wall -openmp
CFLAGS += -Wunknown-pragmas -Wunused-variable #-Wunused-function

#CFLAGS += -D_MYGA_

INC = -I./             \
      -I../include     \
      -I${MPI_INCDIR}  \
      -I${BLAS_INCDIR} \
      -I${GA_INCDIR}

LDFLAGS = -L${MPI_LIBDIR} -L${BLAS_LIBDIR} -L${GA_LIBDIR} -L${galoc}/lib

LIBS = ${MPI_LIBS} ${BLAS_LIBS} ${SCALAPACK_LIBS} ${GA_LIBS} \
       -lifcore -lintlc -limf    \
       ../lib/libpfock.a \
       ../lib/libcint.a  \
       /theoryfs2/ds/richard/lib/intel/mkl/lib/intel64/libmkl_scalapack_ilp64.a\
       ../lib/liberd.a   \
       ../lib/liboed.a 


all: ${EXEC} 

$(EXEC): $(OBJS) Makefile
	$(MPICC) ${CFLAGS} ${LDFLAGS} $(OBJS) -o $(EXEC) ${LIBS}

%.o : %.c Makefile
	$(MPICC) ${CFLAGS} ${INC} -c $< -o $@ 

clean:
	rm -f *.o *.s *.d *~ ${EXEC}
