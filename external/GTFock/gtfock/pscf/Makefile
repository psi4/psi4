include ../make.in

EXEC = scf
SRC = $(wildcard *.c)
OBJS := $(addsuffix .o, $(basename $(SRC)))

ifeq "${SCF_TIMMING}" "1"
CFLAGS += -D__SCF_TIMING__
endif

ifeq "${SCF_OUT}" "1"
CFLAGS += -D__SCF_OUT__
endif


INC = -I./                  \
      -I../install          \
      -I${ERD_INCDIR}       \
      -I${MPI_INCDIR}       \
      -I${SCALAPACK_INCDIR} \
      -I${BLAS_INCDIR}       

LDFLAGS = -L${MPI_LIBDIR}       \
          -L${GA_LIBDIR}        \
          -L${ARMCI_LIBDIR}     \
          -L${SCALAPACK_LIBDIR} \
          -L${BLAS_LIBDIR} \
          -L${ERD_LIBDIR}

PFOCK_LIBS = ../install/libpfock.a

LIBS = ${MPI_LIBS}       \
       ${BLAS_LIBS}      \
       ${SCALAPACK_LIBS} \
       ${PFOCK_LIBS}     \
       ${ERD_LIBS}       \
       ${GA_LIBS}        \
       ${ARMCI_LIBS}     \
       -lifcore -limf

all: ${EXEC} 

$(EXEC): $(OBJS) Makefile ${CINT_LIBS} ${PFOCK_LIBS}
	$(MPICC) ${CFLAGS} ${LDFLAGS} $(OBJS) -o $(EXEC) ${LIBS}

%.o : %.c Makefile
	$(MPICC) ${CFLAGS} ${INC} -c $< -o $@ 

clean:
	rm -f *.o *.s *.d *~ ${EXEC}
