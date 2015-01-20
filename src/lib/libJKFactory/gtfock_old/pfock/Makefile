include ../make.in

LIBPFOCK = libpfock.a

SRC = $(wildcard *.c)
OBJS := $(addsuffix .o, $(basename $(SRC)))

ifeq "${DYNAMIC}" "1"
CFLAGS += -D__DYNAMIC__
endif

ifeq "${DO_SCF}" "1"
CFLAGS += -D__SCF__
endif

INC = -I./                  \
      -I../include          \
      -I${BLAS_INCDIR}      \
      -I${MPI_INCDIR}       \
      -I${GA_INCDIR}        \
      -I${ERD_INCDIR}

all: ${LIBPFOCK}
	mkdir -p ../install
	cp -r ${LIBPFOCK} ../install/
	cp -r pfock.h ../install/        

${LIBPFOCK}: ${OBJS}
	${AR} $@ $^

%.o : %.c Makefile
	$(MPICC) ${CFLAGS} ${INC} -c $< -o $@ 

clean:
	rm -f *.o *.s *.d *~ *.a
