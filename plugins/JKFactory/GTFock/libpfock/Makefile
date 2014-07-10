include ../make.in

SCF=$(DOSCF)
DYNAMIC=1

LIBPFOCK = libpfock.a

SRC = $(wildcard *.c)
OBJS := $(addsuffix .o, $(basename $(SRC)))

CFLAGS = -O3 -Wall -openmp -w2 -fPIC
CFLAGS += -Wunknown-pragmas -Wunused-variable #-Wunused-function

ifeq "${DYNAMIC}" "1"
CFLAGS += -D_DYNAMIC_
endif

ifeq "${SCF}" "1"
CFLAGS += -D_SCF_
endif

INC = -I./                  \
      -I../include          \
      -I${MPI_INCDIR}       \
      -I${BLAS_INCDIR}      \
      -I${GA_INCDIR}

all: ${LIBPFOCK}

${LIBPFOCK}: ${OBJS}
	${AR} $@ $^
	cp -r $@ ../lib/
	sh genheader.sh ${SCF}

%.o : %.c Makefile
	$(CC) ${CFLAGS} ${INC} -c $< -o $@ 

clean:
	rm -f *.o *.s *.d *~ ${LIBPFOCK}
