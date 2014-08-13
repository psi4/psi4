include ../make.in

SRC = $(wildcard *.c)
CFLAGS = -O3 -Wall -w2 -openmp -fPIC
CFLAGS += -Wunknown-pragmas -Wunused-variable
CFLAGS += ${OPTFLAGS}

ifeq "${arch}" "mic"
LIBCINT = libcint_mic.a
OBJS0 := $(addsuffix .o, $(basename $(SRC)))
OBJS := $(addprefix mic/, $(OBJS0))
else
LIBCINT = libcint.a
OBJS := $(addsuffix .o, $(basename $(SRC)))
endif

INC_ = -I./ -I../include

#INC += -I${MPI_INCDIR}

all: ${LIBCINT}

${LIBCINT}: ${OBJS}
	${AR} $@ $^
	cp -r $@ ../lib/
	sh genheader.sh

%.o : %.c Makefile
	$(CC) ${CFLAGS} ${INC} -c $< -o $@ 

mic/%.o : %.c Makefile
	$(CC) ${CFLAGS} ${INC} -c $< -o $@

clean:
	rm -f *.o mic/* *.s *.d *~ ${LIBCINT}
