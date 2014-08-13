include ../../make.in

SRC = $(wildcard *.f)
SRC += $(wildcard *.F)

ifeq "${arch}" "mic"
OEDLIB = liboed_mic.a
OBJS0 := $(addsuffix .o, $(basename $(SRC)))
OBJS := $(addprefix mic/, $(OBJS0))
else
OEDLIB = liboed.a
OBJS := $(addsuffix .o, $(basename $(SRC)))
endif

FFLAGS = -O3 -fPIC
FFLAGS += -reentrancy threaded -recursive
FFLAGS += ${OPTFLAGS}

CFLAGS = -O3 -fPIC
CFLAGS += ${OPTFLAGS}

all: ${OEDLIB}

${OEDLIB}: ${OBJS}
	${AR} $@ ${OBJS}
	cp -r $@ ../../lib/

%.o: %.c
	${CC} ${CFLAGS} -c $< -o $@

%.o: %.f
	${FC} ${FFLAGS} -c $< -o $@

%.o: %.F
	${FC} ${FFLAGS} -c $< -o $@

mic/%.o: %.c
	${CC} ${CFLAGS} -c $< -o $@

mic/%.o: %.f
	${FC} ${FFLAGS} -c $< -o $@

mic/%.o: %.F
	${FC} ${FFLAGS} -c $< -o $@

clean:
	rm -f *.o mic/*.o *~ ${OEDLIB}
