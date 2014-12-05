TOP= ../..

include $(TOP)/config.inc

CFLAGS= $(MYCFLAGS)
FFLAGS= $(MYFFLAGS)

LIBOPT_A= libopt.a

ALL_C= opt.c
ALL_F= blas.f lbfgsb.f linpack.f timer.f
ALL_O= opt.o blas.o lbfgsb.o linpack.o timer.o

all: $(LIBOPT_A)

$(LIBOPT_A): $(ALL_O)
	$(AR) $@ $(ALL_O)
	$(RANLIB) $@

clean:
	$(RM) $(LIBOPT_A) $(ALL_O) $(TAGS)

$(TAGS): $(ALL_C) $(ALL_F)
	$(CTAGS) -f $(TAGS) $(ALL_C) $(ALL_F)

.PHONY: all clean
