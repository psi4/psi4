TOP= ../..

include $(TOP)/config.inc

CFLAGS= -I$(TOP)/src -DFRAGLIB_PATH="\"$(FRAGLIB)\"" $(MYCFLAGS)
FFLAGS= $(MYFFLAGS)
LDFLAGS= -L$(TOP)/src $(MYLDFLAGS)
LIBS= $(LAPACK) -lefp -lm

PROG= efpmd
ALL_C= cfg.c clapack.c common.c grad.c hess.c main.c \
       md.c opt.c optimizer.c parse.c rand.c sp.c
ALL_O= cfg.o clapack.o common.o grad.o hess.o main.o \
       md.o opt.o optimizer.o parse.o rand.o sp.o \
       blas.o lbfgsb.o linpack.o timer.o

all: $(PROG)

$(PROG): $(ALL_O)
	$(FC) -o $@ $(FFLAGS) $(LDFLAGS) $(ALL_O) $(LIBS)

clean:
	$(RM) $(PROG) $(ALL_O) $(TAGS)

$(TAGS): $(ALL_C)
	$(CTAGS) -f $(TAGS) $(ALL_C)

depend:
	@$(CC) $(CFLAGS) -MM $(ALL_C)

.PHONY: all clean depend

blas.o: lbfgsb/blas.f
	$(FC) -c $(FFLAGS) -o $@ $<

lbfgsb.o: lbfgsb/lbfgsb.f
	$(FC) -c $(FFLAGS) -o $@ $<

linpack.o: lbfgsb/linpack.f
	$(FC) -c $(FFLAGS) -o $@ $<

timer.o: lbfgsb/timer.f
	$(FC) -c $(FFLAGS) -o $@ $<
