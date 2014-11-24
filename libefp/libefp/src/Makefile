TOP= ..

include $(TOP)/config.inc

CFLAGS= $(MYCFLAGS)

LIBEFP_A= libefp.a
ALL_C= ai_disp.c balance.c bvec.c clapack.c disp.c efp.c elec.c elec_terms.c \
       int.c log.c parse.c pol.c pol_direct.c stream.c swf.c util.c xr.c
ALL_O= ai_disp.o balance.o bvec.o clapack.o disp.o efp.o elec.o elec_terms.o \
       int.o log.o parse.o pol.o pol_direct.o stream.o swf.o util.o xr.o

all: $(LIBEFP_A)

$(LIBEFP_A): $(ALL_O)
	$(AR) $@ $(ALL_O)
	$(RANLIB) $@

clean:
	$(RM) $(LIBEFP_A) $(ALL_O) $(TAGS)

$(TAGS): $(ALL_C)
	$(CTAGS) -f $(TAGS) $(ALL_C)

depend:
	@$(CC) $(CFLAGS) -MM $(ALL_C)

.PHONY: all clean depend
