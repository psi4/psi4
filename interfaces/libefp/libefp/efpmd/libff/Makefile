TOP= ../..

include $(TOP)/config.inc

CFLAGS= $(MYCFLAGS)

LIBFF_A= libff.a
ALL_C= ff.c
ALL_O= ff.o

all: $(LIBFF_A)

$(LIBFF_A): $(ALL_O)
	$(AR) $@ $(ALL_O)
	$(RANLIB) $@

clean:
	$(RM) $(LIBFF_A) $(ALL_O) $(TAGS)

$(TAGS): $(ALL_C)
	$(CTAGS) -f $(TAGS) $(ALL_C)

.PHONY: all clean
