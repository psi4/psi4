TOP= ..

include $(TOP)/config.inc

CFLAGS= -I$(TOP)/src $(MYCFLAGS)
LDFLAGS= -L$(TOP)/src $(MYLDFLAGS)
LIBS= $(LAPACK) -lefp -lcheck -lm

PROG= test
ALL_C= disp_1a.c disp_1b.c disp_1c.c disp_1d.c disp_2a.c disp_2b.c disp_3a.c \
       disp_3b.c elec_1a.c elec_1b.c elec_1c.c elec_2a.c elec_2b.c elec_3a.c \
       elec_3b.c pol_1.c pol_2.c pol_3.c qm_1.c qm_2.c test.c total_1a.c \
       total_2a.c total_3a.c total_4a.c total_4b.c total_4c.c total_5a.c \
       total_5b.c total_5c.c total_6a.c total_6b.c total_6c.c xr_1a.c xr_1b.c \
       xr_2.c xr_3.c
ALL_O= disp_1a.o disp_1b.o disp_1c.o disp_1d.o disp_2a.o disp_2b.o disp_3a.o \
       disp_3b.o elec_1a.o elec_1b.o elec_1c.o elec_2a.o elec_2b.o elec_3a.o \
       elec_3b.o pol_1.o pol_2.o pol_3.o qm_1.o qm_2.o test.o total_1a.o \
       total_2a.o total_3a.o total_4a.o total_4b.o total_4c.o total_5a.o \
       total_5b.o total_5c.o total_6a.o total_6b.o total_6c.o xr_1a.o xr_1b.o \
       xr_2.o xr_3.o

all: $(PROG)

$(PROG): $(ALL_O)
	$(CC) -o $@ $(CFLAGS) $(LDFLAGS) $(ALL_O) $(LIBS)

clean:
	$(RM) $(PROG) $(ALL_O) $(TAGS)

$(TAGS): $(ALL_C)
	$(CTAGS) -f $(TAGS) $(ALL_C)

depend:
	@$(CC) $(CFLAGS) -MM $(ALL_C)

check: $(PROG)
	@cd .. && ./tests/$(PROG)

.PHONY: all clean depend check
