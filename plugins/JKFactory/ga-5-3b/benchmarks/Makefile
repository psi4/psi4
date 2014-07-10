  GA_CONFIG := # e.g. /path/to/ga/install/bin/ga-config
ifneq ($(GA_CONFIG),)
CPPFLAGS := $(shell $(GA_CONFIG) --cppflags)
  FFLAGS := $(shell $(GA_CONFIG) --fflags)
  CFLAGS := $(shell $(GA_CONFIG) --cflags)
    LIBS := $(shell $(GA_CONFIG) --libs)
   FLIBS := $(shell $(GA_CONFIG) --flibs)
 LDFLAGS := $(shell $(GA_CONFIG) --ldflags)
      CC := $(shell $(GA_CONFIG) --cc)
     F77 := $(shell $(GA_CONFIG) --f77)
endif

dist:
	rm -rf ga-benchmarks
	mkdir ga-benchmarks
	cp -f ../global/testing/mp3.h      ./ga-benchmarks/mp3.h
	cp -f ../global/testing/mp3.fh     ./ga-benchmarks/mp3.fh
	cp -f ../global/testing/ga_shift.F ./ga-benchmarks/ga_shift.F
	cp -f ../global/testing/perf.F     ./ga-benchmarks/ga_ptp.F
	cp -f ../global/testing/perf2.c    ./ga-benchmarks/ga_perf.c
	cp -f ../armci/testing/perf.c      ./ga-benchmarks/armci_perf.c
	cp -f ./Makefile                   ./ga-benchmarks/Makefile
	cp -f ./config.h                   ./ga-benchmarks/config.h
	cp -f ./config.fh                  ./ga-benchmarks/config.fh
	cp -f ./util.c                     ./ga-benchmarks/util.c
	cp -f ./testutil.fh                ./ga-benchmarks/testutil.fh
	cp -f ./README                     ./ga-benchmarks/README
	sed -i '13,31d'                    ./ga-benchmarks/Makefile
	rm -f ga-benchmarks.tgz
	tar -czf ga-benchmarks.tgz ./ga-benchmarks

.SUFFIXES:

benchmarks: ga_shift.x ga_ptp.x ga_perf.x armci_perf.x

ga_shift.x: ga_shift.o util.o
	$(F77) -o $@ $^ $(LDFLAGS) $(LIBS)
ga_shift.o: ga_shift.F mp3.fh config.fh
	$(F77) -o $@ -c $< -DHAVE_CONFIG_H -I. $(CPPFLAGS) $(FFLAGS)

ga_ptp.x: ga_ptp.o util.o
	$(F77) -o $@ $^ $(LDFLAGS) $(LIBS)
ga_ptp.o: ga_ptp.F mp3.fh config.fh
	$(F77) -o $@ -c $< -DHAVE_CONFIG_H -I. $(CPPFLAGS) $(FFLAGS)

ga_perf.x: ga_perf.o util.o
	$(CC) -o $@ $^ $(LDFLAGS) $(LIBS) $(FLIBS)
ga_perf.o: ga_perf.c mp3.h config.h
	$(CC) -o $@ -c $< -DHAVE_CONFIG_H -I. $(CPPFLAGS) $(CFLAGS)

armci_perf.x: armci_perf.o util.o
	$(CC) -o $@ $^ $(LDFLAGS) $(LIBS) $(FLIBS)
armci_perf.o: armci_perf.c mp3.h config.h
	$(CC) -o $@ -c $< -DHAVE_CONFIG_H -I. $(CPPFLAGS) $(CFLAGS)

util.o: util.c
	$(CC) -o $@ -c $< -DHAVE_CONFIG_H -I. $(CPPFLAGS) $(CFLAGS)

config.h:
config.fh:
mp3.h:
mp3.fh:

clean:
	rm -f *.x *.o
