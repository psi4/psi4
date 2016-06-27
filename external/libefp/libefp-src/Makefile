V= 1.2.1

include config.inc

all: efpmd

efpmd: libefp
	cd efpmd/libff && $(MAKE)
	cd efpmd/libopt && $(MAKE)
	cd efpmd/src && $(MAKE)

libefp:
	cd src && $(MAKE)

tags:
	cd src && $(MAKE) $@
	cd efpmd/src && $(MAKE) $@

clean:
	cd src && $(MAKE) $@
	cd tests && $(MAKE) $@
	cd efpmd/libff && $(MAKE) $@
	cd efpmd/libopt && $(MAKE) $@
	cd efpmd/src && $(MAKE) $@

check check-omp check-mpi: efpmd
	cd tests && $(MAKE) $@

install: all
	install -d $(PREFIX)/bin
	install -d $(PREFIX)/include
	install -d $(PREFIX)/lib
	install -d $(FRAGLIB)
	install -m 0755 efpmd/src/efpmd $(PREFIX)/bin
	install -m 0755 efpmd/tools/cubegen.pl $(PREFIX)/bin
	install -m 0755 efpmd/tools/trajectory.pl $(PREFIX)/bin
	install -m 0644 src/efp.h $(PREFIX)/include
	install -m 0644 src/libefp.a $(PREFIX)/lib
	install -m 0644 fraglib/* $(FRAGLIB)

dist:
	git archive --format=tar.gz --prefix=libefp-$V/ -o libefp-$V.tar.gz HEAD

.PHONY: all efpmd libefp tags clean check check-omp check-mpi install dist
