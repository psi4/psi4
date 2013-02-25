V= 0.9.8-beta

include config.inc

all: efpmd

efpmd: libefp
	cd efpmd/src && $(MAKE)

tests: libefp
	cd tests && $(MAKE)

libefp:
	cd src && $(MAKE)

tags clean:
	cd src && $(MAKE) $@
	cd efpmd/src && $(MAKE) $@
	cd tests && $(MAKE) $@

check: tests
	@./tests/test

install: all
	install -d $(PREFIX)/bin
	install -d $(PREFIX)/include
	install -d $(PREFIX)/lib
	install -d $(FRAGLIB)
	install -m 0755 efpmd/src/efpmd $(PREFIX)/bin
	install -m 0644 src/efp.h $(PREFIX)/include
	install -m 0644 src/libefp.a $(PREFIX)/lib
	install -m 0644 fraglib/* $(FRAGLIB)

dist:
	git archive --format=tar.gz --prefix=libefp-$V/ -o libefp-$V.tar.gz HEAD

.PHONY: all efpmd tests libefp tags clean check install dist
