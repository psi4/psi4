all: libs

libs:
	make -C external
	make -C libcint
	make -C libpfock

scf:
	make -C pscf

test:
	make -C testprog

clean:
	make -C libcint clean
	make -C libpfock clean
	make -C testprog clean
	make -C pscf clean

cleanall: clean
	make -C external clean
