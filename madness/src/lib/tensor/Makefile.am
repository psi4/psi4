include $(top_srcdir)/config/MakeGlobal.am

bin_PROGRAMS = test oldtest test_mtxmq test_Zmtxmq jimkernel test_scott
lib_LIBRARIES = libMADtensor.a

thisincludedir = $(includedir)/tensor
thisinclude_HEADERS = aligned.h     mxm.h     tensorexcept.h  tensoriter_spec.h  type_data.h \
                        basetensor.h  tensor.h        tensor_macros.h    vector_factory.h \
                        mtxmq.h     slice.h   tensoriter.h    tensor_spec.h vmath.h

test_SOURCES = test.cc
test_LDADD = libMADtensor.a $(LIBMISC) $(LIBWORLD) 

oldtest_SOURCES = oldtest.cc
oldtest_LDADD = libMADtensor.a $(LIBMISC) $(LIBWORLD) 

test_scott_SOURCES = test_scott.cc
test_scott_LDADD = libMADtensor.a $(LIBMISC) $(LIBWORLD) 

test_mtxmq_SOURCES = test_mtxmq.cc
test_mtxmq_LDADD = libMADtensor.a $(LIBWORLD)
test_mtxmq_CPPFLAGS = $(AM_CPPFLAGS) -DTIME_DGEMM

jimkernel_SOURCES = jimkernel.cc
jimkernel_LDADD = libMADtensor.a $(LIBMISC) $(LIBWORLD)


test_Zmtxmq_SOURCES = test_Zmtxmq.cc
test_Zmtxmq_LDADD = libMADtensor.a $(LIBWORLD)
test_Zmtxmq_CPPFLAGS = $(AM_CPPFLAGS) -DTIME_DGEMM

libMADtensor_a_SOURCES = tensor.cc tensoriter.cc basetensor.cc mtxmq.cc vmath.cc \
                         \
                        aligned.h     mxm.h     tensorexcept.h  tensoriter_spec.h  type_data.h \
                        basetensor.h  tensor.h        tensor_macros.h    vector_factory.h \
                        mtxmq.h     slice.h   tensoriter.h    tensor_spec.h vmath.h

if USE_X86_64_ASM
  libMADtensor_a_SOURCES += mtxmq_asm.S mtxm_gen.h genmtxm.py
mtxmq_asm.o:	mtxmq_asm.S mtxm_gen.h
	$(CC) $(CXXFLAGS) $(CPPFLAGS) $(AM_CPPFLAGS) -I. -DX86_64 -c -o $@ $<

mtxm_gen.h:	genmtxm.py 
	python $(srcdir)/genmtxm.py > $@
endif

if USE_X86_32_ASM
  libMADtensor_a_SOURCES += mtxmq_asm.S mtxm_gen.h genmtxm.py
mtxmq_asm.o:	mtxmq_asm.S mtxm_gen.h
	$(CC) $(CXXFLAGS) $(CPPFLAGS) $(AM_CPPFLAGS) -I. -DX86_32 -c -o $@ $<

mtxm_gen.h:	genmtxm.py 
	python $(srcdir)/genmtxm.py > $@
endif

tensoriter_spec.h tensor_spec.h:	tempspec.py
	python $(srcdir)/tempspec.py




#TESTPROGS = test test_mtxmq
#TEST1SRC = test.o
#TEST2SRC = test_mtxmq.o
#test_mtxmq:     test_mtxmq.o mtxmq.o mtxmq_asm.o
#        $(CXX) -o $@  $^

libraries:	$(lib_LIBRARIES)
	for dir in $(SUBDIRS) ; do $(MAKE) -C $$dir $@ ; done
