#SUBDIRS = lib apps
SUBDIRS = lib

libraries:	$(lib_LIBRARIES)
	for dir in $(SUBDIRS) ; do $(MAKE) -C $$dir $@ ; done
