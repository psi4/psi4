
# SUBDIRS should catch all directories one level down except .git and .
ROOT_DIRECTORY=.
SUBDIRS := ${shell find ${ROOT_DIRECTORY} -maxdepth 1 -type d -print | grep -v git | grep -v '^.$$'}

subdirs:
	for dir in $(SUBDIRS); do \
	$(MAKE) -C $$dir; \
	done

clean:
	for dir in $(SUBDIRS); do \
    $(MAKE) clean -C $$dir; \
    done

