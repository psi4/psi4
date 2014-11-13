include(GenericMacros)
configure_script(cdash.py.in ${CMAKE_CURRENT_SOURCE_DIR}/cdash.py)
configure_script(cdash.run.in ${CMAKE_CURRENT_SOURCE_DIR}/cdash.run)

# vim:et:sw=4:ts=4:
