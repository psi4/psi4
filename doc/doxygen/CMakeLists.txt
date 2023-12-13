find_package(Doxygen QUIET)
if(NOT TARGET Doxygen::doxygen)
   message(STATUS "No Doxygen, no docs.")
else()
    message(STATUS "Documentation targets available: doxyman (html)")

    # Configuration of the Doxygen configuration file
    set(HAVE_QHELP "NO")
    set(QHELP_GENERATOR "")
    set(HAVE_DOT "NO")
    set(DOT_PATH "")
    configure_file(doxygen.cfg.in doxygen.cfg @ONLY)

    add_custom_target(doxyman
                      VERBATIM
                      COMMAND Doxygen::doxygen ${CMAKE_CURRENT_BINARY_DIR}/doxygen.cfg
                      COMMENT "Build Doxygen html documentation ...")
endif()

