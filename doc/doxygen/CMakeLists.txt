# Setup Doxygen documentation
if(DOXYGEN_FOUND)
   # Configuration of the Doxygen configuration file
   #set(PSI_VERSION "${PSI_VERSION}") # this is already defined
   set(top_scrdir ${CMAKE_SOURCE_DIR})
   set(HAVE_QHELP "NO")
   set(QHELP_GENERATOR "")
   set(HAVE_DOT "NO")
   set(DOT_PATH "")
   set(DOXY_OUT_DIR ${PROJECT_BINARY_DIR}/doc/doxygen)
   set(DOXY_OUT_FILE ${DOXY_OUT_DIR}/doxygen.cfg)
   configure_file(
       ${PROJECT_SOURCE_DIR}/doc/doxygen/doxygen.cfg.in
       ${DOXY_OUT_FILE})

 
   add_custom_target(doc ${DOXYGEN_EXECUTABLE}
                         ${DOXY_OUT_FILE}
                         WORKING_DIRECTORY ${DOXY_OUT_DIR}
                         COMMENT "Generating Doxygen code for Psi 4" VERBATIM
                    )
endif(DOXYGEN_FOUND)

