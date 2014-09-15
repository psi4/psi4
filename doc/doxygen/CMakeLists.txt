#Here we define the expressions that need replaced in the doxygen file
set(STRINGS2FIND "\@PSI_VERSION\@"
                 "\@top_srcdir\@"
                 "\@HAVE_QHELP\@"
                 "\@QHELPGENERATOR\@"
                 "\@HAVE_DOT\@"
                 "\@DOT_PATH\@"
)

#What they should be replaced with
set(STRINGS2REPLACE "${PSI_VERSION}"
                    "${CMAKE_SOURCE_DIR}"
                    "NO"
                    ""
                    "NO"
                    ""
)

#Read in the doxygen file, and replace its contents
FILE(READ doxygen.cfg.in FILE_CONTENT)
set(index 0)
foreach(mssg ${STRINGS2FIND})
   LIST(GET STRINGS2REPLACE ${index} REPLACE_TEXT)
   STRING(REGEX REPLACE "${mssg}" "${REPLACE_TEXT}"
          MODIFIED_FILE_CONTENT "${FILE_CONTENT}")
   set(FILE_CONTENT ${MODIFIED_FILE_CONTENT})
   MATH(EXPR index "${index}+1")
endforeach(mssg)

#Write the new contents into a new doxygen file
set(DOXY_OUT_DIR ${CMAKE_BINARY_DIR}/doc/doxygen)
FILE(MAKE_DIRECTORY ${DOXY_OUT_DIR})
set(DOXY_OUT_FILE ${DOXY_OUT_DIR}/doxygen.cfg)
FILE(WRITE ${DOXY_OUT_FILE} ${FILE_CONTENT})

find_package(Doxygen QUIET)
if(DOXYGEN_FOUND)
   add_custom_target(doc ${DOXYGEN_EXECUTABLE}
                         ${DOXY_OUT_FILE}
                         WORKING_DIRECTORY ${DOXY_OUT_DIR}
                         COMMENT "Generating Doxygen code for Psi 4" VERBATIM
)
add_dependencies(doc ${CMAKE_CURRENT_SOURCE_DIR}/doxygen.cfg.in)
endif(DOXYGEN_FOUND)

