#Here we define the expressions that need replaced in the doxygen file
set(STRINGS2FIND "\@PSI_VERSION\@"
                 "\@top_srcdir\@"
                 "\@HAVE_QHELP\@"
                 "\@QHELPGENERATOR\@"
                 "\@HAVE_DOT\@"
                 "\@DOT_PATH\@"
)

#What they should be replace with
set(STRINGS2REPLACE "${PSI_VERSION}"
                    "${CMAKE_SOURCE_DIR}"
                    "NO"
                    ""
                    "NO"
                    ""
)

#Read in the doxygen file, and repalce it's contents
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
FILE(WRITE doxygen.cfg ${FILE_CONTENT})

find_package(Doxygen QUIET)
if(DOXYGEN_FOUND)
   add_custom_target(doc ${DOXYGEN_EXECUTABLE}
                         ${CMAKE_CURRENT_SOURCE_DIR}/doxygen.cfg
                         WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/doc/doxygen
                         COMMENT "Generating Doxygen code for Psi 4" VERBATIM
)
   add_dependencies(doc ${CMAKE_CURRENT_SOURCE_DIR}/doxygen.cfg.in)
endif(DOXYGEN_FOUND)

