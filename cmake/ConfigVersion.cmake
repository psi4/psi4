# Configuration of Psi version
file(READ "${CMAKE_SOURCE_DIR}/VERSION" PSI_VERSION)
string(STRIP "${PSI_VERSION}" PSI_VERSION)
set(PSI_VERSION "\"${PSI_VERSION}\"")

# reset GIT_REVISION
set(GIT_REVISION)

# if GIT_HASH exists then this is exported code
# in this case we read git hash from this file and set DEVELOPMENT_CODE to false
if(EXISTS "${CMAKE_SOURCE_DIR}/cmake/GIT_HASH")
   file(READ "${CMAKE_SOURCE_DIR}/cmake/GIT_HASH" GIT_REVISION)
   string(STRIP "${GIT_REVISION}" GIT_REVISION)
   set(DEVELOPMENT_CODE FALSE) 
   set(PACKAGE_VERSION "${PSI_VERSION}")
else()
   set(DEVELOPMENT_CODE TRUE)
   set(PACKAGE_VERSION "${PSI_VERSION}git")
endif()

set(PACKAGE_NAME "Psi")
set(PACKAGE_STRING "${PACKAGE_NAME} ${PACKAGE_VERSION}")
set(PACKAGE_BUGREPORT "psicode@users.sourceforge.net")
