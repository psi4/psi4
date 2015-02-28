## This file should be placed in the root directory of your project.
## Then modify the CMakeLists.txt file in the root directory of your
## project to incorporate the testing dashboard.
## # The following are required to uses Dart and the Cdash dashboard
##   ENABLE_TESTING()
##   INCLUDE(CTest)
set(CTEST_PROJECT_NAME "Psi")
set(CTEST_NIGHTLY_START_TIME "00:00:00 CEST")

set(CTEST_DROP_METHOD "http")
set(CTEST_DROP_SITE "repo.theochem.kth.se")
set(CTEST_DROP_LOCATION "/cdash/submit.php?project=Psi")
set(CTEST_DROP_SITE_CDASH TRUE)

if(DEFINED ENV{CTEST_MAKE_NUM_PROCS})
    set(MAKECOMMAND "make -j$ENV{CTEST_MAKE_NUM_PROCS}" CACHE STRING "Custom make command")
endif()

# total allowed time for all tests in seconds
# see http://www.itk.org/Wiki/CTest:Using_CTEST_and_CDASH_without_CMAKE#General_CTest_settings
set(CTEST_TIMEOUT           "21600")
