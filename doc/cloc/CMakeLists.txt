# Setup use of cloc perl script for counting lines of code
if(PERL_FOUND)
   add_custom_target(
       count_lines
       COMMAND "${PYTHON_EXECUTABLE}" "${PROJECT_BINARY_DIR}/bin/counter.py")
   add_custom_target(
       bar_charts 
       COMMAND "${PYTHON_EXECUTABLE}" "${PROJECT_BINARY_DIR}/bin/create_bar_charts.py")
   add_dependencies(bar_charts count_lines)
else()
   message(STATUS "Perl missing. You won't be able to create cloc stats.")
endif()
