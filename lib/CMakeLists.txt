set(Folders2Copy "basis" "databases" "grids" "plugin" "python" "quadratures" "efpfrag")
configure_file(${CMAKE_CURRENT_LIST_DIR}/python/inputparser.py.in ${CMAKE_CURRENT_LIST_DIR}/python/inputparser.py)
INSTALL(DIRECTORY ${Folders2Copy} DESTINATION
${CMAKE_INSTALL_PREFIX}/share/psi)
