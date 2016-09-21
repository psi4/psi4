set(Folders2Copy "basis" "databases" "grids" "plugin" "python" "quadratures" "efpfrag" "fsapt" "scripts")
INSTALL(
  DIRECTORY ${Folders2Copy} 
  DESTINATION ${CMAKE_INSTALL_PREFIX}/share/psi4
  USE_SOURCE_PERMISSIONS
  PATTERN "*.pyc" EXCLUDE
)

