set(Files2get "basis/*.gbs" "databases/*" "grids/*.grid" "plugin/*"
              "python/*.py" "python/qcdb/*.py" "python/p4const/*.py"
              "python/p4util/*.py" "quadratures/1_sqrtx/*" "quadratures/1_x/*")
set(WhereDayGo "basis" "databases" "grids" "plugin" "python" "python/qcdb"
               "python/p4const" "python/p4util" "quadratures/1_sqrtx"
               "quadratures/1_x")

set(counter 0)
foreach(FileGlob ${Files2get})
   FILE(GLOB things "${CMAKE_CURRENT_SOURCE_DIR}/${FileGlob}")
   LIST(GET WhereDayGo ${counter} destination)
   set(finaldest ${CMAKE_INSTALL_PREFIX}/share/psi/${destination})
   FILE(MAKE_DIRECTORY ${finaldest})
   INSTALL(FILES ${things} DESTINATION ${finaldest})
   MATH(EXPR counter "${counter}+1")
endforeach(FileGlob)

