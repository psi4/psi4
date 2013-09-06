#
# Figure out which files libint will produce, and compile them
#
include(../CMakeLists.libint.txt)
set(LD_SRC "init_libderiv.cc" "d1hrr_order_0000.cc" "d1vrr_order_0000.cc" "d12hrr_order_0000.cc" "d12vrr_order_0000.cc")
build_libderiv_files(LD_SRC1 ${LIBDERIV_OPT_AM1} TRUE)
foreach (f ${LD_SRC1})
    set(LD_SRC ${LD_SRC} "d1${f}")
endforeach()
build_libderiv_files(LD_SRC2 ${LIBDERIV_OPT_AM2} FALSE)
foreach (f ${LD_SRC2})
    set(LD_SRC ${LD_SRC} "d12${f}")
endforeach()
set(ABC A B C D)
set(XYZ X Y Z)
math(EXPR MAXAM ${LIBDERIV_NEW_AM1}+1)
foreach(abcval ${ABC})
    foreach(xyzval ${XYZ})
        foreach(amval RANGE ${MAXAM})
            list(GET AM ${amval} amlabel)
            set(LD_SRC ${LD_SRC} "deriv_build_${abcval}${xyzval}_${amlabel}.cc")
        endforeach()
    endforeach()
endforeach()
add_library(deriv STATIC ${LD_SRC})

#
# Build the compiler
#
set(COMPILER_SRC build_libderiv.c emit_d1hrr_build_macro.c emit_deriv1_managers.c emit_deriv_build_macro.c
                 emit_d1hrr_build.c emit_deriv12_managers.c emit_deriv_build.c mem_man.c)
add_executable(libderiv_compiler ${COMPILER_SRC})
add_dependencies(libderiv_compiler int)

#
# Generate the libderiv files, after building the compiler, only if libderiv.compiled doesn't exist
#
add_custom_command(
    COMMAND libderiv_compiler
    COMMAND touch libderiv.compiled
    DEPENDS libderiv_compiler
    OUTPUT libderiv.compiled
    OUTPUT ${LD_SRC}
)


#
# This target is always executed
#
#add_custom_target(
#    libderiv ALL
#    DEPENDS libderiv.compiled
#)
