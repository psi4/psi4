
function(am_to_shell am repeat output)
    #                       1    1    2
    #             0    5    0    5    0
    set(_letters "spdfghiklmnoqrtuvwxyz")
    string(LENGTH "${_letters}" _length)
    if(am LESS 0 OR am GREATER_EQUAL _length)
        message(FATAL_ERROR "Invalid angular momentum: ${am}")
    endif()
    string(SUBSTRING "${_letters}" "${am}" 1 _letter)
    string(REPEAT "${_letter}" "${repeat}" _shell)
    set(${output} "${_shell}" PARENT_SCOPE)
endfunction()
