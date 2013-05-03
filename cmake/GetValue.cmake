# A macro to read options in from the command line and environmental variables.

macro(get_value myVALUE myNAME myDEFAULT)
    # If either -DmyNAME or ENV{myNAME} were set, myVALUE is assigned that value
    # Start by grabbing the default
    set(${myVALUE} ${myDEFAULT})
    # Now try the environmental variable
    set(TMPVAL $ENV{${myNAME}})
    if(NOT ${TMPVAL} STREQUAL "")
        set(${myVALUE} ${TMPVAL})
    endif()
    # The -D overrides any environmental variable
    set(TMPVAL ${${myNAME}})
    if(NOT ${TMPVAL} STREQUAL "")
        set(${myVALUE} ${TMPVAL})
    endif()
    # Accumulate a list of all options
    set(ALL_USEROPTS "${ALL_USEROPTS};${myNAME}=${${myVALUE}}") 
endmacro(get_value)
