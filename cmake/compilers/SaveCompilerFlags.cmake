# Take care of updating the cache for fresh configurations

macro(save_compiler_flags lang)

if (${lang} STREQUAL C++)
    set (_lang CXX)
elseif(${lang} STREQUAL CXX)
    set (_lang CXX)
elseif(${lang} STREQUAL C)
    set (_lang C)
elseif(${lang} STREQUAL Fortran)
    set (_lang Fortran)
else()
    message(WARNING "Unknown language: ${lang}")
endif()

if (NOT DEFINED DEFAULT_${_lang}_FLAGS_SET)
    mark_as_advanced(DEFAULT_${_lang}_FLAGS_SET)
    set (DEFAULT_${_lang}_FLAGS_SET ON
        CACHE INTERNAL
        "Flag that the default ${_lang} compiler flags have been set.")

    set(CMAKE_${_lang}_FLAGS "${CMAKE_${_lang}_FLAGS}"
        CACHE STRING
        "Flags used by the compiler during all builds." FORCE)

    set(CMAKE_${_lang}_FLAGS_DEBUG "${CMAKE_${_lang}_FLAGS_DEBUG}"
        CACHE STRING
        "Flags used by the compiler during debug builds." FORCE)

    set(CMAKE_${_lang}_FLAGS_RELEASE "${CMAKE_${_lang}_FLAGS_RELEASE}"
        CACHE STRING
        "Flags used by the compiler during release builds." FORCE)
endif()
endmacro()
