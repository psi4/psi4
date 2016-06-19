include(GNUInstallDirs)

find_package(HDF5 QUIET)

set(_Ambit_SEARCHES)

# Search AMBIT_ROOT first if it is set
if (AMBIT_ROOT)
    set(_Ambit_SEARCH_ROOT PATHS ${AMBIT_ROOT} NO_DEFAULT_PATH)
    list(APPEND _Ambit_SEARCHES _Ambit_SEARCH_ROOT)
endif ()

# Normal search
set(_Ambit_SEARCH_NORMAL
        PATHS "[HKEY_LOCAL_MACHINE\\SOFTWARE\\GnuWin32\\Ambit;InstallPath]"
        "$ENV{PROGRAMFILES}/ambit"
        "/usr"
        )
list(APPEND _Ambit_SEARCHES _Ambit_SEARCH_NORMAL)

set(Ambit_NAMES ambit)

# Try each search configuration.
foreach (search ${_Ambit_SEARCHES})
    find_path(Ambit_INCLUDE_DIR
            NAMES tensor.h
            ${${search}}
            PATH_SUFFIXES include include/ambit)
    find_library(Ambit_LIBRARY
            NAMES ${Ambit_NAMES}
            ${${search}}
            PATH_SUFFIXES lib lib64 ${CMAKE_INSTALL_LIBDIR})
endforeach ()

include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Ambit
        FOUND_VAR Ambit_FOUND
        REQUIRED_VARS Ambit_LIBRARY Ambit_INCLUDE_DIR
        HDF5_LIBRARIES HDF5_INCLUDE_DIRS)

if (Ambit_FOUND)
    set(Ambit_INCLUDE_DIRS ${PCMSolver_INCLUDE_DIR} ${HDF5_INCLUDE_DIRS})
    set(Ambit_LIBRARIES ${PCMSolver_LIBRARY} ${HDF5_LIBRARIES})

    if (NOT TARGET Ambit::Ambit)
        add_library(Ambit::Ambit UNKNOWN IMPORTED)
        set_target_properties(Ambit::Ambit PROPERTIES
                IMPORTED_LOCATION ${Ambit_LIBRARY}
                INTERFACE_LINK_LIBRARIES "${Ambit_LIBRARIES}"
                INTERFACE_INCLUDE_DIRECTORIES "${Ambit_INCLUDE_DIRS}")
    endif ()
    include_directories(SYSTEM "${Ambit_INCLUDE_DIRS}")
endif ()
