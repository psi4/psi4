/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

///File is stolen from pulsar all thanks to bennyp

#ifndef PULSAR_GUARD_PULSAR__PRAGMA_H_
#define PULSAR_GUARD_PULSAR__PRAGMA_H_

#ifdef __cplusplus
extern "C" {
#endif



#if defined(__ICC) || defined(__INTEL_COMPILER)

    // pragmas for Intel
    #define PRAGMA_WARNING_POP                                 _Pragma("warning(pop)")
    #define PRAGMA_WARNING_PUSH                                _Pragma("warning(push)")
    #define PRAGMA_WARNING_IGNORE_UNUSED_PARAMETERS            _Pragma("warning(disable:869)")
    #define PRAGMA_WARNING_IGNORE_UNUSED_VARIABLES             _Pragma("warning(disable:177)")
    #define PRAGMA_WARNING_IGNORE_FP_EQUALITY                  _Pragma("warning(disable:1572)")
    #define PRAGMA_WARNING_IGNORE_FP_CONVERT                   _Pragma("warning(disable:264 173)")
    #define PRAGMA_WARNING_IGNORE_CONVERT                      _Pragma("warning(disable:2259)")
    #define PRAGMA_WARNING_IGNORE_SWITCH_MISSING_DEFAULT       _Pragma("warning(disable:2338)")
    #define PRAGMA_WARNING_IGNORE_POINTLESS_COMPARISON_UINT_0  _Pragma("warning(disable:186)")
    #define PRAGMA_WARNING_IGNORE_STATEMENT_UNREACHABLE        _Pragma("warning(disable:111)")
    #define PRAGMA_WARNING_IGNORE_SHADOW                       _Pragma("warning(disable:1599)")
    #define PRAGMA_WARNING_IGNORE_SHADOW_MEMBER                _Pragma("warning(disable:3280)")
    #define PRAGMA_WARNING_IGNORE_EXTRA_SEMICOLON              // does not have warning for intel
    #define PRAGMA_WARNING_IGNORE_REDECLARED_INLINE            _Pragma("warning(disable:522)")
    #define PRAGMA_WARNING_IGNORE_UNUSED_LOCAL_TYPEDEFS        //! \todo add me
    #define PRAGMA_WARNING_IGNORE_GCC_PRAGMA                   _Pragma("warning(disable:2282)")
    #define PRAGMA_WARNING_IGNORE_NONVIRTUAL_DTOR              _Pragma("warning(disable:444)")
    #define PRAGMA_WARNING_IGNORE_UNUSED_FUNCTION              //! \todo add me
    #define PRAGMA_WARNING_IGNORE_UNRECOGNIZED_PRAGMA          _Pragma("warning(disable:161)")
    #define PRAGMA_WARNING_IGNORE_DEPRECATED_DECLARATIONS      _Pragma("warning(disable:1478)")

#elif defined(__GNUC__) || defined(__GNUG__)

    // pragmas for GCC
    #define PRAGMA_WARNING_PUSH                                _Pragma("GCC diagnostic push")
    #define PRAGMA_WARNING_POP                                 _Pragma("GCC diagnostic pop")
    #define PRAGMA_WARNING_IGNORE_UNUSED_PARAMETERS            _Pragma("GCC diagnostic ignored \"-Wunused-parameter\"")
    #define PRAGMA_WARNING_IGNORE_UNUSED_VARIABLES             _Pragma("GCC diagnostic ignored \"-Wunused-variable\"")
    #define PRAGMA_WARNING_IGNORE_FP_EQUALITY                  _Pragma("GCC diagnostic ignored \"-Wfloat-equal\"")
    #define PRAGMA_WARNING_IGNORE_FP_CONVERT                   _Pragma("GCC diagnostic ignored \"-Wfloat-conversion\"")
    #define PRAGMA_WARNING_IGNORE_CONVERT                      _Pragma("GCC diagnostic ignored \"-Wconversion\"")
    #define PRAGMA_WARNING_IGNORE_SWITCH_MISSING_DEFAULT       _Pragma("GCC diagnostic ignored \"-Wswitch-default\"")
    #define PRAGMA_WARNING_IGNORE_POINTLESS_COMPARISON_UINT_0  _Pragma("GCC diagnostic ignored \"-Wtype-limits\"")
    #define PRAGMA_WARNING_IGNORE_STATEMENT_UNREACHABLE        //! \todo Is this a warning in GCC?
    #define PRAGMA_WARNING_IGNORE_SHADOW                       _Pragma("GCC diagnostic ignored \"-Wshadow\"")
    #define PRAGMA_WARNING_IGNORE_SHADOW_MEMBER                //! \todo doesn't seem to warn in GCC, or may be a part of -Wshadow 
    #define PRAGMA_WARNING_IGNORE_EXTRA_SEMICOLON              _Pragma("GCC diagnostic ignored \"-Wpedantic\"")
    #define PRAGMA_WARNING_IGNORE_REDECLARED_INLINE            // does not have warning for GCC
    #define PRAGMA_WARNING_IGNORE_UNUSED_LOCAL_TYPEDEFS        _Pragma("GCC diagnostic ignored \"-Wunused-local-typedefs\"")
    #define PRAGMA_WARNING_IGNORE_GCC_PRAGMA                   // uh... not a warning in gcc
    #define PRAGMA_WARNING_IGNORE_NONVIRTUAL_DTOR              // Doesn't seem to warn in GCC
    #define PRAGMA_WARNING_IGNORE_UNUSED_FUNCTION              _Pragma("GCC diagnostic ignored \"-Wunused-function\"")
    #define PRAGMA_WARNING_IGNORE_UNRECOGNIZED_PRAGMA          _Pragma("GCC diagnostic ignored \"-Wunknown-pragmas\"")
    #define PRAGMA_WARNING_IGNORE_DEPRECATED_DECLARATIONS      _Pragma("GCC diagnostic ignored \"-Wdeprecated-declarations\"")
#endif





#ifdef __cplusplus
}
#endif

#endif
