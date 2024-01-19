/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2024 The Psi4 Developers.
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

/// File is stolen from pulsar all thanks to bennyp

#ifndef PULSAR_GUARD_PULSAR__PRAGMA_H_
#define PULSAR_GUARD_PULSAR__PRAGMA_H_

// clang-format off
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
#define PRAGMA_WARNING_IGNORE_OVERLOADED_VIRTUAL

#elif defined(__clang__)  // Do clang before GNU because clang defines __GNUC__, too.

#define PRAGMA_WARNING_PUSH                                _Pragma("clang diagnostic push")
#define PRAGMA_WARNING_POP                                 _Pragma("clang diagnostic pop")
#define PRAGMA_WARNING_IGNORE_UNUSED_PARAMETERS            _Pragma("clang diagnostic ignored \"-Wunused-parameter\"")
#define PRAGMA_WARNING_IGNORE_UNUSED_VARIABLES             _Pragma("clang diagnostic ignored \"-Wunused-variable\"")
#define PRAGMA_WARNING_IGNORE_FP_EQUALITY                  _Pragma("clang diagnostic ignored \"-Wfloat-equal\"")
#define PRAGMA_WARNING_IGNORE_FP_CONVERT                   _Pragma("clang diagnostic ignored \"-Wfloat-conversion\"")
#define PRAGMA_WARNING_IGNORE_CONVERT                      _Pragma("clang diagnostic ignored \"-Wconversion\"")
#define PRAGMA_WARNING_IGNORE_SWITCH_MISSING_DEFAULT       _Pragma("clang diagnostic ignored \"-Wswitch-default\"")
#define PRAGMA_WARNING_IGNORE_POINTLESS_COMPARISON_UINT_0  _Pragma("clang diagnostic ignored \"-Wtype-limits\"")
#define PRAGMA_WARNING_IGNORE_STATEMENT_UNREACHABLE        //! \todo Is this a warning in clang?
#define PRAGMA_WARNING_IGNORE_SHADOW                       _Pragma("clang diagnostic ignored \"-Wshadow\"")
#define PRAGMA_WARNING_IGNORE_SHADOW_MEMBER                //! \todo doesn't seem to warn in clang, or may be a part of -Wshadow
#define PRAGMA_WARNING_IGNORE_EXTRA_SEMICOLON              _Pragma("clang diagnostic ignored \"-Wpedantic\"")
#define PRAGMA_WARNING_IGNORE_REDECLARED_INLINE            // does not have warning for clang
#define PRAGMA_WARNING_IGNORE_UNUSED_LOCAL_TYPEDEFS        _Pragma("clang diagnostic ignored \"-Wunused-local-typedefs\"")
#define PRAGMA_WARNING_IGNORE_GCC_PRAGMA                   // uh... not a warning in gcc
#define PRAGMA_WARNING_IGNORE_NONVIRTUAL_DTOR              // Doesn't seem to warn in clang
#define PRAGMA_WARNING_IGNORE_UNUSED_FUNCTION              _Pragma("clang diagnostic ignored \"-Wunused-function\"")
#define PRAGMA_WARNING_IGNORE_UNRECOGNIZED_PRAGMA          _Pragma("clang diagnostic ignored \"-Wunknown-pragmas\"")
#define PRAGMA_WARNING_IGNORE_DEPRECATED_DECLARATIONS      _Pragma("clang diagnostic ignored \"-Wdeprecated-declarations\"")
#define PRAGMA_WARNING_IGNORE_OVERLOADED_VIRTUAL           _Pragma("clang diagnostic ignored \"-Woverloaded-virtual\"")

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
#define PRAGMA_WARNING_IGNORE_OVERLOADED_VIRTUAL

#elif defined(_MSC_VER)

// pragmas for Microsoft Visual Compiler (MSVC)
#define PRAGMA_WARNING_PUSH
#define PRAGMA_WARNING_POP
#define PRAGMA_WARNING_IGNORE_UNUSED_PARAMETERS
#define PRAGMA_WARNING_IGNORE_UNUSED_VARIABLES
#define PRAGMA_WARNING_IGNORE_FP_EQUALITY
#define PRAGMA_WARNING_IGNORE_FP_CONVERT
#define PRAGMA_WARNING_IGNORE_CONVERT
#define PRAGMA_WARNING_IGNORE_SWITCH_MISSING_DEFAULT
#define PRAGMA_WARNING_IGNORE_POINTLESS_COMPARISON_UINT_0
#define PRAGMA_WARNING_IGNORE_STATEMENT_UNREACHABLE
#define PRAGMA_WARNING_IGNORE_SHADOW
#define PRAGMA_WARNING_IGNORE_SHADOW_MEMBER
#define PRAGMA_WARNING_IGNORE_EXTRA_SEMICOLON
#define PRAGMA_WARNING_IGNORE_REDECLARED_INLINE
#define PRAGMA_WARNING_IGNORE_UNUSED_LOCAL_TYPEDEFS
#define PRAGMA_WARNING_IGNORE_GCC_PRAGMA
#define PRAGMA_WARNING_IGNORE_NONVIRTUAL_DTOR
#define PRAGMA_WARNING_IGNORE_UNUSED_FUNCTION
#define PRAGMA_WARNING_IGNORE_UNRECOGNIZED_PRAGMA
#define PRAGMA_WARNING_IGNORE_DEPRECATED_DECLARATIONS
#define PRAGMA_WARNING_IGNORE_OVERLOADED_VIRTUAL

#endif
// clang-format on

// The following is adapted from https://gcc.gnu.org/wiki/Visibility the step-by-step guide at the very bottom
// Visibility macros
#if defined _WIN32 || defined __CYGWIN__
#define PSI_HELPER_SO_EXPORT __declspec(dllexport)
#define PSI_HELPER_SO_LOCAL
#else
#if __GNUC__ >= 4
#define PSI_HELPER_SO_EXPORT __attribute__((visibility("default")))
#define PSI_HELPER_SO_LOCAL __attribute__((visibility("hidden")))
#else
#define PSI_HELPER_SO_EXPORT
#define PSI_HELPER_SO_LOCAL
#endif
#endif

// Use generic helper definitions to define PSI_API and PSI_LOCAL
// PSI_API is used for the public API symbols.
// PSI_LOCAL is used for non-API symbols.
#define PSI_API PSI_HELPER_SO_EXPORT
#define PSI_LOCAL PSI_HELPER_SO_LOCAL

// Use in the header file as follows:
// PSI_DEPRECATED("extremely unsafe, use 'combust' instead!!!") void explode(void);
// will produce this kind output when compiling:
//    warning: 'explode' is deprecated: extremely unsafe, use 'combust' instead!!!
// The macro can similarly be used to deprecate variables.
// The implementation uses the standard attribute available in C++14
#define PSI_DEPRECATED(msg) [[deprecated(msg)]]

#endif
