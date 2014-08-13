/** @file
 * Private header file containing symbolic constants and type declarations
 * for the table module.
 *
 * This file should only be included by internal C files.
 */
#ifndef _table_h
#define _table_h

#include "matypes.h"

/**
 ** constants
 **/

/* invalid handle */
#define TABLE_HANDLE_NONE (Integer)(-1)

/**
 ** types
 **/

/* type of data in each table entry */
typedef char * TableData;

/**
 ** function types
 **/

extern Integer ma_table_allocate(TableData data);
extern void ma_table_deallocate(Integer handle);
extern TableData ma_table_lookup(Integer handle);
extern Integer ma_table_lookup_assoc(TableData data);
extern Boolean ma_table_verify(Integer handle, char *caller);

#endif /* _table_h */
