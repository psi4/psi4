#if HAVE_CONFIG_H
#   include "config.h"
#endif

/*
 * $Id: testc.c,v 1.5 1999-05-27 16:31:16 d3h325 Exp $
 */

/*
 * Test harness for portable dynamic memory allocator.
 */

#include "macdecls.h"
#include "string-util.h"
#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif

/**
 ** constants
 **/

#define NUM_COMMANDS    (int)(sizeof(commands) / sizeof(char *))

/**
 ** types
 **/

typedef enum
{
    C_MA_alloc_get = 0,
    C_MA_allocate_heap,
    C_MA_chop_stack,
    C_MA_free_heap,
    C_MA_get_index,
    C_MA_get_next_memhandle,
    C_MA_get_pointer,
    C_MA_init,
    C_MA_init_memhandle_iterator,
    C_MA_inquire_avail,
    C_MA_inquire_heap,
    C_MA_inquire_stack,
    C_MA_pop_stack,
    C_MA_print_stats,
    C_MA_push_get,
    C_MA_push_stack,
    C_MA_set_auto_verify,
    C_MA_set_error_print,
    C_MA_set_hard_fail,
    C_MA_sizeof,
    C_MA_sizeof_overhead,
    C_MA_summarize_allocated_blocks,
    C_MA_verify_allocator_stuff,
    C_bye,
    C_exit,
    C_help,
    C_quit,
    C_questionmark
} Command;

/**
 ** variables
 **/

static char *commands[] =
{
    "MA_alloc_get",
    "MA_allocate_heap",
    "MA_chop_stack",
    "MA_free_heap",
    "MA_get_index",
    "MA_get_next_memhandle",
    "MA_get_pointer",
    "MA_init",
    "MA_init_memhandle_iterator",
    "MA_inquire_avail",
    "MA_inquire_heap",
    "MA_inquire_stack",
    "MA_pop_stack",
    "MA_print_stats",
    "MA_push_get",
    "MA_push_stack",
    "MA_set_auto_verify",
    "MA_set_error_print",
    "MA_set_hard_fail",
    "MA_sizeof",
    "MA_sizeof_overhead",
    "MA_summarize_allocated_blocks",
    "MA_verify_allocator_stuff",
    "bye",
    "exit",
    "help",
    "quit",
    "?"
};

static char *help[] =
{
    "MA_alloc_get(datatype, nelem, name, memhandle, index)",
    "MA_allocate_heap(datatype, nelem, name, memhandle)",
    "MA_chop_stack(memhandle)",
    "MA_free_heap(memhandle)",
    "MA_get_index(memhandle, index)",
    "MA_get_next_memhandle(ithandle, memhandle)",
    "MA_get_pointer(memhandle, pointer)",
    "MA_init(datatype, nominal_stack, nominal_heap)",
    "MA_init_memhandle_iterator(ithandle)",
    "MA_inquire_avail(datatype)",
    "MA_inquire_heap(datatype)",
    "MA_inquire_stack(datatype)",
    "MA_pop_stack(memhandle)",
    "MA_print_stats()",
    "MA_push_get(datatype, nelem, name, memhandle, index)",
    "MA_push_stack(datatype, nelem, name, memhandle)",
    "MA_set_auto_verify(value)",
    "MA_set_error_print(value)",
    "MA_set_hard_fail(value)",
    "MA_sizeof(datatype1, nelem1, datatype2)",
    "MA_sizeof_overhead(datatype)",
    "MA_summarize_allocated_blocks()",
    "MA_verify_allocator_stuff()",
    "bye",
    "exit",
    "help",
    "quit",
    "?"
};

/**
 ** public routines
 **/

/* ------------------------------------------------------------------------- */
/*
 * Interactive interface to MA routines.
 */
/* ------------------------------------------------------------------------- */

int main(argc, argv)
    int        argc;
    char    *argv[];
{
    char        s[81];  /* string buffer */
    int         c;      /* index of command */
    int         ii1;    /* in int buffer */
    int         ii2;    /* in int buffer */
    int         ii3;    /* in int buffer */
    Integer     oi1;    /* out int buffer */
    MA_AccessIndex oi2; /* out int buffer */
    Pointer     op;     /* out Pointer buffer */
    int         value;  /* return value buffer */

    while (1)
    {
        printf("testc> ");

        if (scanf("%80s", s) != 1)
        {
            printf("*** Input read failed; exiting.\n");
            exit(1);
        }

        c = str_match(s, commands, NUM_COMMANDS);
        if (c == SM_MANY)
        {
            printf("*** '%s' is ambiguous\n", s);
            continue;
        }
        else if (c == SM_NONE)
        {
            printf("*** Unrecognized command '%s'\n", s);
            continue;
        }

        switch ((Command)c)
        {
            case C_MA_alloc_get:
                if (scanf("%d %d %80s", &ii1, &ii2, s) != 3)
                    printf("*** Input read failed for %s\n", commands[c]);
                else
                {
                    value = MA_alloc_get(ii1, ii2, s, &oi1, &oi2);
                    printf("%s=%d memhandle=%ld index=%ld\n",
                        commands[c], value, (long)oi1, (long)oi2);
                }
                break;
            case C_MA_allocate_heap:
                if (scanf("%d %d %80s", &ii1, &ii2, s) != 3)
                    printf("*** Input read failed for %s\n", commands[c]);
                else
                {
                    value = MA_allocate_heap(ii1, ii2, s, &oi1);
                    printf("%s=%d memhandle=%ld\n",
                            commands[c], value, (long)oi1);
                }
                break;
            case C_MA_chop_stack:
                if (scanf("%d", &ii1) != 1)
                    printf("*** Input read failed for %s\n", commands[c]);
                else
                {
                    value = MA_chop_stack(ii1);
                    printf("%s=%d\n", commands[c], value);
                }
                break;
            case C_MA_free_heap:
                if (scanf("%d", &ii1) != 1)
                    printf("*** Input read failed for %s\n", commands[c]);
                else
                {
                    value = MA_free_heap(ii1);
                    printf("%s=%d\n", commands[c], value);
                }
                break;
            case C_MA_get_index:
                if (scanf("%d", &ii1) != 1)
                    printf("*** Input read failed for %s\n", commands[c]);
                else
                {
                    value = MA_get_index(ii1, &oi2);
                    printf("%s=%d index=%ld\n", commands[c], value, (long)oi2);
                }
                break;
            case C_MA_get_next_memhandle:
                if (scanf("%d", &ii1) != 1)
                    printf("*** Input read failed for %s\n", commands[c]);
                else
                {
                    Integer tmp = ii1;
                    value = MA_get_next_memhandle(&tmp, &oi1);
                    printf("%s=%d ithandle=%ld memhandle=%ld\n",
                        commands[c], value, (long)tmp, (long)oi1);
                }
                break;
            case C_MA_get_pointer:
                if (scanf("%d", &ii1) != 1)
                    printf("*** Input read failed for %s\n", commands[c]);
                else
                {
                    value = MA_get_pointer(ii1, &op);
                    printf("%s=%d pointer=%lu\n", commands[c], value,
                        (unsigned long)op);
                }
                break;
            case C_MA_init:
                if (scanf("%d %d %d", &ii1, &ii2, &ii3) != 3)
                    printf("*** Input read failed for %s\n", commands[c]);
                else
                {
                    value = MA_init(ii1, ii2, ii3);
                    printf("%s=%d\n", commands[c], value);
                }
                break;
            case C_MA_init_memhandle_iterator:
                value = MA_init_memhandle_iterator(&oi1);
                printf("%s=%d ithandle=%ld\n", commands[c], value, (long)oi1);
                break;
            case C_MA_inquire_avail:
                if (scanf("%d", &ii1) != 1)
                    printf("*** Input read failed for %s\n", commands[c]);
                else
                {
                    oi1 = MA_inquire_avail(ii1);
                    printf("%s=%ld\n", commands[c], (long)oi1);
                }
                break;
            case C_MA_inquire_heap:
                if (scanf("%d", &ii1) != 1)
                    printf("*** Input read failed for %s\n", commands[c]);
                else
                {
                    oi1 = MA_inquire_heap(ii1);
                    printf("%s=%ld\n", commands[c], (long)oi1);
                }
                break;
            case C_MA_inquire_stack:
                if (scanf("%d", &ii1) != 1)
                    printf("*** Input read failed for %s\n", commands[c]);
                else
                {
                    oi1 = MA_inquire_stack(ii1);
                    printf("%s=%ld\n", commands[c], (long)oi1);
                }
                break;
            case C_MA_pop_stack:
                if (scanf("%d", &ii1) != 1)
                    printf("*** Input read failed for %s\n", commands[c]);
                else
                {
                    value = MA_pop_stack(ii1);
                    printf("%s=%d\n", commands[c], value);
                }
                break;
            case C_MA_print_stats:
                MA_print_stats(1);
                break;
            case C_MA_push_get:
                if (scanf("%d %d %80s", &ii1, &ii2, s) != 3)
                    printf("*** Input read failed for %s\n", commands[c]);
                else
                {
                    value = MA_push_get(ii1, ii2, s, &oi1, &oi2);
                    printf("%s=%d memhandle=%ld index=%ld\n",
                        commands[c], value, (long)oi1, (long)oi2);
                }
                break;
            case C_MA_push_stack: 
                if (scanf("%d %d %80s", &ii1, &ii2, s) != 3)
                    printf("*** Input read failed for %s\n", commands[c]);
                else
                {
                    value = MA_push_stack(ii1, ii2, s, &oi1);
                    printf("%s=%d memhandle=%ld\n",
                            commands[c], value, (long)oi1);
                }
                break;
            case C_MA_set_auto_verify:
                if (scanf("%d", &ii1) != 1)
                    printf("*** Input read failed for %s\n", commands[c]);
                else
                {
                    value = MA_set_auto_verify(ii1);
                    printf("%s=%d\n", commands[c], value);
                }
                break;
            case C_MA_set_error_print:
                if (scanf("%d", &ii1) != 1)
                    printf("*** Input read failed for %s\n", commands[c]);
                else
                {
                    value = MA_set_error_print(ii1);
                    printf("%s=%d\n", commands[c], value);
                }
                break;
            case C_MA_set_hard_fail:
                if (scanf("%d", &ii1) != 1)
                    printf("*** Input read failed for %s\n", commands[c]);
                else
                {
                    value = MA_set_hard_fail(ii1);
                    printf("%s=%d\n", commands[c], value);
                }
                break;
            case C_MA_sizeof:
                if (scanf("%d %d %d", &ii1, &ii2, &ii3) != 3)
                    printf("*** Input read failed for %s\n", commands[c]);
                else
                {
                    oi1 = MA_sizeof(ii1, ii2, ii3);
                    printf("%s=%ld\n", commands[c], (long)oi1);
                }
                break;
            case C_MA_sizeof_overhead:
                if (scanf("%d", &ii1) != 1)
                    printf("*** Input read failed for %s\n", commands[c]);
                else
                {
                    oi1 = MA_sizeof_overhead(ii1);
                    printf("%s=%ld\n", commands[c], (long)oi1);
                }
                break;
            case C_MA_summarize_allocated_blocks:
                MA_summarize_allocated_blocks();
                break;
            case C_MA_verify_allocator_stuff:
                value = MA_verify_allocator_stuff();
                printf("%s=%d\n", commands[c], value);
                break;
            case C_bye:
                /* fall through */
            case C_exit:
                /* fall through */
            case C_quit:
                exit(0);
            case C_help:
                /* fall through */
            case C_questionmark:
                for (ii1 = 0; ii1 < NUM_COMMANDS; ii1++)
                    printf(" %s\n", help[ii1]);
                break;
            default:
                printf("*** Unrecognized case '%d' in switch\n", c);
        }
    }

    return 0;
}
