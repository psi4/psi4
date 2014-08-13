/**
 * @file ga++.h
 *
 * @author Manoj Kumar Krishnan, PNNL.
 * @author Jeff Daily, PNNL.
 */

/**
 * @mainpage
 *
 * @author Manoj Kumar Krishnan, PNNL.
 * @author Jeff Daily, PNNL.
 *
 * The GA Toolkit
 *
 * The Global Arrays (GA) toolkit provides an efficient and portable
 * “shared-memory” programming interface for distributed-memory computers.
 * Each process in a MIMD parallel program can asynchronously access logical
 * blocks of physically distributed dense multi-dimensional arrays, without
 * need for explicit cooperation by other processes. Unlike other
 * shared-memory environments, the GA model exposes to the programmer the
 * non-uniform memory access (NUMA) characteristics of the high performance
 * computers and acknowledges that access to a remote portion of the shared
 * data is slower than to the local portion. The locality information for the
 * shared data is available, and a direct access to the local portions of
 * shared data is provided.
 *
 * Global Arrays have been designed to complement rather than substitute for
 * the message-passing programming model. The programmer is free to use both
 * the shared-memory and message-passing paradigms in the same program, and to
 * take advantage of existing message-passing software libraries. Global
 * Arrays are compatible with the Message Passing Interface (MPI).
 *
 * The Global Arrays toolkit has been in the public domain since 1994. It has
 * been actively supported and employed in several large codes since then.
 */
#ifndef _GAPP_H
#define _GAPP_H

#include "ga.h"
#include "macdecls.h"

#define GANbhdl ga_nbhdl_t   

#include "init_term.h"
#include "services.h"
#include "PGroup.h"
#include "GlobalArray.h"
#include "GAServices.h"

#endif // _GAPP_H
