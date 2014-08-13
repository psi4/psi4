/**
 * @file init_term.h
 *
 * Ga Initialize and Terminate calls.
 */
#ifndef _INITTERM_H
#define _INITTERM_H

namespace GA {

/**
 * Initialize Global Arrays.
 * Allocate and initialize internal data structures in Global Arrays.
 * The limit is per process: it is the amount of memory that the given 
 * processor can contribute to collective allocation of global arrays. 
 * It does not include temporary storage that GA might be allocating (and 
 * releasing) during execution of a particular operation. 
 * limit < 0 means "allow unlimited memory usage" in which case this 
 * operation is equivalent to GA_initialize. This is a collective operation. 
 * @param argc,argv - command line argument lists.
 * @param limit - amount of memory in bytes per process [input]
 */
void Initialize(int argc, char *argv[], size_t limit = 0);


/**
 *Initialize Global Arrays.
 * Allocate and initialize internal data structures in Global Arrays.
 * The limit is per process: it is the amount of memory that the given 
 * processor can contribute to collective allocation of global arrays. 
 * It does not include temporary storage that GA might be allocating (and 
 * releasing) during execution of a particular operation. 
 * limit < 0 means "allow unlimited memory usage" in which case this 
 * operation is equivalent to GA_initialize. This is a collective operation. 
 * @param argc,argv - command line argument lists.
 * @param limit - amount of memory in bytes per process [input]
 * @param heapSize, stackSize - all of the dynamically allocated local memory 
 * @param type - data type.
 * in GA comes from its companion library,  the Memory Allocator (MA) library.
 * MA  allocates and manages local memory using stack and heap disciplines.
 * [refer section 3.2 of GA USer manual for more info]
 */
void Initialize(int argc, char *argv[], unsigned long heapSize, 
	   unsigned long stackSize, int type, size_t limit = 0);


/**
 * Delete all active arrays and destroy internal data structures. 
 * This is a collective operation. 
 */
void Terminate();

}

#endif /* _INITTERM_H */
