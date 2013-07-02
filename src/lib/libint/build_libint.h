/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

/*! \file
    \ingroup INT
    \brief Enter brief description of file here 
*/
#define MAX_AM 22
#define DEFAULT_NEW_AM 8
#define DEFAULT_OPT_AM 8
#define DEFAULT_MAX_CLASS_SIZE 785

typedef struct {

  /* Twice the maximum AM for which manager routines need to be generated */
  int new_am;

  /* Twice the maximum AM for which VRR workers need to be generated */
  int opt_am;

  /* Twice the AM for which manager routines are already
     generated and stored in an existing library */
  int old_am;

  /* Max number of integrals to be processed in one VRR worker. If a class
   is larger than this then split the worker into several. */ 
  int max_class_size;

  /* The maximum AM for which VRR workers will be made inline functions */
  int max_am_to_inline_vrr_worker;

  /* The maximum AM of the VRR managers which will inline VRR workers */
  int max_am_manager_to_inline_vrr_worker;

  /* The maximum AM for which HRR workers will be made inline functions */
  int max_am_to_inline_hrr_worker;

  /* The maximum AM of the HRR managers which will inline HRR workers */
  int max_am_manager_to_inline_hrr_worker;

  /* The maximum AM for which VRR managers will be inlined into HRR managers */
  int max_am_to_inline_vrr_manager;

} LibintParams_t;

