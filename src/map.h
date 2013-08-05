/* Declarations for System V style searching functions.
   Copyright (C) 1995-2013 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, see
   <http://www.gnu.org/licenses/>.  */

// Modified -- Josh Buckner -- 5 Aug 2013

#include <R.h>
#include <Rinternals.h>

/* Action which shall be performed in the call the hsearch.  */
typedef enum
  {
    FIND,
    ENTER
  }
ACTION;

typedef struct entry
  {
    char *key;
    void *data;
  }
ENTRY;

/* Opaque type for internal use.  */
struct _ENTRY;

/* Data type for reentrant functions.  */
struct hsearch_data
  {
    struct _ENTRY *table;
    unsigned int size;
    unsigned int filled;
  };

/* Reentrant versions which can handle multiple hashing tables at the
   same time.  */
extern int hsearch_r (ENTRY __item, ACTION __action, ENTRY **__retval,
		      struct hsearch_data *__htab);
extern int hcreate_r (size_t __nel, struct hsearch_data *__htab);
extern void hdestroy_r (struct hsearch_data *__htab);

// new convenient interfaces into hsearch_r and friends
// Josh Buckner 05 Aug 2013

/* the MAP type

   This struct was implemented to encapsulate the details of a hash map
   from a rowname label to a row index for an R dataframe, array or matrix.

   The GNU extension hash map was used so that the space for the hash map
   could be alloced from R's heap. So one needs to eventually Free the 
   entries as well as the hash table. This is why the struct is needed.
 */
typedef struct map {
  struct hsearch_data * hash_tab;
  ENTRY * entries;
  size_t n_entries;
} MAP;

/* function: create_map
   This function accepts a character vector SEXP and returns a hash map
   of a string in the character vector to it's position in the vector strs.
   The position is converted to a string for storage. The storage for the
   entries and the map are allocated using R's Calloc. The storage must be
   Freed using the delete_map function. The delete map function will Free
   each entry as well as the hash table.

   If SEXP strs is not a character vector, bad things will happen.
*/
MAP * create_map(SEXP strs);

/* function: delete_map
   Free's each entry of the MAP pointer strpos as well as the hash table.
   The storage must have been allocated with R's Calloc or bad things will
   happen.
 */
void delete_map(MAP * strpos);

/* function: get_pos
   consumes: a string to_find, a MAP pointer strpos allocated by
     create_map; if strpos has not been eaten by create_map first,
     bad things will happen
   returns: an integer resulting from applying the hash map defined
     by strpos to to_find. The integer will be returned from the GNU
     hsearch_r function as a string and will need to be converted to
     an integer before this function returns
 */
int get_pos(const char * to_find, MAP * strpos);

// end new stuff
