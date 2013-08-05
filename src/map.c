/* Copyright (C) 1993-2013 Free Software Foundation, Inc.
   This file is part of the GNU C Library.
   Contributed by Ulrich Drepper <drepper@gnu.ai.mit.edu>, 1993.

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

// let R handle all memory allocation
#include<R.h>
#include<Rinternals.h>

#include <errno.h>
#include <string.h>

#include "map.h"
#include "cuseful.h"

// new convenient interfaces into hsearch_r and friends
// Josh Buckner 05 Aug 2013

/* function: create_map
   This function accepts a character vector SEXP and returns a hash map
   of a string in the character vector to it's position in the vector strs.
   The position is converted to a string for storage. The storage for the
   entries and the map are allocated using R's Calloc. The storage must be
   Freed using the delete_map function. The delete map function will Free
   each entry as well as the hash table.

   If SEXP strs is not a character vector, bad things will happen.
*/
MAP * create_map(SEXP strs) {

  // TODO: check to see that strs is a character vector

  int
    n_strs = length(strs),
    n_map = ceil((4.0 * (double) n_strs) / 3.0);
  // hash must be >= 20% empty for efficient look up

  MAP * strpos = Calloc(1, MAP);
  strpos->hash_tab = Calloc(n_map, struct hsearch_data);
  if( 0 == hcreate_r(n_map, strpos->hash_tab) )
    error("In create_strpos: Failed to create hash table. Out of memory?");

  strpos->entries = Calloc(n_strs, ENTRY);
  strpos->n_entries = n_strs;

  ENTRY * inserted;
  for(int i = 0; i < n_strs; i++) {
    strpos->entries[i].key = (char *) CHAR(STRING_ELT(strs, i));
    strpos->entries[i].data = Calloc(1 + digits(i), char);
    sprintf(strpos->entries[i].data, "%d", i);
    if( 0 == hsearch_r(strpos->entries[i], ENTER, &inserted, strpos->hash_tab) )
	error("In create_strpos: Can't insert key. Table full?");
  }
  return strpos;
}

/* function: delete_map
   Free's each entry of the MAP pointer strpos as well as the hash table.
   The storage must have been allocated with R's Calloc or bad things will
   happen.
 */
void delete_map(MAP * strpos) {

  // TODO: check for NULL pointers (strpos members too)

  hdestroy_r(strpos->hash_tab);
  for(int i = 0; i < strpos->n_entries; i++)
    Free(strpos->entries[i].data);

  Free(strpos->entries);
  Free(strpos->hash_tab);
  Free(strpos);
}

/* function: get_pos
   consumes: a string to_find, a MAP pointer strpos allocated by
     create_map; if strpos has not been eaten by create_map first,
     bad things will happen
   returns: an integer resulting from applying the hash map defined
     by strpos to to_find. The integer will be returned from the GNU
     hsearch_r function as a string and will need to be converted to
     an integer before this function returns
 */
int get_pos(const char * to_find, MAP * strpos) {

  // TODO: check for NULL pointers (strpos members too)

  ENTRY to_find_e, * found;

  // ENTRY's key is not const but R character vectors are always const
  // so we need the cast to avoid compiler warnings.
  to_find_e.key = (char *) to_find;
  if( 0 == hsearch_r(to_find_e, FIND, &found, strpos->hash_tab) )
    error("In get_pos: String not found.");

  // convert hashed string to long and return
  return strtol(found->data, NULL, 0);
}

// end new stuff

/* [Aho,Sethi,Ullman] Compilers: Principles, Techniques and Tools, 1986
   [Knuth]            The Art of Computer Programming, part 3 (6.4)  */


/* The reentrant version has no static variables to maintain the state.
   Instead the interface of all functions is extended to take an argument
   which describes the current status.  */
typedef struct _ENTRY
{
  unsigned int used;
  ENTRY entry;
}
_ENTRY;


/* For the used double hash method the table size has to be a prime. To
   correct the user given table size we need a prime test.  This trivial
   algorithm is adequate because
   a)  the code is (most probably) called a few times per program run and
   b)  the number is small because the table must fit in the core  */
static int
isprime (unsigned int number)
{
  /* no even number will be passed */
  unsigned int div = 3;

  while (div * div < number && number % div != 0)
    div += 2;

  return number % div != 0;
}


/* Before using the hash table we must allocate memory for it.
   Test for an existing table are done. We allocate one element
   more as the found prime number says. This is done for more effective
   indexing as explained in the comment for the hsearch function.
   The contents of the table is zeroed, especially the field used
   becomes zero.  */
int
hcreate_r (nel, htab)
     size_t nel;
     struct hsearch_data *htab;
{
  /* Test for correct arguments.  */
  if (htab == NULL)
      return 0;

  /* There is still another table active. Return with error. */
  if (htab->table != NULL)
    return 0;

  /* We need a size of at least 3.  Otherwise the hash functions we
     use will not work.  */
  if (nel < 3)
    nel = 3;
  /* Change nel to the first prime number not smaller as nel. */
  nel |= 1;      /* make odd */
  while (!isprime (nel))
    nel += 2;

  htab->size = nel;
  htab->filled = 0;

  /* allocate memory and zero out */
  htab->table = (_ENTRY *) Calloc (htab->size + 1, _ENTRY);
  if (htab->table == NULL)
    return 0;

  /* everything went alright */
  return 1;
}

/* After using the hash table it has to be destroyed. The used memory can
   be freed and the local static variable can be marked as not used.  */
void
hdestroy_r (htab)
     struct hsearch_data *htab;
{
  /* Test for correct arguments.  */
  if (htab == NULL)
      return;

  /* Free used memory.  */
  Free (htab->table);

  /* the sign for an existing table is an value != NULL in htable */
  htab->table = NULL;
}

/* This is the search function. It uses double hashing with open addressing.
   The argument item.key has to be a pointer to an zero terminated, most
   probably strings of chars. The function for generating a number of the
   strings is simple but fast. It can be replaced by a more complex function
   like ajw (see [Aho,Sethi,Ullman]) if the needs are shown.

   We use an trick to speed up the lookup. The table is created by hcreate
   with one more element available. This enables us to use the index zero
   special. This index will never be used because we store the first hash
   index in the field used where zero means not used. Every other value
   means used. The used field can be used as a first fast comparison for
   equality of the stored and the parameter value. This helps to prevent
   unnecessary expensive calls of strcmp.  */
int
hsearch_r (item, action, retval, htab)
     ENTRY item;
     ACTION action;
     ENTRY **retval;
     struct hsearch_data *htab;
{
  unsigned int hval;
  unsigned int count;
  unsigned int len = strlen (item.key);
  unsigned int idx;

  /* Compute an value for the given string. Perhaps use a better method. */
  hval = len;
  count = len;
  while (count-- > 0)
    {
      hval <<= 4;
      hval += item.key[count];
    }
  if (hval == 0)
    ++hval;

  /* First hash function: simply take the modul but prevent zero. */
  idx = hval % htab->size + 1;

  if (htab->table[idx].used)
    {
      /* Further action might be required according to the action value. */
      if (htab->table[idx].used == hval
	  && strcmp (item.key, htab->table[idx].entry.key) == 0)
	{
	  *retval = &htab->table[idx].entry;
	  return 1;
	}

      /* Second hash function, as suggested in [Knuth] */
      unsigned int hval2 = 1 + hval % (htab->size - 2);
      unsigned int first_idx = idx;

      do
	{
	  /* Because SIZE is prime this guarantees to step through all
             available indeces.  */
          if (idx <= hval2)
	    idx = htab->size + idx - hval2;
	  else
	    idx -= hval2;

	  /* If we visited all entries leave the loop unsuccessfully.  */
	  if (idx == first_idx)
	    break;

            /* If entry is found use it. */
          if (htab->table[idx].used == hval
	      && strcmp (item.key, htab->table[idx].entry.key) == 0)
	    {
	      *retval = &htab->table[idx].entry;
	      return 1;
	    }
	}
      while (htab->table[idx].used);
    }

  /* An empty bucket has been found. */
  if (action == ENTER)
    {
      /* If table is full and another entry should be entered return
	 with error.  */
      if (htab->filled == htab->size)
	{
	  *retval = NULL;
	  return 0;
	}

      htab->table[idx].used  = hval;
      htab->table[idx].entry = item;

      ++htab->filled;

      *retval = &htab->table[idx].entry;
      return 1;
    }

  *retval = NULL;
  return 0;
}
