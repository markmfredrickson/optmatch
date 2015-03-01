#include <limits.h>
#include "cuseful.h"

/* function: digits
   computes the number characters required to store an integer base 10
   as a string. It will include a character for the negative sign. It
   will not include a character for the string termination symbol.

   This function should be much faster than the usual call to a
   logarithm function.
 */
int digits(int n) {
  if(n == INT_MIN) return 11;
  if(n < 0) return(1 + digits(-n));

  if(n >= 10000) {
    if(n >= 10000000) {
      if(n >= 100000000) {
	if(n >= 1000000000) return 10;
	return 9;
      }
      return 8;
    }
    if(n >= 100000) {
      if(n >= 1000000) return 7;
      return 6;
    }
    return 5;
  }
  if(n >= 100) {
    if(n >= 1000) return 4;
    return 3;
  }
  if(n >= 10) return 2;
  return 1;
}
