/*------------------------------------------------------------------------------*\
** ymlEstNrE.c
**
** Estimate number of elements needed in the main memory structure tbl[].
**
** The main memory structure tbl[] must be large enough for all the key:value
** pairs from the raw YAML file.  In order to allocate memory for tbl[], we
** need the number of elements multiplied by the size of each element.  The
** size of each element is given by the size of the yTB data structure.  But
** there is no way of knowing the number of elements until the YAML file is
** parsed.  At the same time, this software is what does the parsing.  It is
** a paradox.
**
** This subroutine makes a quick estimate of the number of key:value elements
** that will have to be stored in tbl[].  It does it by counting the number of
** value-delimiting characters in the raw YAML text file.  Those characters are
** colons (":"), commas (","), dashes ('-'), and square brackets ('[' and ']').
** Counting both square brackets will yield a slightly higher count.  For example,
** "[ a, b, c]" will yield a total of 4: two commas and two brackets.  A count
** of 3 would have sufficed because there are only three values.  Also, comments
** are not excluded when counting, so any of the delimiting characters in
** comments will add to the overall total.  And any of the delimiting characters
** within key or value strings will also add more than needed to the total.
** It would not be uncommon to see a YAML entry like this:
**
**     some-one-key:  "for example, a value three: 3"
**
** Such an entry would yield a count of 5 while 1 would have sufficed for the
** single value.  However, the algorithm is simple, the over-estimates will
** never be large, and allocating a few dozen bytes of additional memory will
** be better than not having enough.
**
**   LDx = level of detail
**   raw = pointer to raw YAML text
**   nrp = pointer to number of elements in tbl[]
\*------------------------------------------------------------------------------*/

#include "ymlPkgInc.h"			/* everything common to this package	*/

/*------------------------------------------------------------------------------*/

int ymlEstNrE(int LDx, char* raw, int* nrp)
{
    char* x = raw;			/* start of text that will be parsed	*/
    int   n = 0;			/* running count of desired characters	*/

    while ( *x )			/* loop until NULL terminator reached	*/
    {
       if ( *x == ':' |
            *x == '[' |
            *x == ']' |
            *x == ',' |
            *x == '-' ) { n++; }	/* found one char, increment count	*/
       x++;				/* move on to the next character	*/
    }

    *nrp = n;				/* return the count			*/

    return(0);				/* no errors can happen here		*/
}

/*------------------------------------------------------------------------------*/
