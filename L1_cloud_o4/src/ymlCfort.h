/*------------------------------------------------------------------------------*\
** ymlCfort.h -- not for users
**
** Things relevant to the C/Fortran interface.
\*------------------------------------------------------------------------------*/

#define   KEY_LEN   260		/* max key   string length (language-agnostic)	*/
#define   VAL_LEN   260		/* max value string length (language-agnostic)	*/
#define   STR_LEN   260		/* generic   string length (language-agnostic)	*/

#define C_KEY_LEN (KEY_LEN + 1)	/* max key   string length (C: end null added)	*/
#define C_VAL_LEN (VAL_LEN + 1)	/* max value string length (C: end null added)	*/
#define C_STR_LEN (STR_LEN + 1)	/* generic   string length (C: end null added)	*/

#define F_KEY_LEN (KEY_LEN)	/* max key   string length (Fortran)		*/
#define F_VAL_LEN (VAL_LEN)	/* max value string length (Fortran)		*/
#define F_STR_LEN (STR_LEN)	/* generic   string length (Fortran)		*/

/*------------------------------------------------------------------------------*/
