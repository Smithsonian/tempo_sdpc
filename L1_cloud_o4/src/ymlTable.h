/*------------------------------------------------------------------------------*\
** ymlTable.h -- not for users
**
** To call any functions in this package, use "ymlForApps.h" (not this file).
\*------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------*\
**
** The main memory structure tbl[] is defined here.  When used properly, the
** memory space for tbl[] becomes part of the user's memory space.  As a result,
** a user can have several of them, each for a different YAML file.
**
** A user should not access this memory structure directly.  Instead, use any of
** the ymlGet_X_() functions (see ymlForApps.h).
**
** The memory structure is organized as follows:
**
**    lvl     blk     usd        key                 val                 stk[]
** +-------+-------+-------+--------------+-------------------------+------------+
** |   0   |   0   |   0   |  Some Key 1  | Value 1 for Some Key 1  | key parts  |
** +-------+-------+-------+--------------+-------------------------+------------+
** |   0   |   0   |   1   |  Some Key 1  | Value 2 for Some Key 1  |    '\0'    |
** +-------+-------+-------+--------------+-------------------------+------------+
** |   0   |   0   |   4   |  Some Key 1  | Value 3 for Some Key 1  |    '\0'    |
** +-------+-------+-------+--------------+-------------------------+------------+
** |   1   |   1   |   0   |  Some Key 2  | Value 1 for Some Key 2  | key parts  |
** +-------+-------+-------+--------------+-------------------------+------------+
** |   2   |   2   |   1   |  Some Key 3  | Value 1 for Some Key 3  | key parts  |
** +-------+-------+-------+--------------+-------------------------+------------+
** |   2   |   2   |   0   |  Some Key 3  | Value 2 for Some Key 3  |    '\0'    |
** +-------+-------+-------+--------------+-------------------------+------------+
**
** The "key" column are null-terminated strings that represent multi-level tokens
** like "Input Files TC_L1A_EV_NASA" or "Runtime Parameters Start Time".
**
** The stk[] array holds the original unchanged parts of "key".  White space and
** capitalization are preserved.  For example, stk[0] might be "Runtime Parameters"
** and stk[1] might be "Start Time".
**
** Blank spaces and capitalization are preserved in the stk[] portion of the
** table, but not in the "key" column.  When a user calls one of the ymlGet_X_()
** functions, the keys are "normalized" so that blank spaces and capitalization
** in the key string are ignored.  So, for example, setting key to "InPuTFiles
** FoO bAr" is equivalent to "Input Files Foo Bar" or "INPUTFILESFOOBAR".
**
** The "val" column are null-terminated strings that are the values associated
** with each key.
**
** The "lvl" is the YAML block level number (showing level of indentation,
** starting and ending with 0).
**
** The "blk" is the YAML block number (starting with 0 and always increasing).
**
** The "usd" counts the number of times an item was retrieved (used) by user.
**
** tbl[] is not initialized until the user makes a call to ymlOpen() or one of
** the ymlGet_X_() functions.
\*------------------------------------------------------------------------------*/

#define YML_MAX_STK 10			/* max YAML block levels		*/

typedef struct				/* table for YAML data			*/
{
   int   lvl;				/* YAML level number			*/
   int   blk;				/* YAML block number			*/
   int   usd;				/* item use count			*/
   char  key[C_KEY_LEN];		/* parameter name			*/
   char  val[C_VAL_LEN];		/* parameter value			*/
   char* stk[YML_MAX_STK];		/* pointers to original key parts	*/
}  yTB;

/*------------------------------------------------------------------------------*/
