/*------------------------------------------------------------------------------*\
** GetConfig.h -- Everything for using the CONFIG_READER
**
** The follwoing are functions that a user may call:
**
**   GetConfig_xxx
**   MetaConfig
**   ShareArgsWithFortran
**   ResetConfig
**
** GetConfig_xxx represents the following functions:
**
**   GetConfigString
**   GetConfigInteger
**   GetConfigReal
**   GetConfigDouble
**   GetConfigRealArr
**   GetConfigDoubleArr
**
** All the ymlGet_xxx functions are identical except for the parameter type that
** is returned (char, int, float, double).  Arrays (Arr) are one-dimensional
** lists of the same data types.
**
** Parameters:
**
**   in  flg = "E" or "e" specifies whether no value is deemed as an error
**   in  key = specifies what caller wants out of the control file
**   in  cnt = (only for Arr) count of elements in the array
**   out val = value returned to caller
**
** Return Codes:
**
**     1 - there is no value to retrieve (and flg was not "E" or "e")
**     0 - everything is OK, val has retrieved value
**    -1 - key string length exceeds CFG_KEY_LEN
**    -2 - key not found in control file
**    -3 - unable to open control file
**    -4 - control file opened, read wrong number of bytes
**    -5 - control file opened, but read failed
**    -6 - control file appears empty
**    -7 - control file has TAB characters (TABs are illegal in YAML)
**    -8 - YAML parser had problem initializing
**    -9 - YAML parser had problem while parsing
**   -10 - internal overflow of stk[KEY_STACK_LEN]
**   -11 - (not used)
**   -12 - internal overflow of tmp[CFG_KEY_LEN]
**   -13 - internal overflow of cfg[].val[CFG_VAL_LEN] (truncated; never returned)
**   -14 - there is no value to retrieve (and flg was "E" or "e")
**   -15 - failed to allocate memory for control file ingest
**   -16 - internal table not initialized; call ymlOpen() or ymlRdFile() first
**   -17 - file name is blank; there is no file to open
**   -18 - calling parameter out of valid range
**
** A negative return code indicates an error that is generally fatal.
**
** A positive return code indicates an unusual (non-nominal) condition that is
** generally not an error.
**
** The sizes given in this include for key and value (CFG_KEY_LEN and CFG_VAL_LEN)
** include the terminating null required by C.  The Fortran include file has
** CFG_KEY_LEN and CFG_VAL_LEN set to one less than the C include file so that key
** and value may be the same maximum size whether processed by C or Fortran.
**
** Each call to ymlGet_xxx writes the value, specified by key, into the caller's
** variable(s).  The caller is responsible for making the variable(s) large enough
** and of the correct data type.
**
** The first parameter in the call to GetConfig_xxx is a 1-character flg that
** may be either "E" (or "e") or something else (including a null string).  When
** it is "E" or "e" and there is no value being returned by GetConfig_xxx, a
** fatal error is raised.  When it is not "E" or "e" and no value is being
** returned, no error is raised.
**
** Each call to GetConfig_xxx writes the value, specified by key, into the
** caller's variable.  Multiple calls with the same key return subsequent values
** in a list for the key.  If the end of the list is reached, a null is written
** as the first string character instead of the value.
**
** The control file name (possibly a full path) needs to be the first command-line
** parameter when executing a program that uses GetConfig_xxx.  The control file
** is assumed to be in YAML format.
**
** Example 1:
**
**   #include "GetConfig.h"
**
**   char val[CFG_VAL_LEN];
**
**   GetConfigString("E", "Some Key Name", val);
**
** Example 2:
**
**   #include "GetConfig.h"
**
**   char key [CFG_KEY_LEN] = "Runtime Parameters Orbit Number";
**   char val [CFG_VAL_LEN];
**
**   int ier = GetConfigString("0", key, val);
**
**
**   if ( ier < 0 ) {
**      printf("  Error %d getting config <%s>\n", ier, key);
**   } else {
**      if ( ! *val ) { printf("    OK: no config value for <%s>\n", key); }
**      else          { printf("    OK: <%s> = <%s>\n", key, val); }
**   }
**
\*------------------------------------------------------------------------------*/

#include <hdf5.h>
#include "ymlForApps.h"

#define CFG_KEY_LEN YML_KEY_LEN	/* max key   string length (incl. end null)	*/
#define CFG_VAL_LEN YML_VAL_LEN	/* max value string length in control file	*/

extern char* ControlFile;	/* control file (null-terminated string)	*/

int GetConfigString    (char* flg, char* key, char*   val);
int GetConfigInteger   (char* flg, char* key, int*    val);
int GetConfigReal      (char* flg, char* key, float*  val);
int GetConfigDouble    (char* flg, char* key, double* val);
int GetConfigRealArr   (char* flg, char* key, float*  vls, const int cnt);
int GetConfigDoubleArr (char* flg, char* key, double* vls, const int cnt);
int MetaConfig         (hid_t* fid);
int GetConfigClose     ();
int ResetConfig        ();

/*------------------------------------------------------------------------------*/
