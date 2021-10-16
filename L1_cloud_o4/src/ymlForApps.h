/*------------------------------------------------------------------------------*\
** ymlForApps.h -- everything for using YAML in APPs
\*------------------------------------------------------------------------------*/

#include "ymlCfort.h"
#include "ymlTable.h"

/*------------------------------------------------------------------------------*\
** The sizes given in this include for key and value (YML_KEY_LEN and YML_VAL_LEN)
** include the terminating null required by C.  The Fortran include file has
** YML_KEY_LEN and YML_VAL_LEN set to one less than the C include file so that key
** and value may be the same maximum size whether processed by C or Fortran.
\*------------------------------------------------------------------------------*/

#define YML_KEY_LEN C_KEY_LEN		/* max key   str length (with end null)	*/
#define YML_VAL_LEN C_VAL_LEN		/* max value str length (with end null)	*/

/*------------------------------------------------------------------------------*\
** Each call to ymlGet...() writes the value, specified by key, into the caller's
** variable or array.  The caller is responsible for making the variable (or
** array) large enough and of the correct data type.
\*------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------*\
** ymlOpen()
**
** Reads a file and YAML-parses it into user-supplied table tbl[].
**
**   fnm = file name to read and parse
**   nrp = pointer to number of elements in tbl[]
**   ptr = pointer to tbl[], table for YAML data
\*------------------------------------------------------------------------------*/

int ymlOpen(char* fnm, int* nrp, yTB** ptr);

/*------------------------------------------------------------------------------*\
** ymlGetOneStr()
**
** Returns one string.
**
**   nry = number of elements in tbl[]
**   tbl = table for YAML data
**   key = desired parameter name
**   cnt = which of many values -- 0="next", 1=first, 2=second, etc.
**   val = string value (returned)
\*------------------------------------------------------------------------------*/

int ymlGetOneStr(int nry, yTB tbl[], char* key, int cnt, char* val);

/*------------------------------------------------------------------------------*\
** ymlGetArrStr()
**
** Returns an array of strings.
**
**   nry = number of elements in tbl[]
**   tbl = table for YAML data
**   key = desired parameter name
**   cnt = which of many values -- 0="next", 1=first, 2=second, etc.
**   inb = minimum block number from which to get the array
**   blk = actual block number from which array came
**   arm = maximum number of elements in arr[]
**   arr = array of strings (returned)
\*------------------------------------------------------------------------------*/

int ymlGetArrStr(int nry, yTB tbl[], char* key, int cnt, int inb, int* blk,
                                            int arm, char arr[arm][YML_VAL_LEN]);

/*------------------------------------------------------------------------------*\
** ymlGetKeys()
**
** Returns an array of unique strings that are somewhat like perl hash keys,
** although each includes all previous levels.  For example, suppose that the
** original YAML file included the following:
**
**  fee:
**    fi: 1
**    fo: 2
**  fum:  3
**
** Setting lvl=1 would return the following:
**
**   arr[0] = "FEEFI"           (i.e., fee -> fi)
**   arr[1] = "FEEFO"           (i.e., fee -> fo)
**   arr[2] = '\0'
** --
**   nry = number of elements in tbl[]
**   tbl = table with YAML data
**   lvl = desired level number (0, 1, 2, ...)
**   arm = maximum number of elements in arr[]
**   arr = array of strings (returned)
\*------------------------------------------------------------------------------*/

int ymlGetKeys(int nry, yTB tbl[], int lvl, int arm, char arr[arm][YML_VAL_LEN]);

/*------------------------------------------------------------------------------*\
** ymlClrUseFlg()
**
** Where tbl[].key matches the specified key, re-initializes the tbl[].usd flag
** so the same values may be read again.
**
**   nry = number of elements in tbl[]
**   tbl = table of YAML data
**   key = element to reset
\*------------------------------------------------------------------------------*/

int ymlClrUsdFlg(int nry, yTB tbl[], char* key);

/*------------------------------------------------------------------------------*\
**
** Return Codes:
**
**   0 - everything is OK, something meaningful is being returned
**
** A positive return code indicates an unusual (non-nominal) condition that is
** generally not an error:
**
**   1 - the key was found in the file, but there is no value to retrieve
**   2 - buf[YML_KEY_LEN] overflow; was truncated
**   3 - yTB[] overflow; error message says how much to increase
**   4 - string truncated; error message says how much to increase
**   5 - stk[] overflow; error message says how much to increase
**   6 - cnt > nry; searching past the table's end will return nothing
**   7 - arr[] overflow; error message says how much to increase
**   8 - yTB[].val[] overflow, error message says how much, was truncated
**   9 - Maximum Level Of Detail is not an integer
**
** A negative return code indicates an error that is generally fatal:
**
**  -1 - file name is blank; there is no file to open
**  -2 - failed to allocate memory for YAML file ingest
**  -3 - unable to open YAML file
**  -4 - YAML file opened, read wrong number of bytes
**  -5 - YAML file opened, but read failed
**  -6 - YAML file appears empty
**  -7 - YAML file has TAB characters (TABs are illegal in YAML)
**  -8 - YAML parser had problem initializing
**  -9 - YAML parser had problem while parsing
** -10 - table not initialized; call ymlOpen() first
** -11 - atexit() failed to register a cleanup routine
** -12 - key not found in YAML file
** -13 - count is negative; only positive values make sense
** -14 - max > YML_MAX_STK and that should never happen (coding error)
** -15 - parsing mode parameter must be either T or E
** -16 - YAML TOKEN not recognized; this should never happen (coding error)
** -17 - YAML EVENT not recognized; this should never happen (coding error)
** -18 - parsing went bad; YAML level got out of sync somehow
** -19 - stat() (getting size of YAML file) failed
**
\*------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------*\
** Example 1:
**
**   #include "ymlForApps.h"
**
**   main()
**   {
**      char val[YML_VAL_LEN];
**      int  nry;
**      yTB* tbl;
**
**      ymlOpen("control.txt", &nry, &tbl);
**      ymlGetOneStr(nry, tbl, "Input Files ESDT", 0, val);
**
**      printf("%s\n", val);
**
**      exit(0);
**   }
**
\*------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------*\
** Example 2:
**
**   #include "ymlForApps.h"
**
**   main()
**   {
**      char fn1[] = "file_1.yml";
**      char fn2[] = "file_2.yml";
**
**      char k1[YML_KEY_LEN] = "key for first file";
**      char v1[YML_VAL_LEN];
**
**      char k2[YML_KEY_LEN] = "some key for file 2";
**      char v2[YML_VAL_LEN];
**
**      int  nr1;
**      int  nr2;
**
**      yTB* tb1;
**      yTB* tb2;
**
**      int  rtc;
**
**      rtc = ymlOpen(fn1, &nr1, &tb1);
**
**      if ( rtc != 0 ) { printf("  Error %d opening <%s>\n", rtc, fn1); }
**      else            { printf("     OK opened and parsed %s\n", fn1); }
**
**      rtc = ymlOpen(fn2, &nr2, &tb2);
**
**      if ( rtc != 0 ) { printf("  Error %d opening <%s>\n", rtc, fn2); }
**      else            { printf("     OK opened and parsed %s\n", fn2); }
**
**      rtc = ymlGetOneStr(nr1, tb1, k1, 0, v1);
**
**      if ( rtc != 0 ) { printf("  Error %d getting value for <%s>\n", rtc, k1); }
**      else            { printf("     OK got value for <%s> -- <%s>\n", k1, v1); }
**
**      rtc = ymlGetOneStr(nr2, tb2, k2, 0, v2);
**
**      if ( rtc != 0 ) { printf("  Error %d getting value for <%s>\n", rtc, k2); }
**      else            { printf("     OK got value for <%s> -- <%s>\n", k2, v2); }
**
**      exit(rtc);
**   }
**
\*------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------*\
** Example 3:
**
**   #include "ymlForApps.h"
**
**   main()
**   {
**      char fnm[] = "file_3.yml";
**      int  nry;
**      yTB* tbl;
**
**      int rtc = ymlOpen(fnm, &nry, &tbl);
**
**      if ( rtc != 0 ) { printf("  Error %d opening <%s>\n", rtc, fnm); }
**      else            { printf("     OK opened and parsed %s\n", fnm); }
**
**      char key[YML_KEY_LEN] = "three values";
**      int  cnt = 0;
**      int  inb = 0;
**      int  blk;
**      int  arm = 3;
**      char arr[arm][YML_VAL_LEN];
**
**      rtc = ymlGetArrStr(nry, tbl, key, cnt, inb, &blk, arm, arr);
**
**      if ( rtc != 0 ) { printf("  Error %d getting value for <%s>\n", rtc, key); }
**      else
**      {
**         int idx; for( idx = 0; idx < arm; idx++ )
**         {
**             printf("     OK got value for <%s> -- <%s>\n", key, arr[idx]);
**         }
**      }
**
**      exit(rtc);
**   }
**
\*------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------*\
** Example 4:
**
**   #include "ymlForApps.h"
**
**   main()
**   {
**      char fnm[] = "file_4.yml";
**      int  nry;
**      yTB* tbl;
**
**      int rtc = ymlOpen(fnm, &nry, &tbl);
**
**      if ( rtc != 0 ) { printf("  Error %d opening <%s>\n", rtc, fnm); }
**      else            { printf("     OK opened and parsed %s\n", fnm); }
**
**      int  lvl = 1;
**      int  arm = 3;
**      char arr[arm][YML_VAL_LEN];
**
**      rtc = ymlGetKeys(nry, tbl, lvl, arm, arr);
**
**      if ( rtc != 0 ) { printf("  Error %d getting keys for <%s>\n", rtc, fnm); }
**      else
**      {
**         int idx; for( idx = 0; idx < arm; idx++ )
**         {
**             if ( arr[idx][0] == '\0' ) { break; }
**             printf("     OK got key <%s>\n", arr[idx]);
**         }
**      }
**
**      exit(rtc);
**   }
**
\*------------------------------------------------------------------------------*/
