/*------------------------------------------------------------------------------*\
** ymlPkgInc.h -- not for users
**
** To call any functions in this package, use "ymlForApps.h" (not this file).
\*------------------------------------------------------------------------------*/

#define _FILE_OFFSET_BITS 64	/* for 32-bit executable on a 64-bit system	*/

#include <sys/stat.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <errno.h>
#include "yaml.h"
#include "Messages.h"
#include "O3_PEATE_Common.h"
#include "ymlForApps.h"

#define Err(rco) if(rco < 0){Error(msg);return(rco);}else{Info(LDx,"OK %s",msg);}

/*------------------------------------------------------------------------------*/

extern int  mss;				/* error msg[] buffer size	*/
extern char msg[];				/* error or progress message	*/

extern int  cfx;				/* index for writing tbl[]	*/
extern int  stx;				/* index for writing stk[]	*/

/*------------------------------------------------------------------------------*/
/* prototypes for internal functions not available to users of this package	*/

int ParUseEvt(int LDx, int LDy, int nry, yTB tbl[], yaml_parser_t parser);
int ParUseTok(int LDx, int LDy, int nry, yTB tbl[], yaml_parser_t parser);
int ResetPars(int LDx, int* lvl, int* blk, int* vnx, int* knx, int* arv);
int ShowEvent(int LDx, yaml_event_type_t ytp, char* ynm);
int ShowFlags(int LDx, int lvl, int blk, int knx, int vnx, int arv);
int ShowToken(int LDx, yaml_token_type_t ytp, char* ynm);
int Show_TorE(int LDx, yaml_event_type_t ytp, char* ynm);
int ymlAddKey(int LDx, int nry, yTB tbl[], int lvl, char* stk[], int max);
int ymlAddVal(int LDx, int nry, yTB tbl[], int lvl, int blk, char* val);
int ymlChkTxt(int LDx, char buf[]);
int ymlDmpStk(int LDx, yTB tbl[]);
int ymlDmpTbl(int LDx, int nry, yTB tbl[]);
int ymlErrors(int LDx, yaml_parser_t parser);
int ymlEstNrE(int LDx, char* mem, int* nrp);
int ymlIniTbl(int LDx, int* nrp, yTB** ptr);
int ymlMaxDtl(int LDx, char* key, char* val);
int ymlNrmKey(int LDx, char a[], char b[]);
int ymlParTxt(int LDx, int toe, char* mem, int* nrp, yTB** ptr);
int ymlRdFile(int LDx, char* fnm, char** mem);
int ymlStkPsh(int LDx, char* str, int lvl, char* stk[], int max);

/*------------------------------------------------------------------------------*/
