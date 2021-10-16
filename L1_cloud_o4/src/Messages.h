/*------------------------------------------------------------------------------*\
* Messages.h -- Error Message Display System (for C)
*
* Example 1:
*
*   #include "Messages.h"
*
*   int Lx;
*
*   SetMaximumLevelOfDetail(L3);
*   GetMaximumLevelOfDetail(&Lx);
*
*   Info(L0, "Level of Detail is set at %d", Lx);
*
*   Warn(L2, "warning number: %5d   and another: %d", 17, 42);
*
*   Info(L5, "info    number: %5d   and another: %d", 17, 42);
*
*   Error(   "error   number: %5d   and another: %d", 17, 42);
*   Info(L0, "continuation of an %s message", "error");
*
\*------------------------------------------------------------------------------*/

enum LevelOfDetail { L0=0, L1, L2, L3, L4, L5, L6, L7, L8, L9, L10, L11, L12 };
  
#define EWI_LEN 500

extern char ewibuf[EWI_LEN];

#define error(...) \
	snprintf(ewibuf, sizeof(ewibuf)-1, __VA_ARGS__); \
	MsgErr(__FILE__, __LINE__, ewibuf);
#define ERROR error
#define Error error

#define warn(Lx, ...) \
	snprintf(ewibuf, sizeof(ewibuf)-1, __VA_ARGS__); \
	MsgWrn(Lx, __FILE__, __LINE__, ewibuf);
#define WARN warn
#define Warn warn

#define info(Lx, ...) \
	snprintf(ewibuf, sizeof(ewibuf)-1, __VA_ARGS__); \
	MsgInf(Lx, __FILE__, __LINE__, ewibuf);
#define INFO info
#define Info info

int               SetMaximumLevelOfDetail(int  Lx);
#define setmaxlod SetMaximumLevelOfDetail
#define SetMaxLoD SetMaximumLevelOfDetail
#define SetMaxLOD SetMaximumLevelOfDetail
#define SETMAXLOD SetMaximumLevelOfDetail

int               GetMaximumLevelOfDetail(int *Lx);
#define getmaxlod GetMaximumLevelOfDetail
#define GetMaxLoD GetMaximumLevelOfDetail
#define GetMaxLOD GetMaximumLevelOfDetail
#define GETMAXLOD GetMaximumLevelOfDetail

/*------------------------------------------------------------------------------*/
