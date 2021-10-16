!//---------------------------------------------------------------------
!//   Messages.inc -- Error Message Display System (for F77 and F90)
!//
!//   Example 1:
!//
!//     |#include "Messages.inc"
!//     |
!//     |      INTEGER Lx
!//     |
!//     |      SetMaximumLevelOfDetail(L3)
!//     |      GetMaximumLevelOfDetail(Lx)
!//     |
!//     |      Info(L0, "Level of Detail is set at ", Lx)
!//     |
!//     |      Warn(L2, "warning number: ", 17, "   and another: ", 42)
!//     |
!//     |      Info(L5, "info    number: ", 17, "   and another: ", 42)
!//     |
!//     |      Error(   "error   number: ", 17, "   and another: ", 42)
!//     |      Info(L0, "continuation of an error message")
!//
!//   F77 Note: When using this include file in F77 code, make sure that
!//   the "#include" statement begins in column 1 (while the rest of the
!//   Fortran statements begin in column 7 or beyond).
!//
!//   F90 Note: When using this include file in F90 modules, declare all
!//   symbols PRIVATE and then declare PUBLIC only those symbols that are
!//   truly needed to be visible.  For example:
!//
!//      MODULE SomethingModule	! start of module
!//      USE HDF5		! USE statements, if any, go here
!//      IMPLICIT NONE
!//      PRIVATE		! all symbols default to PRIVATE
!//      #include "Messages.inc"
!//      PUBLIC :: SubName	! declare PUBLIC only what is needed
!//      CONTAINS
!//      SUBROUTINE SubName ()	! start of subroutine or function
!//      [...]
!//      Info(L3, "some progress message to display")
!//      [...]
!//
!//---------------------------------------------------------------------

      INTEGER L0, L1, L2, L3, L4, L5, L6, L7, L8, L9, L10, L11, L12

      PARAMETER ( L0  =  0 )
      PARAMETER ( L1  =  1 )
      PARAMETER ( L2  =  2 )
      PARAMETER ( L3  =  3 )
      PARAMETER ( L4  =  4 )
      PARAMETER ( L5  =  5 )
      PARAMETER ( L6  =  6 )
      PARAMETER ( L7  =  7 )
      PARAMETER ( L8  =  8 )
      PARAMETER ( L9  =  9 )
      PARAMETER ( L10 = 10 )
      PARAMETER ( L11 = 11 )
      PARAMETER ( L12 = 12 )

      INTEGER EWI_LEN

      PARAMETER ( EWI_LEN = 500 )

      CHARACTER EWIBUF * ( EWI_LEN + 1 )

#define error(...) \
        WRITE(EWIBUF, *) __VA_ARGS__ newline\\
        CALL MsgErr(__FILE__, __LINE__, EWIBUF)
#define ERROR error
#define Error error

#define warn(Lx, ...) \
        WRITE(EWIBUF, *) __VA_ARGS__ newline\\
        CALL MsgWrn(Lx, __FILE__, __LINE__, EWIBUF)
#define WARN warn
#define Warn warn

#define info(Lx, ...) \
        WRITE(EWIBUF, *) __VA_ARGS__ newline\\
        CALL MsgInf(Lx, __FILE__, __LINE__, EWIBUF)
#define INFO info
#define Info info

#define SetMaximumLevelOfDetail(Lx) CALL SetMaximumLevelOfDetail(Lx)
#define setmaxlod SetMaximumLevelOfDetail
#define SetMaxLoD SetMaximumLevelOfDetail
#define SetMaxLOD SetMaximumLevelOfDetail
#define SETMAXLOD SetMaximumLevelOfDetail

#define GetMaximumLevelOfDetail(Lx) CALL GetMaximumLevelOfDetail(Lx)
#define getmaxlod GetMaximumLevelOfDetail
#define GetMaxLoD GetMaximumLevelOfDetail
#define GetMaxLOD GetMaximumLevelOfDetail
#define GETMAXLOD GetMaximumLevelOfDetail

!//---------------------------------------------------------------------
