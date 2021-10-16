!//---------------------------------------------------------------------
!//   MessagesSpecial.inc -- Error Message Display System (F77 and F90)
!//
!//   This include file is "special" (compared to Messages.inc) because
!//   it will work without a preprocessor.  It involves only standard
!//   Fortran and so it should work even on non-Unix systems.
!//
!//   One biggest difference in this version of the include file is
!//   that it does not support a variable number of arguments in the
!//   calls to the MESSAGES routines.  So to embed variables in a
!//   message, the programmer has to use the "buffer" method (see
!//   Example 2 below).
!//
!//   Another big difference in using this include file is that the
!//   source code file name and line number are not automatically
!//   generated and presented at the higher levels of displayed detail.
!//   Instead, the programmer has to supply the information as part of
!//   the calls to the MESSAGES routines.  And since it is very
!//   difficult for a programmer to supply the correct source code line
!//   number, that number takes on a different definition.  It is no
!//   longer treated as a line number, but rather as an arbitrary
!//   unique integer code.  Such code might be unique only to the one
!//   source code file, or it could be unique across the whole project
!//   (that decision is left to the programmer).
!//
!//   Example 1 (static message):
!//
!//  |      include "MessagesSpecial.inc"
!//  |
!//  |      CHARACTER fnm*50
!//  |      fnm = "SourceFile.f"
!//  |
!//  |      INTEGER uc123 = 123   ! such "unique codes" could all be
!//  |      INTEGER uc456 = 456   ! in a project-wide include file
!//  |      INTEGER uc078 =  78   ! where all codes (0...) are defined
!//  |
!//  |      CALL MsgErr (    fnm, uc123, "static error   message")
!//  |      CALL MsgWrn (L2, fnm, uc456, "static warning message")
!//  |      CALL MsgInf (L5, fnm, uc078, "static info    message")
!//
!//
!//   Example 2 (message with variables):
!//
!//  |      include "MessagesSpecial.inc"
!//  |
!//  |      CHARACTER fnm*50
!//  |      CHARACTER lsf*50
!//  |      INTEGER acc
!//  |      CHARACTER buf * ( EWI_LEN )
!//  |      fnm = "SourceFile.f"
!//  |      lsf = "leapsec.dat"
!//  |      acc = 440
!//  |
!//  |      INTEGER   uc123 = 123   ! such "unique codes" could all be
!//  |      INTEGER   uc456 = 456   ! in a project-wide include file
!//  |      INTEGER   uc078 =  78   ! where all codes (0...) are defined
!//  |
!//  |      WRITE (buf, '(A, A11, A, I3, A)')
!//  |     +                 "could not open <", lsf,
!//  |     +                 "> because file access is <", acc, ">"
!//  |
!//  |      CALL MsgErr (    fnm, uc123, buf)
!//  |      CALL MsgWrn (L2, fnm, uc456, buf)
!//  |      CALL MsgInf (L5, fnm, uc078, buf)
!//
!//
!//   Example 3 (set and get maximum level of displayed detail):
!//
!//  |      include "MessagesSpecial.inc"
!//  |
!//  |      INTEGER Lx
!//  |
!//  |      CALL SetMaximumLevelOfDetail(L3)
!//  |      CALL GetMaximumLevelOfDetail(Lx)
!//
!//   F90 Note: When using this include file in F90 modules, declare all
!//   symbols PRIVATE and then declare PUBLIC only those symbols that are
!//   truly needed to be visible.  For example:
!//
!//   MODULE SomethingModule	! start of module
!//   USE HDF5			! USE statements, if any, go here
!//   IMPLICIT NONE
!//   PRIVATE			! all symbols default to PRIVATE
!//   include "MessagesSpecial.inc"
!//   PUBLIC :: SubName		! declare PUBLIC only what is needed
!//   CONTAINS
!//   SUBROUTINE SubName ()	! start of subroutine or function
!//   [...]
!//   Info(L5, fnm, uc078, "some progress message to display")
!//   [...]
!//
!//---------------------------------------------------------------------

      INTEGER L0, L1, L2, L3, L4, L5, L6, L7, L8, L9

      PARAMETER ( L0 = 0 )
      PARAMETER ( L1 = 1 )
      PARAMETER ( L2 = 2 )
      PARAMETER ( L3 = 3 )
      PARAMETER ( L4 = 4 )
      PARAMETER ( L5 = 5 )
      PARAMETER ( L6 = 6 )
      PARAMETER ( L7 = 7 )
      PARAMETER ( L8 = 8 )
      PARAMETER ( L9 = 9 )

      INTEGER EWI_LEN

      PARAMETER ( EWI_LEN = 500 )

      CHARACTER EWIBUF * ( EWI_LEN + 1 )

!//---------------------------------------------------------------------
