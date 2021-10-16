!-----------------------------------------------------------------------
! GetConfig.inc -- Everything for using GetConfig_X_() functions
!
! GetConfig_X_() represents the following functions:
!
!   GetConfigString()
!   GetConfigInteger()
!   GetConfigReal()
!   GetConfigDouble()
!   GetConfigRealArr()
!   GetConfigDoubleArr()
!
! All are identical except for the parameter type that is returned
! (CHARACTER, INTEGER, REAL, DOUBLE, or an array of REALs or DOUBLEs).
!
! Parameters:
!   in  flg = specifies whether no value is deemed as an error
!   in  key = specifies what caller wants out of the control file
!   in  cnt = (only for Arr) count of elements in the array
!   out val = value returned to caller
!   out err = return code:
!           1 - there is no value to retrieve (and flg was not "E" or "e")
!           0 - everything is OK, val has retrieved value
!          -1 - key string length exceeds CFG_KEY_LEN
!          -2 - key not found in control file
!          -3 - unable to open control file
!          -4 - control file opened, read wrong number of bytes
!          -5 - control file opened, but read failed
!          -6 - control file appears empty
!          -7 - control file has TAB characters (TABs are illegal in YAML)
!          -8 - YAML parser had problem initializing
!          -9 - YAML parser had problem while parsing
!         -10 - internal overflow of stk[KEY_STACK_LEN]
!         -11 - (not used)
!         -12 - internal overflow of tmp[CFG_KEY_LEN]
!         -13 - internal overflow of cfg[].val[CFG_VAL_LEN]
!         -14 - there is no value to retrieve (and flg was "E" or "e")
!         -15 - failed to allocate memory for control file ingest
!
! A negative return code indicates an error that is generally fatal.
! A positive return code indicates an unusual condition that may not
! be an error.
!
! key and val are simple strings.
!
! The first parameter in the call to GetConfig_X_() is a 1-character
! flg that may be either "E" (or "e") or something else (including a
! null string).  When it is "E" or "e" and there is no value being
! returned by GetConfig_X_(), an error is displayed.  When it is not
! "E" or "e" and no value is being returned, no error is displayed.
!
! Each call to GetConfig_X_() writes the value, specified by key, into
! the caller's variable.  Multiple calls with the same key return subsequent
! values in the list for the key.  If the end of the list is reached, a null
! is written as the first string character instead of the value.
!
! The control file name (possibly a full path) needs to be the first
! command-line parameter when executing a program that uses GetConfig_X_().
! The control file is assumed to be in YAML format.
!
! Example 1:
!
!        INCLUDE "GetConfig.inc"
!
!        CHARACTER VAL*(CFG_VAL_LEN)
!
!        GetConfigString("E", "Some Key Name", VAL)
!
! Example 2:
!
!        INCLUDE "GetConfig.inc"
!
!        CHARACTER KEY*(CFG_KEY_LEN)
!        CHARACTER VAL*(CFG_VAL_LEN)
!
!        INTEGER IER
!
!        KEY = "Runtime Parameters Orbit Number"
!
!        IER = GetConfigString("0", KEY, VAL)
!
!        IF ( IER .LT. 0 ) THEN
!           WRITE(*, 100) IER, KEY
!        ELSE
!           IF ( VAL .EQ. "" ) THEN
!              WRITE(*, 200) KEY
!           ELSE
!              WRITE(*, 300) KEY, VAL
!           ENDIF
!        ENDIF
!
!  100   FORMAT(1X, "  Error ", I5, " getting config <", A, ">)
!  200   FORMAT(1X, "    OK: no config value for <", A, ">")
!  300   FORMAT(1X, "    OK: <", A, "> = <", A, ">")
!
!-----------------------------------------------------------------------

      !INTEGER CFG_KEY_LEN
      INTEGER CFG_VAL_LEN

      !PARAMETER ( CFG_KEY_LEN = 260 )
      PARAMETER ( CFG_VAL_LEN = 260 )

      INTEGER GetConfigString
      !INTEGER GetConfigInteger
      !INTEGER GetConfigReal
      !INTEGER GetConfigDouble
      !INTEGER GetConfigRealArr
      !INTEGER GetConfigDoubleArr
      !INTEGER MetaConfig
      !INTEGER GetConfigClose

!-----------------------------------------------------------------------
