MODULE error_module

  USE parameters_module

  IMPLICIT NONE

  TYPE ErrorType
    LOGICAL                :: DoDebug
    LOGICAL                :: Fatal
    LOGICAL                :: Pixel
    INTEGER                :: xPixelIndex
    INTEGER                :: yPixelIndex
    INTEGER(KIND=2)        :: Code
    INTEGER(KIND=8)        :: StartTime
    CHARACTER(LEN=maxChar) :: LoggerFileName
    INTEGER                :: VerbosityLevel

    !CHARACTER(LEN=maxChar) :: message(maxChar)
    INTEGER                :: nWarning
    INTEGER                :: nMessage
    INTEGER                :: nPixelError
    INTEGER                :: LoggerUnit
  ENDTYPE ErrorType
  
  PUBLIC  :: InitErrorLogger,   &
             RaiseFatalError,   &
             RaisePixelError,   &
             ResetPixelError,   &
             CheckError
  
  CONTAINS

  SUBROUTINE InitErrorLogger( DoDebug, ErrorFileName, VerbosityLevel, Error )

    ! --------------------
    ! Subroutine Arguments
    ! --------------------
    LOGICAL,          INTENT(IN)    :: DoDebug
    CHARACTER(LEN=*), INTENT(IN)    :: ErrorFileName
    INTEGER,          INTENT(IN)    :: VerbosityLevel
    TYPE(ErrorType),  INTENT(INOUT) :: Error

    ! ---------------
    ! local variables 
    ! ---------------
    CHARACTER(LEN=8)  :: caldate
    CHARACTER(LEN=10) :: caltime
    CHARACTER(LEN=5)  :: calzone

    ! =====================================================================
    ! InitErrorLogger Starts here
    ! =====================================================================

    ! Save the logger file name
    Error%LoggerFileName = ErrorFileName
    Error%VerbosityLevel = VerbosityLevel
    Error%DoDebug           = DoDebug
    Error%Fatal             = .FALSE.
    Error%Pixel             = .FALSE.
    Error%xPixelIndex       = 0
    Error%yPixelIndex       = 0
    Error%code              = 0
    Error%nWarning          = 0
    Error%nMessage          = 0
    Error%nPixelError       = 0
    Error%StartTime         = 0.0 !jbak time() ! s since Jan 1, 1970
    
    
    IF(Error%DoDebug) THEN

      ! Get Calender date/time
      CALL date_and_time(DATE=caldate,TIME=caltime,ZONE=calzone)

      ! Set Logger Unit
      Error%LoggerUnit = errunit

      ! Open Error Log
      OPEN(UNIT=errunit,                     &
           FILE=TRIM(ADJUSTL(ErrorFileName)),&
           STATUS='REPLACE',                 &
           ACTION='WRITE'                    )

      ! Write Header to Error Log
      WRITE(Error%LoggerUnit,'(A)') '########################################################'
      WRITE(Error%LoggerUnit,'(A)') '#               SPLAT DEBUG ERROR LOGGER               #'
      WRITE(Error%LoggerUnit,'(A)') '########################################################'
      WRITE(Error%LoggerUnit,'(A)')
      WRITE(Error%LoggerUnit,'(A)') 'LOG Start: ' // caldate // ' ' // caltime // ' UTC'// calzone
      WRITE(Error%LoggerUnit,'(A)')

    ELSE

      ! If the logger is not set print warning to screen
      Error%LoggerUnit = 6

    ENDIF
    
    ! Test the fatal error logger
    ! CALL RaisePixelError( Error, 1, 'error_module', 'InitErrorLogger',  &
    !                      'Raising pixel error to test logger', 'nothing')
    ! CALL RaiseWarning( Error, 1, 'error_module', 'InitErrorLogger',  &
    !                   'Raising warning  to test logger', 'nothing')
    ! CALL RaiseFatalError( Error, 1, 'error_module', 'InitErrorLogger',  &
    !                      'Raising fatal error to test logger', 'nothing')

  END SUBROUTINE InitErrorLogger
  
  SUBROUTINE WriteLoggerMessage( Error, code, ModuleName, SubroutineName, &
                                 ThisTime, Header, Message, Action        )

    ! --------------------
    ! Subroutine Arguments
    ! --------------------
    TYPE(ErrorType),            INTENT(INOUT) :: error
    INTEGER(KIND=2),            INTENT(IN)    :: code
    CHARACTER(LEN=*),           INTENT(IN)    :: ModuleName
    CHARACTER(LEN=*),           INTENT(IN)    :: SubroutineName
    INTEGER(KIND=8),            INTENT(IN)    :: ThisTime
    CHARACTER(LEN=*),           INTENT(IN)    :: Header
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN)    :: Message
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN)    :: Action

    ! ---------------
    ! Local Variables
    ! ---------------
    CHARACTER(LEN=maxChar) :: nMessStr, TimeStr, CodeStr, xPixStr, yPixStr
    REAL(KIND=8)           :: ThisTime_min

    ! =====================================================================
    ! WriteLoggerMessage Starts here
    ! =====================================================================

    ! Increment the # of messages
    Error%nMessage = Error%nMessage + 1

    ! Time in Minutes since start of simulation
    ThisTime_min = REAL(ThisTime-Error%StartTime,KIND=8) / 60.0d0
    
    ! Do some string processing
    WRITE(nMessStr,'(I100)') Error%nMessage
    WRITE(CodeStr,'(I100)') Code
    WRITE(TimeStr,'(F16.3)') ThisTime_min
    WRITE(xPixStr,'(I100)') Error%xPixelIndex
    WRITE(yPixStr,'(I100)') Error%yPixelIndex

    ! Write Error Information to logger
    WRITE(Error%LoggerUnit,'(A)') '========================================================='
    WRITE(Error%LoggerUnit,'(A)') 'MESSAGE NUMBER:' // TRIM(ADJUSTL(nMessStr))
    WRITE(Error%LoggerUnit,'(A)') TRIM(ADJUSTL(Header))
    WRITE(Error%LoggerUnit,'(A)') '  --->    MODULE: ' // TRIM(ADJUSTL(ModuleName))
    WRITE(Error%LoggerUnit,'(A)') '  --->SUBROUTINE: ' // TRIM(ADJUSTL(SubroutineName))
    WRITE(Error%LoggerUnit,'(A)') '  ---> TIME[min]: ' // TRIM(ADJUSTL(TimeStr))
    WRITE(Error%LoggerUnit,'(A)') '  --->      CODE: ' // TRIM(ADJUSTL(CodeStr))
    WRITE(Error%LoggerUnit,'(A)') '  --->PIXEL[X,Y]: [' // TRIM(ADJUSTL(xPixStr)) // ',' &
                                             // TRIM(ADJUSTL(yPixStr)) // ']'
    WRITE(Error%LoggerUnit,'(A)') '  ===>   MESSAGE: ' // TRIM(ADJUSTL(Message))
    WRITE(Error%LoggerUnit,'(A)') '  ===>    ACTION: ' // TRIM(ADJUSTL(Action))
    WRITE(Error%LoggerUnit,'(A)')

  END SUBROUTINE WriteLoggerMessage

  SUBROUTINE RaiseFatalError( Error, code, ModuleName, SubroutineName, &
                              Message_in, Action_in )

    ! --------------------
    ! Subroutine Arguments
    ! --------------------
    TYPE(ErrorType),            INTENT(INOUT) :: error
    INTEGER(KIND=2),            INTENT(IN)    :: code
    CHARACTER(LEN=*),           INTENT(IN)    :: ModuleName
    CHARACTER(LEN=*),           INTENT(IN)    :: SubroutineName
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN)    :: Message_in
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN)    :: Action_in

    ! ---------------
    ! Local Variables
    ! ---------------
    CHARACTER(LEN=maxChar) :: Message, Action, Header
    INTEGER(KIND=8)        :: ThisTime

    ! =====================================================================
    ! RaiseFatalError Starts here
    ! =====================================================================

    ! Get time 
    ThisTime = 0.0 !jbak time()
    
    ! Raise error flag/numerical code
    Error%Fatal = .TRUE.
    Error%Code  = Code

    ! Set Defaults
    Message = '' ; IF(PRESENT(Message_in)) Message = Message_in
    Action = ''  ; IF(PRESENT(Action_in) ) Action  = Action_in

    ! Write Message to log
    Header = '=====================> FATAL ERROR <====================='
    CALL WriteLoggerMessage( Error, code, ModuleName, SubroutineName, &
                             ThisTime, Header, Message, Action        )

    ! Close Logger
    CALL CloseErrorLogger(Error)

    ! Print message to screen
    print*,'======> FATAL ERROR HAS BEEN RAISED <====='
    print*,'Check the error log file for more information'
    print*,'->',TRIM(ADJUSTL(Error%LoggerFileName))
    print*,''
    print*,'Here is the traceback...'

    ! Print traceback
    CALL ABORT()

  END SUBROUTINE RaiseFatalError

  SUBROUTINE RaisePixelError(Error, Code, ModuleName, SubroutineName, &
                             Message_in, Action_in                    )
    
    ! --------------------
    ! Subroutine Arguments
    ! --------------------
    TYPE(ErrorType),            INTENT(INOUT) :: error
    INTEGER(KIND=2),            INTENT(IN)    :: code
    CHARACTER(LEN=*),           INTENT(IN)    :: ModuleName
    CHARACTER(LEN=*),           INTENT(IN)    :: SubroutineName
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN)    :: Message_in
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN)    :: Action_in

    ! ---------------
    ! Local Variables
    ! ---------------
    CHARACTER(LEN=maxChar) :: Message, Action, Header
    INTEGER(KIND=8)        :: ThisTime
    
    ! =====================================================================
    ! RaisePixelError Starts here
    ! =====================================================================
    
    ! Raise error flag/numerical code
    Error%Pixel = .TRUE.
    Error%code  = code
    
    ! Increment Number of pixel errors
    Error%nPixelError = Error%nPixelError + 1

    ! Get time 
    ThisTime = 0.0 !jbak time()
    
    ! Set Defaults
    Message = '' ; IF(PRESENT(Message_in)) Message = Message_in
    Action = ''  ; IF(PRESENT(Action_in) ) Action  = Action_in

    ! Set Header
    Header = '=====================> PIXEL ERROR <====================='
    CALL WriteLoggerMessage( Error, code, ModuleName, SubroutineName, &
                             ThisTime, Header, Message, Action        )

  END SUBROUTINE RaisePixelError
  
  
  ! Case for normal skipping e.g. Too High SZA - we dont want to log an error
  ! but want to use the built in facility to skip to next simulation
  SUBROUTINE SkipPixel(Error)
    
    ! --------------------
    ! Subroutine Arguments
    ! --------------------
    TYPE(ErrorType),            INTENT(INOUT) :: Error
    
    ! =====================================================================
    ! SkipPixel Starts here
    ! =====================================================================
    
    ! Raise error flag/numerical code
    Error%Pixel = .TRUE.
    
  END SUBROUTINE SkipPixel
  
  SUBROUTINE ResetPixelError( Error, xPixel, yPixel )
    
    ! --------------------
    ! Subroutine Arguments
    ! --------------------
    TYPE(ErrorType),  INTENT(INOUT) :: Error
    INTEGER,          INTENT(IN)    :: xPixel, yPixel
    ! =====================================================================
    ! ResetPixelError Starts here
    ! =====================================================================
    
    ! Raise error flag/numerical code
    Error%Pixel = .FALSE.
    Error%code  = 0

    ! Save the new pixel
    Error%xPixelIndex = xPixel
    Error%yPixelIndex = yPixel
    
  END SUBROUTINE ResetPixelError
  
  SUBROUTINE RaiseWarning( Error, Code, ModuleName, SubroutineName, &
                           Message_in, Action_in )
    
    ! --------------------
    ! Subroutine Arguments
    ! --------------------
    TYPE(ErrorType),            INTENT(INOUT) :: error
    INTEGER(KIND=2),            INTENT(IN)    :: code
    CHARACTER(LEN=*),           INTENT(IN)    :: ModuleName
    CHARACTER(LEN=*),           INTENT(IN)    :: SubroutineName
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN)    :: Message_in
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN)    :: Action_in

    ! ---------------
    ! Local Variables
    ! ---------------
    CHARACTER(LEN=maxChar) :: Message, Action, Header
    INTEGER(KIND=8)        :: ThisTime
    
    ! =====================================================================
    ! RaiseWarning Starts here
    ! =====================================================================
    
    ! Get time 
    ThisTime = 0.0 !jbak time()
    
    ! Set Defaults
    Message = '' ; IF(PRESENT(Message_in)) Message = Message_in
    Action = ''  ; IF(PRESENT(Action_in) ) Action  = Action_in

    ! Increment Number of warnings
    Error%nWarning = Error%nWarning + 1

    ! Set Header
    Header = '=======================> WARNING <======================='
    CALL WriteLoggerMessage( Error, Code, ModuleName, SubroutineName, &
                             ThisTime, Header, Message, Action        )

  END SUBROUTINE RaiseWarning

  LOGICAL FUNCTION CheckError( Error )
    
    ! --------------------
    ! Subroutine Arguments
    ! --------------------
    TYPE(ErrorType),  INTENT(INOUT) :: Error

    ! =====================================================================
    ! CheckError Starts here
    ! =====================================================================

    ! Placeholder for now - stop for everything
    IF( error%Fatal ) STOP 1
    
    ! Set Return value to true if there is a pixel error
    CheckError = Error%Pixel
    
    RETURN
    
  END FUNCTION CheckError

  SUBROUTINE CloseErrorLogger(Error)

    ! --------------------
    ! Subroutine Arguments
    ! --------------------
    TYPE(ErrorType),  INTENT(IN) :: Error

    ! ---------------
    ! local variables
    ! ---------------
    INTEGER(KIND=8)        :: ThisTime
    REAL(KIND=8)           :: ThisTime_min
    CHARACTER(LEN=maxChar) :: TimeStr, PixStr,WarnStr
    CHARACTER(LEN=1)       :: FatalStr

    ! =====================================================================
    ! CloseErrorLogger Starts here
    ! =====================================================================

    ! Get time 
    ThisTime = 0.0 !jbak time() 
     ThisTime_min = REAL(ThisTime-Error%StartTime,KIND=8)/60.0d0

    ! Process some strings
    WRITE(TimeStr,'(F16.3)') ThisTime_min
    FatalStr = 'F' ; IF(Error%Fatal) FatalStr = 'T'
    WRITE(WarnStr,'(I100)') Error%nWarning
    WRITE(PixStr,'(I100)') Error%nPixelError

    ! Just need to close the log file
    IF(Error%DoDebug) THEN

      ! Write Error Summary Error Information
      WRITE(Error%LoggerUnit,'(A)') '########################################################' 
      WRITE(Error%LoggerUnit,'(A)') '===================> ERROR SUMMARY <===================='
      WRITE(Error%LoggerUnit,'(A)') '########################################################' 
      WRITE(Error%LoggerUnit,'(A)') '  --->         EXIT TIME[min]: '// TRIM(ADJUSTL(TimeStr))
      WRITE(Error%LoggerUnit,'(A)') '  --->             FATAL EXIT: '// FatalStr
      WRITE(Error%LoggerUnit,'(A)') '  --->     NUMBER OF WARNINGS: '// TRIM(ADJUSTL(WarnStr))
      WRITE(Error%LoggerUnit,'(A)') '  ---> NUMBER OF PIXEL ERRORS: '// TRIM(ADJUSTL(PixStr))
      WRITE(Error%LoggerUnit,'(A)') '########################################################' 

      ! Close file
      CLOSE(Error%LoggerUnit)

    ENDIF

  END SUBROUTINE CloseErrorLogger

END MODULE error_module
