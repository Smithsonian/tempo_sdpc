!
module m_vlidort_l_inputs

  use m_vlidort_inputs
  use m_vlidort_aux

contains

! ###############################################################
! #                                                             #
! #                    THE VLIDORT  MODEL                       #
! #                                                             #
! #  Vectorized LInearized Discrete Ordinate Radiative Transfer #
! #  -          --         -        -        -         -        #
! #                                                             #
! ###############################################################
!C ###############################################################
! #                                                             #
! #  Author :      Robert. J. D. Spurr                          #
! #                                                             #
! #  Address :      RT Solutions, Inc.                          #
! #            9 Channing Street                                #
! #             Cambridge, MA 02138, USA                        #
! #            Tel: (617) 492 1183                              #
! #                                                             #
! #  Email :      rtsolutions@verizon.net                       #
! #                                                             #
! #  Versions      :   2.0, 2.2, 2.3, 2.4                       #
! #  Release Date :   December 2005 (2.0)                       #
! #  Release Date :   March 2007    (2.2)                       #
! #  Release Date :   October 2007  (2.3)                       #
! #  Release Date :   December 2008  (2.4)                      #
! #                                                             #
! ###############################################################

!    #####################################################
!    #                                                   #
!    #   This Version of VLIDORT comes with a GNU-style  #
!    #   license. Please read the license carefully.     #
!    #                                                   #
!    #####################################################

! ###############################################################
! #                                                             #
! # Subroutines in this Module                                  #
! #                                                             #
! #            VLIDORT_L_INPUT_MASTER (master)                  #
! #            VLIDORT_L_INPUTREAD                              #
! #            VLIDORT_L_INIT_VARS                              #
! #            VLIDORT_L_CHECK_INPUT                            #
! #                                                             #
! ###############################################################

      SUBROUTINE VLIDORT_L_INPUT_MASTER ( FILNAM, EFILNAM, STATUS )

!  Read all control inputs for linearized VLIDORT
!  ----------------------------------------------

!  include file of constants and dimensions

      INCLUDE 'VLIDORT.PARS'

!  include file of INPUT variables

      INCLUDE 'VLIDORT_INPUTS.VARS'

!  Module arguments (input filename, error filename,output status)

      CHARACTER*(*)      FILNAM
      CHARACTER*(*)      EFILNAM
      INTEGER            STATUS

!  local variables

      INTEGER        STATUS_SUB, FILUNIT, LEN_STRING
      CHARACTER*70   MAIL
      EXTERNAL       LEN_STRING

!  initialize status

      STATUS = VLIDORT_SUCCESS

!  initialize error file

      VLIDORT_ERROR_FILENAME = EFILNAM
      VLIDORT_ERROR_INIT     = .TRUE.

!  initialize standard variables
 
      CALL VLIDORT_INIT_CONTROL_VARS
      CALL VLIDORT_INIT_SURFACE_VARS
      CALL VLIDORT_INIT_MODEL_VARS

!  initialize linearization variables

      CALL VLIDORT_L_INIT_VARS

!  Open file

      FILUNIT = VLIDORT_INUNIT
      OPEN(VLIDORT_INUNIT,FILE=FILNAM,ERR=300,STATUS='OLD')

!  read standard inputs

      CALL VLIDORT_INPUTREAD ( STATUS_SUB )
      IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
        STATUS = VLIDORT_SERIOUS
        MAIL = 'Error from VLIDORT_INPUTREAD'
        CALL VLIDORT_ERROR_TRACE ( MAIL, ' ', STATUS )
        CLOSE(FILUNIT)
        RETURN
      ENDIF

!  read additional inputs for linearization

      CALL VLIDORT_L_INPUTREAD ( STATUS_SUB )
      IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
        STATUS = VLIDORT_SERIOUS
        MAIL = 'Error from VLIDORT_L_INPUTREAD'
        CALL VLIDORT_ERROR_TRACE ( MAIL, ' ', STATUS )
        CLOSE(FILUNIT)
        RETURN
      ENDIF

!  normal return

      CLOSE(FILUNIT)
      RETURN
      
!  Open file error

300   CONTINUE
      STATUS = VLIDORT_SERIOUS
      MAIL = 'openfile failure for '//FILNAM(1:LEN_STRING(FILNAM))
      CALL VLIDORT_ERROR_TRACE ( MAIL, ' ', STATUS )
      CLOSE(FILUNIT)

!  Finish

      RETURN
      END subroutine VLIDORT_L_INPUT_MASTER

!

      SUBROUTINE VLIDORT_L_INPUTREAD ( STATUS )

!  Read all linearization control inputs for VLIDORT
!  ------------------------------------------------

!  include file of constants and dimensions

      INCLUDE 'VLIDORT.PARS'

!  include files of input and bookkeeping variables

      INCLUDE 'VLIDORT_INPUTS.VARS'
      INCLUDE 'VLIDORT_BOOKKEEP.VARS'

!  include file to be filled out with input data

      INCLUDE 'VLIDORT_L_INPUTS.VARS'

!  output status

      INTEGER         STATUS

!  local variables

      CHARACTER*9     PREFIX
      PARAMETER       ( PREFIX = 'VLIDORT -' )
      LOGICAL         ERROR
      CHARACTER * 80  PAR_STR
      LOGICAL         GFINDPAR
      EXTERNAL        GFINDPAR
      INTEGER         FILUNIT, LEN_STRING, I, J
      INTEGER         DUM_INDEX, DUM_NPARS, NC
      CHARACTER*10    DUM_NAME
      CHARACTER*70    MAIL
      EXTERNAL        LEN_STRING

!  initialize status

      STATUS = VLIDORT_SUCCESS
      ERROR  = .FALSE.

!  file unit

      FILUNIT = VLIDORT_INUNIT

!  linearization control

      PAR_STR = 'Simulation only?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) DO_SIMULATION_ONLY
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

      PAR_STR = 'Atmospheric profile weighting functions?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) DO_PROFILE_LINEARIZATION
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

      PAR_STR = 'Atmospheric column weighting functions?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) DO_COLUMN_LINEARIZATION
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

      PAR_STR = 'Surface property weighting functions?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) DO_SURFACE_LINEARIZATION
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

!  disabled
!      PAR_STR = 'Surface emission weighting functions?'
!      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
!     &     READ (FILUNIT,*,ERR=998) DO_SURFBB_LINEARIZATION
!      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

!  Atmospheric profile weighting function input control
!  individual layer number does not exceed whole number

      IF ( DO_PROFILE_LINEARIZATION ) THEN

        PAR_STR='Atmospheric profile weighting functions, total number'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR))
     &    READ (FILUNIT,*,ERR=998) N_TOTALPROFILE_WFS
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

!  This code should be prepared, not taken from file-read
!        PAR_STR = 'Atmospheric weighting functions, layer control'
!        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
!          DO I = 1, NLAYERS
!            READ (FILUNIT,*,ERR=998) LAYER_VARY_FLAG(I),
!     &                               LAYER_VARY_NUMBER(I)
!          ENDDO
!        ENDIF
!        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

      ENDIF

!  atmospheric Bulk/column weighting function input control

       IF ( DO_COLUMN_LINEARIZATION ) THEN

         PAR_STR='Atmospheric column weighting functions, total number'
         IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR))
     &    READ (FILUNIT,*,ERR=998) N_TOTALCOLUMN_WFS
         CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

         DO I = 1, NLAYERS
           LAYER_VARY_FLAG(I)   = .TRUE.
           LAYER_VARY_NUMBER(I) = N_TOTALCOLUMN_WFS
         ENDDO

       ENDIF

!  Surface weighting function input control

      IF ( DO_SURFACE_LINEARIZATION ) THEN

!  Lambertian WF flag

       IF ( DO_LAMBERTIAN_SURFACE ) THEN
        PAR_STR = 'Do Lambertian albedo weighting function?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          READ (FILUNIT,*,ERR=998)DO_KERNEL_FACTOR_WFS(1)
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )
        N_TOTALBRDF_WFS = 1
        DO_KERNEL_FACTOR_WFS(1) = .TRUE.
       ENDIF

!  kernel input. Names should agree with earlier entries.

       IF ( .NOT. DO_LAMBERTIAN_SURFACE ) THEN
        PAR_STR = 'Kernels, indices, # pars, Jacobian flags'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, N_BRDF_KERNELS
            READ (FILUNIT,57,ERR=998)
     &       DUM_NAME, DUM_INDEX,DUM_NPARS,DO_KERNEL_FACTOR_WFS(I),
     &           (DO_KERNEL_PARAMS_WFS(I,J),J=1,3)

            IF ( DUM_NAME .NE. BRDF_NAMES(I) ) THEN
              MAIL = 'Input BRDF name not same as earlier list'
              STATUS  = VLIDORT_SERIOUS
              CALL VLIDORT_ERROR_TRACE ( MAIL, ' ', STATUS )
              RETURN
            ENDIF

            IF ( DUM_INDEX .NE. WHICH_BRDF(I) ) THEN
              MAIL = 'Input BRDF index not same as earlier list'
              STATUS  = VLIDORT_SERIOUS
              CALL VLIDORT_ERROR_TRACE ( MAIL, ' ', STATUS )
              RETURN
            ENDIF

            IF ( DUM_NPARS .NE. N_BRDF_PARAMETERS(I) ) THEN
              MAIL = 'Input BRDF parameter no. not same as earlier'
              STATUS  = VLIDORT_SERIOUS
              CALL VLIDORT_ERROR_TRACE ( MAIL, ' ', STATUS )
              RETURN
            ENDIF

!  Compute total number of pars

            IF ( DO_KERNEL_FACTOR_WFS(I) ) NC = NC + 1
            DO J = 1, N_BRDF_PARAMETERS(I)
              IF ( DO_KERNEL_PARAMS_WFS(I,J) ) NC = NC + 1
            ENDDO

          ENDDO
          N_TOTALBRDF_WFS = NC
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )
 57     FORMAT( A10, I3, I2, 1X, L2, 2X, 3L2 ) 

!  Check total number of BRDF weighting functions is not out of bounds

        IF ( N_TOTALBRDF_WFS .GT. MAX_SURFACEWFS ) THEN
         MAIL  = 'Bad error: Insuffient dimensioning for Surface WFs'
         STATUS = VLIDORT_SERIOUS
         CALL VLIDORT_ERROR_TRACE ( MAIL, ' ', STATUS )
         RETURN
        ENDIF
   
       ENDIF
      ENDIF

!  weighting function names

      IF ( DO_PROFILE_LINEARIZATION ) THEN
        PAR_STR='Atmospheric profile Jacobian names (character*31)'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, N_TOTALPROFILE_WFS
            READ (FILUNIT,'(a31)',ERR=998) PROFILEWF_NAMES(I)
          ENDDO
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )
      ENDIF

      IF ( DO_COLUMN_LINEARIZATION ) THEN
        PAR_STR='Atmospheric column Jacobian names (character*31)'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, N_TOTALCOLUMN_WFS
            READ (FILUNIT,'(a31)',ERR=998) COLUMNWF_NAMES(I)
          ENDDO
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )
      ENDIF

      IF ( DO_SURFACE_LINEARIZATION ) THEN
       IF ( .NOT. DO_LAMBERTIAN_SURFACE ) THEN
        PAR_STR = 'Surface Jacobian names (character*22)'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, N_TOTALBRDF_WFS
            READ (FILUNIT,'(a22)',ERR=998) SURFACEWF_NAMES(I)
          ENDDO
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )
       ELSE
        SURFACEWF_NAMES(1) = 'Lambertian Amplitude  '
       ENDIF
      ENDIF

!  weighting function names

!      IF ( DO_ATMOS_LINEARIZATION ) THEN
!        PAR_STR='Atmospheric weighting function names (character*31)'
!        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
!          DO I = 1, N_TOTALATMOS_WFS
!            READ (FILUNIT,'(a31)',ERR=998) ATMOSWF_NAMES(I)
!          ENDDO
!        ENDIF
!        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )
!      ENDIF

!      IF ( DO_SURFACE_LINEARIZATION ) THEN
!        PAR_STR = 'Surface weighting function names (character*31)'
!        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
!          DO I = 1, N_TOTALBRDF_WFS
!            READ (FILUNIT,'(a31)',ERR=998) SURFACEWF_NAMES(I)
!          ENDDO
!        ENDIF
!        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )
!      ENDIF

!  set overall linearization flags

      DO_ATMOS_LINEARIZATION = ( DO_PROFILE_LINEARIZATION.OR.
     &                           DO_COLUMN_LINEARIZATION )

      DO_LINEARIZATION = ( DO_ATMOS_LINEARIZATION   .OR.
     &                     DO_SURFACE_LINEARIZATION .OR.
     &                     DO_SURFBB_LINEARIZATION )

!  normal return

      RETURN

!  line read error - abort immediately

998   CONTINUE
      STATUS = VLIDORT_SERIOUS
      MAIL = 'read failure for '//PAR_STR(1:LEN_STRING(PAR_STR))
      CALL VLIDORT_ERROR_TRACE ( MAIL, ' ', STATUS )

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_L_INPUTREAD

!

      SUBROUTINE VLIDORT_L_INIT_VARS

!  Initialises all linearization control inputs for VLIDORT
!  --------------------------------------------------------

!  include file of constants and dimensions

      INCLUDE 'VLIDORT.PARS'

!  include file to be initialised

      INCLUDE 'VLIDORT_L_INPUTS.VARS'

!  local variables

      INTEGER             B, Q, N

!  initialise overall control inputs

      DO_SIMULATION_ONLY        = .FALSE.
      DO_LINEARIZATION          = .FALSE.

!  initialise atmospheric weighting function control inputs

      DO_PROFILE_LINEARIZATION  = .FALSE.
      DO_COLUMN_LINEARIZATION   = .FALSE.
      DO_ATMOS_LINEARIZATION    = .FALSE.
      N_TOTALCOLUMN_WFS  = 0
      N_TOTALPROFILE_WFS = 0
      DO N = 1, MAXLAYERS
        LAYER_VARY_NUMBER(N) = 0
        LAYER_VARY_FLAG(N)   = .FALSE.
      ENDDO

!  initialise surface weighting function control inputs

      DO_SURFACE_LINEARIZATION = .FALSE.
      DO_SURFBB_LINEARIZATION  = .FALSE.
      N_TOTALBRDF_WFS  = 0
      DO B = 1, MAX_BRDF_KERNELS
        DO_KERNEL_FACTOR_WFS(B) = .FALSE.
        DO Q = 1, MAX_BRDF_PARAMETERS
          DO_KERNEL_PARAMS_WFS(B,Q) = .FALSE.
        ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_L_INIT_VARS

!

      SUBROUTINE VLIDORT_L_CHECK_INPUT (STATUS)

!  include files
!  -------------

!  include file of constants and dimensions

      INCLUDE 'VLIDORT.PARS'

!  Include file of input variables

      INCLUDE 'VLIDORT_INPUTS.VARS'

!  include file with variables to be checked out

      INCLUDE 'VLIDORT_L_INPUTS.VARS'

!  Module output
!  -------------

      INTEGER           STATUS

!  local variables
!  ---------------

      CHARACTER*70      MAIL, TRACE
      CHARACTER*1       NOTRACE
      PARAMETER         ( NOTRACE = ' ')

!  Initialize output status

      STATUS = VLIDORT_SUCCESS

!  Do some checking
!  ================

!  Check formatting; should not exceed 10 entries for the number of
!  layer weighting functions. This is needed to avoid errors with
!  compilers not using variable formatting.

!      IF ( MAX_PARAMETERS .GT. 10 ) THEN
!        MAIL='Max layer variations > 10, get formatting problems!'
!        STATUS = VLIDORT_SERIOUS
!        CALL VLIDORT_ERROR_TRACE ( MAIL, NOTRACE, STATUS )
!      ENDIF

!  Check something is being varied

      IF ( .NOT. DO_SIMULATION_ONLY ) THEN
        IF ( .NOT. DO_ATMOS_LINEARIZATION .AND.
     &         .NOT. DO_SURFACE_LINEARIZATION. AND.
     &         .NOT. DO_SURFBB_LINEARIZATION  ) THEN
          MAIL='Bad input: Nothing being varied or calculated'
          STATUS = VLIDORT_SERIOUS
          CALL VLIDORT_ERROR_TRACE ( MAIL, NOTRACE, STATUS )
        ENDIF
      ENDIF

!  Check and fix AF flags if you want simulation only

      IF ( DO_SIMULATION_ONLY ) THEN
        IF ( DO_ATMOS_LINEARIZATION ) THEN
          MAIL  = 'Input warning; simulation only ---> NO Atmos WFs'
          TRACE = 'Action: LIDORT switched off WF flag for you!'
          STATUS = VLIDORT_WARNING
          CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
          DO_PROFILE_LINEARIZATION = .FALSE.
          DO_COLUMN_LINEARIZATION  = .FALSE.
          DO_ATMOS_LINEARIZATION   = .FALSE.
        ENDIF
        IF ( DO_SURFACE_LINEARIZATION ) THEN
          MAIL  = 'Input warning; simulation only ---> NO Surface WFs'
          TRACE = 'Action: LIDORT switched off WF flag for you!'
          STATUS = VLIDORT_WARNING
          CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
          DO_SURFACE_LINEARIZATION = .FALSE.
        ENDIF
      ENDIF

!  Check total number of BRDF weighting functions is not out of bounds

      IF ( DO_SURFACE_LINEARIZATION ) THEN
       IF ( .NOT.DO_LAMBERTIAN_SURFACE ) THEN
        IF ( N_TOTALBRDF_WFS .GT. MAX_SURFACEWFS ) THEN
         MAIL  = 'Bad error: Insuffient dimensioning for Surface WFs'
         TRACE = 'Action: increase MAX_SURFACE_WFS in LIDORT.PARS!'
         STATUS = VLIDORT_SERIOUS
         CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
        ENDIF
       ENDIF
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_L_CHECK_INPUT

end module m_vlidort_l_inputs
