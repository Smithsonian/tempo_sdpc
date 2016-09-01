!
module m_vlidort_brdf_supplement


contains

! ###############################################################
! #                                                             #
! #                    THE VLIDORT  MODEL                       #
! #                                                             #
! #  Vectorized LInearized Discrete Ordinate Radiative Transfer #
! #  -          --         -        -        -         -        #
! #                                                             #
! ###############################################################

! ###############################################################
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
! #  Versions     :   2.0, 2.2, 2.3, 2.4                        #
! #  Release Date :   December 2005 (2.0)                       #
! #  Release Date :   March 2007    (2.2)                       #
! #  Release Date :   October 2007  (2.3)                       #
! #  Release Date :   December 2008  (2.4)                      #
! #                                                             #
! #       NEW: TOTAL COLUMN JACOBIANS         (2.4)             #
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
! # Subroutines in this Module. The "BRDF SUPPLEMENT"           #
! #                                                             #
! #   Master routines:                                          #
! #                                                             #
! #       VLIDORT_BRDF_CREATOR        (master)                  #
! #       VLIDORT_BRDF_FRTERMS        (master)                  #
! #       VLIDORT_DBCORRECTION        (master)                  #
! #                                                             #
! #    Basic BRDF setup, called by BRDF Creator:                #
! #                                                             #
! #            VLIDORT_BRDF_MAKER                               #
! #            BRDF_QUADRATURE                                  #
! #            BRDF_QUADRATURE_TPZ                              #
! #                                                             #
! #    External calls in VLIDORT_FOURIER_BRDF:                  #
! #                                                             #
! #              BRDF_SURFACE_DIRECTBEAM                        #
! #              BVP_BRDF_SURFACE_HOM                           #
! #              BVP_BRDF_SURFACE_BEAM                          #
! #              BOA_BRDF_SOURCE                                #
! #                                                             #
! ###############################################################

      SUBROUTINE VLIDORT_BRDF_CREATOR 

!  Prepares the bidirectional reflectance functions
!  necessary for VLIDORT.

!  Include files
!  -------------

!  include file of dimensions and numbers

      INCLUDE 'VLIDORT.PARS'

!  Include files of input and bookkeeping variables 

      INCLUDE 'VLIDORT_INPUTS.VARS'
      INCLUDE 'VLIDORT_BOOKKEEP.VARS'

!  Input arguments controlling type of surface
!  output arguments also go in here.

      INCLUDE 'VLIDORT_REFLECTANCE.VARS'

!  BRDF functions
!  --------------

      EXTERNAL       ROSSTHIN_VFUNCTION
      EXTERNAL       ROSSTHICK_VFUNCTION
      EXTERNAL       LISPARSE_VFUNCTION
      EXTERNAL       LIDENSE_VFUNCTION
      EXTERNAL       HAPKE_VFUNCTION
      EXTERNAL       ROUJEAN_VFUNCTION
      EXTERNAL       RAHMAN_VFUNCTION
      EXTERNAL       COXMUNK_VFUNCTION
      EXTERNAL       COXMUNK_VFUNCTION_DB
      EXTERNAL       GISSCOXMUNK_VFUNCTION
      EXTERNAL       GISSCOXMUNK_VFUNCTION_DB

!  Two new Land BRDFs, introduced 1 July 2008

      EXTERNAL       RHERMAN_VFUNCTION
      EXTERNAL       BREON_VFUNCTION

!  Hapke old uses exact DISORT code
!      EXTERNAL       HAPKE_FUNCTION_OLD

!  local arguments

      INTEGER          K, N
      INTEGER          LOCAL_BRDF_NPARS
      DOUBLE PRECISION LOCAL_BRDF_PARS ( MAX_BRDF_PARAMETERS )

!  BRDF quadrature
!  ---------------

!  Save these quantities for efficient coding

      CALL BRDF_QUADRATURE

!  Fill BRDF arrays
!  ----------------

      DO K = 1, N_BRDF_KERNELS

!  Copy parameter variables into local quantities

        LOCAL_BRDF_NPARS = N_BRDF_PARAMETERS(K)
        DO N = 1, MAX_BRDF_PARAMETERS
          LOCAL_BRDF_PARS(N) = BRDF_PARAMETERS(K,N)
        ENDDO

!  Ross thin kernel, (0 free parameters)

        IF ( WHICH_BRDF(K) .EQ. ROSSTHIN_IDX ) THEN
          CALL VLIDORT_BRDF_MAKER
     I          ( K, ROSSTHIN_VFUNCTION, ROSSTHIN_VFUNCTION,
     I            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS )
        ENDIF

!  Ross thick kernel, (0 free parameters)

        IF ( WHICH_BRDF(K) .EQ. ROSSTHICK_IDX ) THEN
          CALL VLIDORT_BRDF_MAKER
     I          ( K, ROSSTHICK_VFUNCTION, ROSSTHICK_VFUNCTION,
     I            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS )
        ENDIF

!  Li Sparse kernel; 2 free parameters

        IF ( WHICH_BRDF(K) .EQ. LISPARSE_IDX ) THEN
          CALL VLIDORT_BRDF_MAKER
     I          ( K, LISPARSE_VFUNCTION, LISPARSE_VFUNCTION,
     I            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS )
        ENDIF

!  Li Dense kernel; 2 free parameters

        IF ( WHICH_BRDF(K) .EQ. LIDENSE_IDX ) THEN
          CALL VLIDORT_BRDF_MAKER
     I          ( K, LIDENSE_VFUNCTION, LIDENSE_VFUNCTION,
     I            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS )
        ENDIF

!  Hapke kernel (3 free parameters)

        IF ( WHICH_BRDF(K) .EQ. HAPKE_IDX ) THEN
          CALL VLIDORT_BRDF_MAKER
     I          ( K, HAPKE_VFUNCTION, HAPKE_VFUNCTION,
     I            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS )
        ENDIF

!  Rahman kernel (3 free parameters)

        IF ( WHICH_BRDF(K) .EQ. RAHMAN_IDX ) THEN
          CALL VLIDORT_BRDF_MAKER
     I          ( K, RAHMAN_VFUNCTION, RAHMAN_VFUNCTION,
     I            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS )
        ENDIF

!  Roujean kernel (0 free parameters)

        IF ( WHICH_BRDF(K) .EQ. ROUJEAN_IDX ) THEN
          CALL VLIDORT_BRDF_MAKER
     I          ( K, ROUJEAN_VFUNCTION, ROUJEAN_VFUNCTION,
     I            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS )
        ENDIF

!  Scalar-only original Cox-Munk kernel: (2 free parameters) 

        IF ( WHICH_BRDF(K) .EQ. COXMUNK_IDX ) THEN
          IF ( DO_COXMUNK_DBMS ) THEN
            CALL VLIDORT_BRDF_MAKER
     I          ( K, COXMUNK_VFUNCTION, COXMUNK_VFUNCTION_DB, 
     I            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS )
          ELSE
            CALL VLIDORT_BRDF_MAKER
     I          ( K, COXMUNK_VFUNCTION, COXMUNK_VFUNCTION, 
     I            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS )
          ENDIF
        ENDIF

!  GISS Cox-Munk kernel: (2 free parameters) 

        IF ( WHICH_BRDF(K) .EQ. GISSCOXMUNK_IDX ) THEN
          IF ( DO_COXMUNK_DBMS ) THEN
            CALL VLIDORT_BRDF_MAKER
     I          ( K, GISSCOXMUNK_VFUNCTION, GISSCOXMUNK_VFUNCTION_DB, 
     I            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS )
          ELSE
            CALL VLIDORT_BRDF_MAKER
     I          ( K, GISSCOXMUNK_VFUNCTION, GISSCOXMUNK_VFUNCTION, 
     I            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS )
          ENDIF
        ENDIF

!  Rondeaux-Herman (vegetation) model (3 free parameters)

        IF ( WHICH_BRDF(K) .EQ. RHERMAN_IDX ) THEN
          CALL VLIDORT_BRDF_MAKER
     I          ( K, RHERMAN_VFUNCTION, RHERMAN_VFUNCTION,
     I            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS )
        ENDIF

!  Breon et al (desert) model (3 free parameters)

        IF ( WHICH_BRDF(K) .EQ. BREON_IDX ) THEN
          CALL VLIDORT_BRDF_MAKER
     I          ( K, BREON_VFUNCTION, BREON_VFUNCTION,
     I            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS )
        ENDIF

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_BRDF_CREATOR

!

      SUBROUTINE BRDF_QUADRATURE

!  Exactly the same as the scalar version

!  Include files
!  =============

!  include file of dimensions and numbers

      INCLUDE 'VLIDORT.PARS'

!  Include file of arguments controlling type of surface

      INCLUDE 'VLIDORT_INPUTS.VARS'

!  output arguments go in here.

      INCLUDE 'VLIDORT_BRDF.VARS'

!  local variables
!  ---------------

      INTEGER             I, I1, K

!  BRDF quadrature (Gauss-Legendre)
!  ---------------

!  Save these quantities for efficient coding

      NBRDF_HALF = NSTREAMS_BRDF / 2
      CALL GAULEG ( ZERO, ONE, X_BRDF, A_BRDF, NBRDF_HALF )
        DO I = 1, NBRDF_HALF
        I1 = I + NBRDF_HALF
          X_BRDF(I1) = - X_BRDF(I)
          A_BRDF(I1) =   A_BRDF(I)
        CXE_BRDF(I) = X_BRDF(I)
        SXE_BRDF(I) = DSQRT(ONE-X_BRDF(I)*X_BRDF(I))
        ENDDO
        DO I = 1, NSTREAMS_BRDF
       X_BRDF(I) = PIE * X_BRDF(I)
       CX_BRDF(I) = DCOS ( X_BRDF(I) )
       SX_BRDF(I) = DSIN ( X_BRDF(I) )
      ENDDO

!  Half space cosine-weight arrays (emission only, non-Lambertian)

      IF ( DO_SURFACE_EMISSION .AND.
     &       .NOT. DO_LAMBERTIAN_SURFACE) THEN
        DO K = 1, NBRDF_HALF
          BAX_BRDF(K) = X_BRDF(K) * A_BRDF(K) / PIE
        ENDDO
      ENDIF

!  Finish

      END SUBROUTINE BRDF_QUADRATURE

!

      SUBROUTINE BRDF_QUADRATURE_TPZ

!  Include files
!  =============

!  include file of dimensions and numbers

      INCLUDE 'VLIDORT.PARS'

!  Include file of arguments controlling type of surface

      INCLUDE 'VLIDORT_INPUTS.VARS'

!  output arguments go in here.

      INCLUDE 'VLIDORT_BRDF.VARS'

!  local variables
!  ---------------

      INTEGER             I, I1
      DOUBLE PRECISION DF1, DEL

!  BRDF quadrature (Trapezium)
!  ---------------

!  Save these quantities for efficient coding

      DF1 = DFLOAT(NSTREAMS_BRDF - 1 )
      DEL = TWO * PIE / DF1
        DO I = 1, NSTREAMS_BRDF
        I1 = I - 1
        X_BRDF(I) = DFLOAT(I1) * DEL - PIE
        X_BRDF(I) = DFLOAT(I1) * DEL
        CX_BRDF(I) = DCOS ( X_BRDF(I) )
        SX_BRDF(I) = DSIN ( X_BRDF(I) )
        CXE_BRDF(I) = CX_BRDF(I)
        SXE_BRDF(I) = SX_BRDF(I)
      ENDDO
        DO I = 2, NSTREAMS_BRDF - 1
        A_BRDF(I)  = DEL / PIE
      ENDDO
      A_BRDF(1)              = DEL * HALF / PIE
      A_BRDF(NSTREAMS_BRDF)  = DEL * HALF / PIE

!  Half space cosine-weight arrays (emission only, non-Lambertian)

!      IF ( DO_SURFACE_EMISSION .AND.
!     &       .NOT. DO_LAMBERTIAN_SURFACE) THEN
!        DO K = 1, NBRDF_HALF
!          BAX_BRDF(K) = X_BRDF(K) * A_BRDF(K) / PIE
!        ENDDO
!      ENDIF

!  Finish

      END SUBROUTINE BRDF_QUADRATURE_TPZ

!

      SUBROUTINE VLIDORT_BRDF_MAKER
     I       ( BRDF_INDEX, BRDF_VFUNCTION, BRDF_VFUNCTION_DB,
     I         LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS )

!  Prepares the bidirectional reflectance scatter matrices

!  Include files
!  -------------

!  include file of dimensions and numbers

      INCLUDE 'VLIDORT.PARS'

!  Include files of input and bookkeeping variables 

      INCLUDE 'VLIDORT_INPUTS.VARS'
      INCLUDE 'VLIDORT_BOOKKEEP.VARS'

!  output arguments go in here.

      INCLUDE 'VLIDORT_BRDF.VARS'

!  Module arguments
!  ----------------

!  Index

      INTEGER        BRDF_INDEX

!  BRDF functions (external calls)

      EXTERNAL       BRDF_VFUNCTION
      EXTERNAL       BRDF_VFUNCTION_DB

!  Local number of parameters and local parameter array

      INTEGER          LOCAL_BRDF_NPARS
      DOUBLE PRECISION LOCAL_BRDF_PARS ( MAX_BRDF_PARAMETERS )
      
!  local variables
!  ---------------

      INTEGER          I, UI, J, K, KE, IB
      DOUBLE PRECISION MUX, SZASURCOS(MAXBEAMS),SZASURSIN(MAXBEAMS)
      DOUBLE PRECISION PHIANG, COSPHI, SINPHI
      DOUBLE PRECISION SZAANG, COSSZA, SINSZA
      DOUBLE PRECISION VZAANG, COSVZA, SINVZA

!  Exact DB calculation
!  --------------------

!  New code. 06 August 2007. RT Solutions Inc.
!   Must use adjusted values of BOA geometry
!   ( Adjusted values = original, if GEOMETRY_SPECHEIGHT = BOA height )

      IF ( DO_DBCORRECTION ) THEN
        DO UI = 1, N_USER_STREAMS
          VZAANG = USER_VZANGLES_ADJUST(UI)
          COSVZA = DCOS(VZAANG*DEG_TO_RAD)
          SINVZA = DSIN(VZAANG*DEG_TO_RAD)        
          DO IB = 1, NBEAMS
            DO K = 1, N_USER_RELAZMS
              SZAANG = SZANGLES_ADJUST(UI,IB,K)
              COSSZA = DCOS(SZAANG*DEG_TO_RAD)
              SINSZA = DSIN(SZAANG*DEG_TO_RAD)
              PHIANG = USER_RELAZMS_ADJUST(UI,IB,K)
              COSPHI = DCOS(PHIANG*DEG_TO_RAD)
              SINPHI = DSIN(PHIANG*DEG_TO_RAD)
              CALL BRDF_VFUNCTION_DB
     C        ( MAX_BRDF_PARAMETERS, N_BRDF_STOKESSQ,
     C          LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,
     G          COSSZA, SINSZA, COSVZA, SINVZA,
     G          PHIANG, COSPHI, SINPHI,
     O          EXACTDB_BRDFUNC(1,BRDF_INDEX,UI,K,IB) )
            ENDDO
          ENDDO
        ENDDO
      ENDIF

      IF ( DO_SSFULL ) RETURN

!  Regular RT calculation
!  ======================

!  Usable solar beams
!  ------------------

      IF ( DO_REFRACTIVE_GEOMETRY ) THEN
        DO IB = 1, NBEAMS
          MUX =  DCOS(SZA_LOCAL_INPUT(NLAYERS,IB)*DEG_TO_RAD)
          SZASURCOS(IB) = MUX
          SZASURSIN(IB) = DSQRT(ONE-MUX*MUX)
        ENDDO
      ELSE
        DO IB = 1, NBEAMS
          SZASURCOS(IB) = COS_SZANGLES(IB)
          SZASURSIN(IB) = SIN_SZANGLES(IB)
        ENDDO
      ENDIF

!  Quadrature outgoing directions
!  ------------------------------

!  Incident Solar beam

      DO IB = 1, NBEAMS 
       DO I = 1, NSTREAMS
        DO K = 1, NSTREAMS_BRDF
         CALL BRDF_VFUNCTION
     C        ( MAX_BRDF_PARAMETERS, N_BRDF_STOKESSQ,
     C          LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,
     G          SZASURCOS(IB), SZASURSIN(IB),
     G          QUAD_STREAMS(I), QUAD_SINES(I),
     G          X_BRDF(K), CX_BRDF(K), SX_BRDF(K),
     O          BRDFUNC_0(1,BRDF_INDEX,I,IB,K) )
        ENDDO
       ENDDO
      ENDDO

!  incident quadrature directions

      DO I = 1, NSTREAMS
        DO J = 1, NSTREAMS
          DO K = 1, NSTREAMS_BRDF
            CALL BRDF_VFUNCTION
     C          ( MAX_BRDF_PARAMETERS, N_BRDF_STOKESSQ,
     C            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,
     G            QUAD_STREAMS(J), QUAD_SINES(J),
     G            QUAD_STREAMS(I), QUAD_SINES(I),
     G            X_BRDF(K), CX_BRDF(K), SX_BRDF(K),
     O            BRDFUNC(1,BRDF_INDEX,I,J,K) )
          ENDDO
        ENDDO
      ENDDO

!  Emissivity (optional) - BRDF quadrature input directions

      IF ( DO_SURFACE_EMISSION ) THEN
        NBRDF_HALF = NSTREAMS_BRDF / 2
        DO I = 1, NSTREAMS
          DO KE = 1, NBRDF_HALF
            DO K = 1, NSTREAMS_BRDF
              CALL BRDF_VFUNCTION
     C            ( MAX_BRDF_PARAMETERS, N_EMISS_STOKESSQ,
     C              LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,
     G              CXE_BRDF(KE), SXE_BRDF(KE),
     G              QUAD_STREAMS(I), QUAD_SINES(I),
     G              X_BRDF(K), CX_BRDF(K), SX_BRDF(K),
     O              EBRDFUNC(1,BRDF_INDEX,I,KE,K) )
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!  User-streams outgoing directions
!  --------------------------------

      IF ( DO_USER_STREAMS ) THEN

!  Incident Solar beam

        DO IB = 1, NBEAMS
         DO UI = 1, N_USER_STREAMS
          DO K = 1, NSTREAMS_BRDF
            CALL BRDF_VFUNCTION
     C          ( MAX_BRDF_PARAMETERS, N_BRDF_STOKESSQ,
     C            LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,
     G            SZASURCOS(IB), SZASURSIN(IB),
     G            USER_STREAMS(UI), USER_SINES(UI),
     G            X_BRDF(K), CX_BRDF(K), SX_BRDF(K),
     O            USER_BRDFUNC_0(1,BRDF_INDEX,UI,IB,K) )
          ENDDO
         ENDDO
        ENDDO

!  incident quadrature directions

        DO UI = 1, N_USER_STREAMS
          DO J = 1, NSTREAMS
            DO K = 1, NSTREAMS_BRDF
              CALL BRDF_VFUNCTION
     C            ( MAX_BRDF_PARAMETERS, N_BRDF_STOKESSQ,
     C              LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,
     G              QUAD_STREAMS(J), QUAD_SINES(J),
     I              USER_STREAMS(UI), USER_SINES(UI),
     G              X_BRDF(K), CX_BRDF(K), SX_BRDF(K),
     O              USER_BRDFUNC(1,BRDF_INDEX,UI,J,K) )
            ENDDO
          ENDDO
        ENDDO

!  Emissivity (optional) - BRDF quadrature input directions

        IF ( DO_SURFACE_EMISSION ) THEN
          DO UI = 1, N_USER_STREAMS
            DO KE = 1, NBRDF_HALF
              DO K = 1, NSTREAMS_BRDF
              CALL BRDF_VFUNCTION
     C            ( MAX_BRDF_PARAMETERS, N_EMISS_STOKESSQ,
     C              LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,
     I              CXE_BRDF(KE), SXE_BRDF(KE),
     I              USER_STREAMS(UI), USER_SINES(UI), 
     G              X_BRDF(K), CX_BRDF(K), SX_BRDF(K),
     O              USER_EBRDFUNC(1,BRDF_INDEX,UI,KE,K) )
              ENDDO
            ENDDO
          ENDDO
        ENDIF

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_BRDF_MAKER

!

      SUBROUTINE VLIDORT_BRDF_FRTERMS
     I    ( DO_INCLUDE_SURFACE,
     I      DO_INCLUDE_SURFEMISS,
     I      FOURIER, DELFAC )

!  Prepares Fourier components of the bidirectional reflectance functions
!    necessary for VLIDORT.

!  Include files
!  -------------

!  include file of dimensions and numbers

      INCLUDE 'VLIDORT.PARS'

!  Include files of input and bookkeeping variables 

      INCLUDE 'VLIDORT_INPUTS.VARS'
      INCLUDE 'VLIDORT_BOOKKEEP.VARS'

!  include file of setup variables

      INCLUDE 'VLIDORT_SETUPS.VARS'

!  Input and output arguments go in here.

      INCLUDE 'VLIDORT_REFLECTANCE.VARS'
      INCLUDE 'VLIDORT_BRDF.VARS'

!  Module arguments
!  ----------------

      LOGICAL             DO_INCLUDE_SURFACE
      LOGICAL             DO_INCLUDE_SURFEMISS
      INTEGER             FOURIER
      DOUBLE PRECISION    DELFAC

!  local variables

      INTEGER             I, UI, J, K, KPHI, KL, Q, IB, O1
      DOUBLE PRECISION    SUM, REFL, HELP, HELP_A
      INTEGER             COSSIN_MASK(16)
      DATA COSSIN_MASK / 1,1,2,0,1,1,2,0,2,2,1,0,0,0,0,1 /

!  Weighted azimuth factor
!  ( Results are stored in commons )

      IF ( FOURIER .NE. 0 ) THEN
        DO K = 1, NSTREAMS_BRDF
          BRDF_COSAZMFAC(K) = A_BRDF(K) * DCOS ( FOURIER * X_BRDF(K) )
          BRDF_SINAZMFAC(K) = A_BRDF(K) * DSIN ( FOURIER * X_BRDF(K) )
        ENDDO
      ELSE
        DO K = 1, NSTREAMS_BRDF
          BRDF_COSAZMFAC(K) = A_BRDF(K)
          BRDF_SINAZMFAC(K) = ZERO
        ENDDO
      ENDIF

!  factor

      HELP = HALF * DELFAC

!  Quadrature outgoing directions
!  ------------------------------

!  Incident Solar beam (direct beam reflections)

      DO KL = 1, N_BRDF_KERNELS
       IF ( .NOT. LAMBERTIAN_KERNEL_FLAG(KL) ) THEN
        DO IB = 1, NBEAMS
         IF ( DO_REFLECTED_DIRECTBEAM(IB) ) THEN
          DO I = 1, NSTREAMS
           DO Q = 1, N_BRDF_STOKESSQ
            SUM = ZERO
            IF ( COSSIN_MASK(Q).EQ.1 ) THEN
             DO K = 1, NSTREAMS_BRDF
              SUM  = SUM + BRDFUNC_0(Q,KL,I,IB,K)*BRDF_COSAZMFAC(K)
             ENDDO
            ELSE IF ( COSSIN_MASK(Q).EQ.2 ) THEN
             DO K = 1, NSTREAMS_BRDF
              SUM  = SUM + BRDFUNC_0(Q,KL,I,IB,K)*BRDF_SINAZMFAC(K)
             ENDDO
            ENDIF
            BIREFLEC_0(Q,KL,I,IB) = SUM * HELP
           ENDDO
          ENDDO
         ENDIF
        ENDDO
       ENDIF
      ENDDO

!  incident quadrature directions (surface multiple reflections)

      IF ( DO_INCLUDE_SURFACE ) THEN
        DO KL = 1, N_BRDF_KERNELS
          IF ( .NOT. LAMBERTIAN_KERNEL_FLAG(KL) ) THEN
            DO I = 1, NSTREAMS
              DO J = 1, NSTREAMS
               DO Q = 1, N_BRDF_STOKESSQ
                SUM = ZERO
                IF ( COSSIN_MASK(Q).EQ.1 ) THEN
                 DO K = 1, NSTREAMS_BRDF
                  SUM  = SUM + BRDFUNC(Q,KL,I,J,K)*BRDF_COSAZMFAC(K)
                 ENDDO
                ELSE IF ( COSSIN_MASK(Q).EQ.2 ) THEN
                 DO K = 1, NSTREAMS_BRDF
                  SUM  = SUM + BRDFUNC(Q,KL,I,J,K)*BRDF_SINAZMFAC(K)
                 ENDDO
                ENDIF
                BIREFLEC(Q,KL,I,J) = SUM * HELP
               ENDDO
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ENDIF

!  debug information

      IF ( DO_DEBUG_WRITE ) THEN
        WRITE(55,'(A)')'BRDF_1 Fourier 0 quad values'
        IF ( FOURIER .EQ. 0 ) THEN
          DO I = 1, NSTREAMS
          WRITE(55,'(1PE12.5,3x,1P10E12.5)')
     &     BIREFLEC_0(1,1,I,1),(BIREFLEC(1,1,I,J),J=1,NSTREAMS)
         ENDDO
        ENDIF
      ENDIF
  
!  albedo check, always calculate the spherical albedo.
!   (Plane albedo calculations are commented out)

      IF ( FOURIER .EQ. 0 ) THEN

        DO KL = 1, N_BRDF_KERNELS

!  ............................... Plane albedo calculations
!          SUM = ZERO
!          DO I = 1, NSTREAMS
!            SUM = SUM + BIREFLEC_0(KL,I) * AX(I)
!          ENDDO
!          SUM = SUM * TWO
!          write(*,*)0,BRDF_FACTORS(KL),sum
!          DO J = 1, NSTREAMS
!            SUM = ZERO
!            DO I = 1, NSTREAMS
!              SUM = SUM + BIREFLEC(KL,I,J) * AX(I)
!            ENDDO
!            SUM = SUM * TWO
!            write(*,*)j,BRDF_FACTORS(KL),sum
!          ENDDO

!  ............................... Spherical albedo calculation

          HELP_A = ZERO
          DO I = 1, NSTREAMS
            SUM = ZERO
            DO J = 1, NSTREAMS
               SUM = SUM + BIREFLEC(1,KL,I,J) * QUAD_STRMWTS(J)
            ENDDO
            HELP_A = HELP_A + SUM * QUAD_STRMWTS(I)
          ENDDO
          SPHERICAL_ALBEDO(KL) = HELP_A*FOUR

        ENDDO
      ENDIF

!  User-streams outgoing directions
!  --------------------------------

      IF ( DO_USER_STREAMS ) THEN

!  Incident Solar beam (direct beam reflections)

       DO KL = 1, N_BRDF_KERNELS
        IF ( .NOT. LAMBERTIAN_KERNEL_FLAG(KL) ) THEN
         DO IB = 1, NBEAMS
          IF ( DO_REFLECTED_DIRECTBEAM(IB) ) THEN
           DO UI = 1, N_USER_STREAMS
            DO Q = 1, N_BRDF_STOKESSQ
             SUM = ZERO
             IF ( COSSIN_MASK(Q).EQ.1 ) THEN
              DO K = 1, NSTREAMS_BRDF
               SUM = SUM+USER_BRDFUNC_0(Q,KL,UI,IB,K)*BRDF_COSAZMFAC(K)
              ENDDO
             ELSE IF ( COSSIN_MASK(Q).EQ.2 ) THEN
              DO K = 1, NSTREAMS_BRDF
               SUM = SUM+USER_BRDFUNC_0(Q,KL,UI,IB,K)*BRDF_SINAZMFAC(K)
              ENDDO
             ENDIF
             USER_BIREFLEC_0(Q,KL,UI,IB) = SUM * HELP
            ENDDO
           ENDDO
          ENDIF
         ENDDO
        ENDIF
       ENDDO

!  incident quadrature directions (surface multiple reflections)

       IF ( DO_INCLUDE_SURFACE ) THEN
        DO KL = 1, N_BRDF_KERNELS
         IF ( .NOT. LAMBERTIAN_KERNEL_FLAG(KL) ) THEN
          DO UI = 1, N_USER_STREAMS
           DO J = 1, NSTREAMS
            DO Q = 1, N_BRDF_STOKESSQ
             SUM = ZERO
             IF ( COSSIN_MASK(Q).EQ.1 ) THEN
              DO K = 1, NSTREAMS_BRDF
               SUM = SUM+USER_BRDFUNC(Q,KL,UI,J,K)*BRDF_COSAZMFAC(K)
              ENDDO
             ELSE IF ( COSSIN_MASK(Q).EQ.2 ) THEN
              DO K = 1, NSTREAMS_BRDF
               SUM = SUM+USER_BRDFUNC(Q,KL,UI,J,K)*BRDF_SINAZMFAC(K)
              ENDDO
             ENDIF
             USER_BIREFLEC(Q,KL,UI,J) = SUM * HELP
            ENDDO
           ENDDO
          ENDDO
         ENDIF
        ENDDO
       ENDIF

      ENDIF

!  Emissivity
!  ----------

!  Assumed to exist only for the total intensity
!    (first element of Stokes Vector) - is this right ?????? ---> YES

      IF ( DO_INCLUDE_SURFEMISS ) THEN

!  Initialise

        DO I = 1, NSTREAMS
          EMISSIVITY(I,1) = ONE
          DO O1 = 2, NSTOKES
            EMISSIVITY(I,O1) = ZERO
          ENDDO
        ENDDO
        IF ( DO_USER_STREAMS ) THEN
          DO UI = 1, N_USER_STREAMS
            USER_EMISSIVITY(UI,1) = ONE
            DO O1 = 2, NSTOKES
              USER_EMISSIVITY(UI,O1) = ZERO
            ENDDO
          ENDDO
        ENDIF

!  Start loop over kernels

        DO KL = 1, N_BRDF_KERNELS

!  Lambertian case

          IF ( LAMBERTIAN_KERNEL_FLAG(KL) ) THEN

            DO I = 1, NSTREAMS
              A_EMISSIVITY(I,KL) = BRDF_FACTORS(KL)
            ENDDO
            IF ( DO_USER_STREAMS ) THEN
              DO UI = 1, N_USER_STREAMS
                A_USER_EMISSIVITY(UI,KL) = BRDF_FACTORS(KL)
              ENDDO
            ENDIF

!  bidirectional reflectance

          ELSE

!  Quadrature polar directions

            DO I = 1, NSTREAMS
              REFL = ZERO
              DO KPHI= 1, NSTREAMS_BRDF
                SUM = ZERO
                DO K = 1, NBRDF_HALF
                  SUM = SUM + EBRDFUNC(1,KL,I,K,KPHI) * BAX_BRDF(K)
                ENDDO
                REFL = REFL + A_BRDF(KPHI) * SUM
              ENDDO
              A_EMISSIVITY(I,KL) = REFL * BRDF_FACTORS(KL)
            ENDDO

!   user-defined polar directions

            IF ( DO_USER_STREAMS ) THEN
              DO UI = 1, N_USER_STREAMS
                REFL = ZERO
                DO KPHI= 1, NSTREAMS_BRDF
                  SUM = ZERO
                  DO K = 1, NBRDF_HALF
                    SUM = SUM+USER_EBRDFUNC(1,KL,UI,K,KPHI)*BAX_BRDF(K)
                  ENDDO
                  REFL = REFL + A_BRDF(KPHI) * SUM
                ENDDO
                A_USER_EMISSIVITY(UI,KL) = REFL*BRDF_FACTORS(KL)
              ENDDO
            ENDIF

          ENDIF

!  end loop over kernels

        ENDDO

!  Total emissivities

        DO KL = 1, N_BRDF_KERNELS
          DO I = 1, NSTREAMS
            EMISSIVITY(I,1) = EMISSIVITY(I,1) - A_EMISSIVITY(I,KL)
          ENDDO
          IF ( DO_USER_STREAMS ) THEN
            DO UI = 1, N_USER_STREAMS
              USER_EMISSIVITY(UI,1) = USER_EMISSIVITY(UI,1) -
     &               A_USER_EMISSIVITY(UI,KL)
            ENDDO
          ENDIF
        ENDDO

!  end emissivity clause

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_BRDF_FRTERMS

!

      SUBROUTINE BRDF_SURFACE_DIRECTBEAM
     I     ( FOURIER_COMPONENT,
     I       DELTA_FACTOR )

!  Computation of direct beam contributions for a BRDF surface

!  include files
!  -------------

!  include file of dimensions and numbers

      INCLUDE 'VLIDORT.PARS'

!  Include files of input and bookkeeping variables

      INCLUDE 'VLIDORT_INPUTS.VARS'
      INCLUDE 'VLIDORT_BOOKKEEP.VARS'

!  Include file of setup variables (Input to the present module)

      INCLUDE 'VLIDORT_SETUPS.VARS'

!  Include file of BRDF variables (input and some output)

      INCLUDE 'VLIDORT_BRDF.VARS'

!  Include file of Direct beam reflectance variables (output goes here)

      INCLUDE 'VLIDORT_REFLECTANCE.VARS'

!  module arguments
!  ----------------

      DOUBLE PRECISION DELTA_FACTOR
      INTEGER          FOURIER_COMPONENT

!  Local variables
!  ---------------

      DOUBLE PRECISION X0_FLUX, X0_BOA, ATTN, REFL_ATTN, SUM
      INTEGER          I, UI, KL, O1, O2, OM, IB

!  Initialize
!  ----------

!   Safety first!  Return if there is no reflection.

      DO IB = 1, NBEAMS
        DO I = 1, NSTREAMS
          DO O1 = 1, NSTOKES
            DIRECT_BEAM(I,IB,O1) = ZERO
          ENDDO
        ENDDO
        IF ( DO_USER_STREAMS ) THEN
          DO UI = 1, N_USER_STREAMS
            DO O1 = 1, NSTOKES
              USER_DIRECT_BEAM(UI,IB,O1) = ZERO
            ENDDO
          ENDDO
        ENDIF
      ENDDO

!  Attenuation of solar beam
!  -------------------------

!  New code to deal with refractive geometry case
!   R. Spurr, 7 May 2005. RT Solutions Inc.

      DO IB = 1, NBEAMS
       IF ( DO_REFLECTED_DIRECTBEAM(IB) ) THEN

        IF ( DO_REFRACTIVE_GEOMETRY ) THEN
         X0_BOA = DCOS(SZA_LOCAL_INPUT(NLAYERS,IB)*DEG_TO_RAD)
        ELSE
         X0_BOA = COS_SZANGLES(IB)
        ENDIF

!  There should be no flux factor here.
!    Bug fixed 18 November 2005. Earlier Italian Job!!
!       X0_FLUX        = FOUR * X0_BOA * FLUX_FACTOR / DELTA_FACTOR

        X0_FLUX        = FOUR * X0_BOA / DELTA_FACTOR
        ATTN           = X0_FLUX * SOLAR_BEAM_OPDEP(IB)
        ATMOS_ATTN(IB) = ATTN

!  Start loop over albedo kernels
!    - reflected along the quadrature directions
!    - reflected along the user-defined directions

        DO KL = 1, N_BRDF_KERNELS

          REFL_ATTN = BRDF_FACTORS(KL) * ATTN

!  For the Lambertian case
!  -----------------------

!  All terms are the same; zero for Fourier not zero

          IF ( LAMBERTIAN_KERNEL_FLAG(KL) ) THEN

!  zero it!

            DO O1 = 1, NSTOKES
              DO I = 1, NSTREAMS
                A_DIRECT_BEAM(I,IB,O1,KL) = ZERO
              ENDDO
            ENDDO
            IF ( DO_USER_STREAMS ) THEN
              DO O1 = 1, NSTOKES
                DO UI = 1, N_USER_STREAMS
                  A_USER_DIRECT_SOURCE(UI,IB,O1,KL) = ZERO
                ENDDO
              ENDDO
            ENDIF

!  Set value for Fourier = 0
!  Natural light only, so only the first Stokes component is nonzero

            IF ( FOURIER_COMPONENT .EQ. 0 ) THEN
              DO I = 1, NSTREAMS
                A_DIRECT_BEAM(I,IB,1,KL) = REFL_ATTN
              ENDDO
              IF ( DO_USER_STREAMS ) THEN
                DO UI = 1, N_USER_STREAMS
                  A_USER_DIRECT_SOURCE(UI,IB,1,KL) = REFL_ATTN
                ENDDO
              ENDIF
            ENDIF

!  For the Bidirectional case
!  --------------------------

          ELSE

            DO I = 1, NSTREAMS
              DO O1 = 1, NSTOKES
                SUM = ZERO
                DO O2 = 1, NSTOKES
                  OM = MUELLER_INDEX(O1,O2)
                  SUM = SUM + FLUXVEC(O2) * BIREFLEC_0(OM,KL,I,IB)
                ENDDO
                A_DIRECT_BEAM(I,IB,O1,KL) = REFL_ATTN * SUM
              ENDDO
            ENDDO
            IF ( DO_USER_STREAMS ) THEN
              DO UI = 1, N_USER_STREAMS
                DO O1 = 1, NSTOKES
                  SUM = ZERO
                  DO O2 = 1, NSTOKES
                    OM = MUELLER_INDEX(O1,O2)
                    SUM = SUM + FLUXVEC(O2)*USER_BIREFLEC_0(OM,KL,UI,IB)
                  ENDDO
                  A_USER_DIRECT_SOURCE(UI,IB,O1,KL) = REFL_ATTN * SUM
                ENDDO
              ENDDO
            ENDIF

          ENDIF

!  Finish Kernel loop

        ENDDO

!  Total contributions (summed over the kernels)
!  -------------------

        DO KL = 1, N_BRDF_KERNELS
          DO I = 1, NSTREAMS
            DO O1 = 1, NSTOKES
              DIRECT_BEAM(I,IB,O1) = DIRECT_BEAM(I,IB,O1) +
     &           A_DIRECT_BEAM(I,IB,O1,KL)
            ENDDO
          ENDDO
          IF ( DO_USER_STREAMS ) THEN
            DO UI = 1, N_USER_STREAMS
              DO O1 = 1, NSTOKES
                USER_DIRECT_BEAM(UI,IB,O1) = 
     &                  USER_DIRECT_BEAM(UI,IB,O1) + 
     &                 A_USER_DIRECT_SOURCE(UI,IB,O1,KL)
              ENDDO
            ENDDO
          ENDIF
        ENDDO

!  end direct beam calculation

       ENDIF
      ENDDO

!  finish

      RETURN
      END SUBROUTINE BRDF_SURFACE_DIRECTBEAM

!

      SUBROUTINE BVP_BRDF_SURFACE_HOM
     I     ( DO_INCLUDE_SURFACE,
     I       FOURIER_COMPONENT,
     I       SURFACE_FACTOR )

!  Bidirectional surface reflectance

!  Include files
!  -------------

!  Additional sums for the final albedo-reflecting layer

!  include file of dimensiopns and numbers

      INCLUDE 'VLIDORT.PARS'

!  include files of input and bookkeeping variables

      INCLUDE 'VLIDORT_INPUTS.VARS'
      INCLUDE 'VLIDORT_BOOKKEEP.VARS'

!  include file of solution stuff (input to this module)

      INCLUDE 'VLIDORT_SOLUTION.VARS'

!  Output include files (surface reflectance, kernel contributions)

      INCLUDE 'VLIDORT_REFLECTANCE.VARS'
      INCLUDE 'VLIDORT_BRDF.VARS'

!  control input
!  -------------

      LOGICAL          DO_INCLUDE_SURFACE
      INTEGER          FOURIER_COMPONENT
      DOUBLE PRECISION SURFACE_FACTOR

!  local variables
!  ---------------

      INTEGER          I,J,O1,O2,OM,KL,K,KO1,K0,K1,K2,NELEMENTS,M,NL
      DOUBLE PRECISION FACTOR_KERNEL
      DOUBLE PRECISION H_1,   H_2,   H_1_S,   H_2_S
      DOUBLE PRECISION H_1_CR,   H_2_CR,   H_1_CI,   H_2_CI
      DOUBLE PRECISION H_1_S_CR, H_2_S_CR, H_1_S_CI, H_2_S_CI

!  Initialization
!  ==============

!  Zero OM

      OM = 0
      M = FOURIER_COMPONENT

!  Zero total reflected contributions

      DO I = 1, NSTREAMS
        DO K = 1, NSTKS_NSTRMS
          DO O1 = 1, NSTOKES
            R2_HOMP(I,O1,K) = ZERO
            R2_HOMM(I,O1,K) = ZERO
          ENDDO
        ENDDO
      ENDDO

!  Return with Zeroed values if albedo flag not set

      IF ( .NOT. DO_INCLUDE_SURFACE ) RETURN

!  Last layer

      NL = NLAYERS

!  Offset

      KO1 = K_REAL(NL) + 1

!  Integrate Downward streams of particular solutions

      DO KL = 1, N_BRDF_KERNELS

        FACTOR_KERNEL = SURFACE_FACTOR * BRDF_FACTORS(KL)
        
!  For Lambertian reflectance, all streams are the same
!  ----------------------------------------------------

        IF ( LAMBERTIAN_KERNEL_FLAG(KL) ) THEN

          IF ( FOURIER_COMPONENT .EQ. 0 ) THEN

!  Zero the elements away from NSTOKES = 1

            DO I = 1, NSTREAMS
              DO O1 = 2, NSTOKES
                DO K = 1, NSTKS_NSTRMS
                  A_DIFFUSE_HOMP(I,O1,K,KL) = ZERO
                  A_DIFFUSE_HOMM(I,O1,K,KL) = ZERO
                ENDDO
              ENDDO
            ENDDO

!  Homogeneous real solutions

            DO K = 1, K_REAL(NL)
              DO I = 1, NSTREAMS
                H_1 = ZERO
                H_2 = ZERO
                DO J = 1, NSTREAMS
                  H_1 = H_1 + QUAD_STRMWTS(J)*SOLA_XPOS(J,1,K,NL)
                  H_2 = H_2 + QUAD_STRMWTS(J)*SOLB_XNEG(J,1,K,NL)
                ENDDO
                A_DIFFUSE_HOMP(I,1,K,KL) = H_1 * FACTOR_KERNEL
                A_DIFFUSE_HOMM(I,1,K,KL) = H_2 * FACTOR_KERNEL
              ENDDO
            ENDDO

!  Homogeneous complex solutions

            DO K = 1, K_COMPLEX(NL)
              K0 = 2*K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              DO I = 1, NSTREAMS
                H_1_CR = ZERO
                H_2_CR = ZERO
                H_1_CI = ZERO
                H_2_CI = ZERO
                DO J = 1, NSTREAMS
                  H_1_CR = H_1_CR + QUAD_STRMWTS(J)*SOLA_XPOS(J,1,K1,NL)
                  H_2_CR = H_2_CR + QUAD_STRMWTS(J)*SOLB_XNEG(J,1,K1,NL)
                  H_1_CI = H_1_CI + QUAD_STRMWTS(J)*SOLA_XPOS(J,1,K2,NL)
                  H_2_CI = H_2_CI + QUAD_STRMWTS(J)*SOLB_XNEG(J,1,K2,NL)
                ENDDO
                A_DIFFUSE_HOMP(I,1,K1,KL) = H_1_CR * FACTOR_KERNEL
                A_DIFFUSE_HOMM(I,1,K1,KL) = H_2_CR * FACTOR_KERNEL
                A_DIFFUSE_HOMP(I,1,K2,KL) = H_1_CI * FACTOR_KERNEL
                A_DIFFUSE_HOMM(I,1,K2,KL) = H_2_CI * FACTOR_KERNEL
              ENDDO
            ENDDO

!  zero these terms otherwise

          ELSE
            DO I = 1, NSTREAMS
              DO O1 = 1, NSTOKES
                DO K = 1, NSTKS_NSTRMS
                  A_DIFFUSE_HOMP(I,O1,K,KL) = ZERO
                  A_DIFFUSE_HOMM(I,O1,K,KL) = ZERO
                ENDDO
              ENDDO
            ENDDO
          ENDIF

!  For bidirectional reflecting surface
!  ------------------------------------

        ELSE

!  help variable

          DO I = 1, NSTREAMS
           DO O1 = 1, NSTOKES
            DO J = 1, NSTREAMS
             DO O2 = 1, NSTOKES
              OM = MUELLER_INDEX(O1,O2)
              AXBID(I,J,OM,KL) = QUAD_STRMWTS(J) * BIREFLEC(OM,KL,J,I)
             ENDDO
            ENDDO
           ENDDO
          ENDDO

!  debug
!          DO O1 = 1, NSTOKES
!            DO O2 = 1, NSTOKES
!              J = MUELLER_INDEX(O1,O2)
!              I  = MUELLER_INDEX(O2,O1)
!            if (m.le.2) write(71,'(3i4,1p2e15.7)')
!     &               m,o1,o2,AXBID(4,5,J,1),AXBID(4,5,I,1)
!            ENDDO
!          ENDDO
!c          if (m.eq.2)pause'axbid'

!  homogeneous real solutions

          DO K = 1, K_REAL(NL)
            DO I = 1, NSTREAMS
              DO O1 = 1, NSTOKES
                H_1 = ZERO
                H_2 = ZERO
                DO J = 1, NSTREAMS
                  H_1_S = ZERO
                  H_2_S = ZERO
                  DO O2 = 1, NSTOKES
                    OM = MUELLER_INDEX(O1,O2)
                    H_1_S = H_1_S + 
     &                  AXBID(I,J,OM,KL) * SOLA_XPOS(J,O2,K,NL)
                    H_2_S = H_2_S +
     &                  AXBID(I,J,OM,KL) * SOLB_XNEG(J,O2,K,NL)
                  ENDDO
                  H_1 = H_1 + H_1_S
                  H_2 = H_2 + H_2_S
                ENDDO
                A_DIFFUSE_HOMP(I,O1,K,KL) = FACTOR_KERNEL * H_1
                A_DIFFUSE_HOMM(I,O1,K,KL) = FACTOR_KERNEL * H_2
              ENDDO
            ENDDO
          ENDDO

!  homogeneous complex solutions

          DO K = 1, K_COMPLEX(NL)
            K0 = 2*K - 2
            K1 = KO1 + K0
            K2 = K1  + 1
            DO I = 1, NSTREAMS
              DO O1 = 1, NSTOKES
                H_1_CR = ZERO
                H_2_CR = ZERO
                H_1_CI = ZERO
                H_2_CI = ZERO
                DO J = 1, NSTREAMS
                  H_1_S_CR = ZERO
                  H_2_S_CR = ZERO
                  H_1_S_CI = ZERO
                  H_2_S_CI = ZERO
                  DO O2 = 1, NSTOKES
                    OM = MUELLER_INDEX(O1,O2)
                    H_1_S_CR = H_1_S_CR +
     &                AXBID(I,J,OM,KL) * SOLA_XPOS(J,O2,K1,NL)
                    H_2_S_CR = H_2_S_CR +
     &                AXBID(I,J,OM,KL) * SOLB_XNEG(J,O2,K1,NL)
                    H_1_S_CI = H_1_S_CI +
     &                AXBID(I,J,OM,KL) * SOLA_XPOS(J,O2,K2,NL)
                    H_2_S_CI = H_2_S_CI +
     &                AXBID(I,J,OM,KL) * SOLB_XNEG(J,O2,K2,NL)
                  ENDDO
                  H_1_CR = H_1_CR + H_1_S_CR
                  H_2_CR = H_2_CR + H_2_S_CR
                  H_1_CI = H_1_CI + H_1_S_CI
                  H_2_CI = H_2_CI + H_2_S_CI
                ENDDO
                A_DIFFUSE_HOMP(I,O1,K1,KL) = FACTOR_KERNEL * H_1_CR
                A_DIFFUSE_HOMM(I,O1,K1,KL) = FACTOR_KERNEL * H_2_CR
                A_DIFFUSE_HOMP(I,O1,K2,KL) = FACTOR_KERNEL * H_1_CI
                A_DIFFUSE_HOMM(I,O1,K2,KL) = FACTOR_KERNEL * H_2_CI
              ENDDO
            ENDDO
          ENDDO

        ENDIF

!  End loop over albedo kernels

      ENDDO

!  Total terms
!  -----------

      DO KL = 1, N_BRDF_KERNELS
       NELEMENTS = NSTOKES
       IF ( LAMBERTIAN_KERNEL_FLAG(KL) ) NELEMENTS = 1
       DO I = 1, NSTREAMS
        DO O1 = 1, NELEMENTS
         DO K = 1, NSTKS_NSTRMS
          R2_HOMP(I,O1,K) = R2_HOMP(I,O1,K) + A_DIFFUSE_HOMP(I,O1,K,KL)
          R2_HOMM(I,O1,K) = R2_HOMM(I,O1,K) + A_DIFFUSE_HOMM(I,O1,K,KL)
         ENDDO
        ENDDO
       ENDDO
      ENDDO

!  debug

!      DO I = 1, NSTREAMS
!        DO O1 = 1, NSTOKES
!          DO K = 1, NSTKS_NSTRMS
!           write(*,'(3i3,1p2e14.6)')
!     &             2,1,K,R2_HOMP(2,1,K),R2_HOMM(2,1,K)
!          ENDDO
!        ENDDO
!      ENDDO
!      pause'r2homp'
      
!  Finish

      RETURN
      END SUBROUTINE BVP_BRDF_SURFACE_HOM

!

      SUBROUTINE BVP_BRDF_SURFACE_BEAM
     I     ( DO_INCLUDE_SURFACE,
     I       FOURIER_COMPONENT, IBEAM,
     I       SURFACE_FACTOR )

!  Include files
!  -------------

!  Additional sums for the final albedo-reflecting layer

!  include file of dimensiopns and numbers

      INCLUDE 'VLIDORT.PARS'

!  include files of input and bookkeeping variables

      INCLUDE 'VLIDORT_INPUTS.VARS'
      INCLUDE 'VLIDORT_BOOKKEEP.VARS'

!  include file of solution stuff (input to this module)

      INCLUDE 'VLIDORT_SOLUTION.VARS'

!  Output include files (surface reflectance, kernel contributions)

      INCLUDE 'VLIDORT_REFLECTANCE.VARS'
      INCLUDE 'VLIDORT_BRDF.VARS'

!  control input
!  -------------

      LOGICAL          DO_INCLUDE_SURFACE
      INTEGER          FOURIER_COMPONENT, IBEAM
      DOUBLE PRECISION SURFACE_FACTOR

!  local variables
!  ---------------

      DOUBLE PRECISION REFL_B, REFL_B_S, FACTOR_KERNEL
      INTEGER          I, J, O1, O2, OM, KL, NL

!  Initialization
!  ==============

!  Zero total reflected contributions

      DO I = 1, NSTREAMS
        DO O1 = 1, NSTOKES
          R2_BEAM(I,O1)    = ZERO
        ENDDO
      ENDDO

!  Return with Zeroed values if albedo flag not set

      IF ( .NOT. DO_INCLUDE_SURFACE ) RETURN

!  Last layer

      NL = NLAYERS

!  ==================================================
!  ==================================================
!  Integrate Downward streams of particular solutions
!  ==================================================
!  ==================================================

      DO KL = 1, N_BRDF_KERNELS

        FACTOR_KERNEL = SURFACE_FACTOR * BRDF_FACTORS(KL)

!  For Lambertian reflectance, all streams are the same
!  ----------------------------------------------------

        IF ( LAMBERTIAN_KERNEL_FLAG(KL) ) THEN

          IF ( FOURIER_COMPONENT .EQ. 0 ) THEN

!  Particular solution

            REFL_B = ZERO
            DO J = 1, NSTREAMS
              REFL_B = REFL_B + QUAD_STRMWTS(J)*WLOWER(J,1,NL)
            ENDDO
            REFL_B = REFL_B * FACTOR_KERNEL
            DO I = 1, NSTREAMS
              A_DIFFUSE_BEAM(I,1,KL) = REFL_B
              DO O1 = 2, NSTOKES
                A_DIFFUSE_BEAM(I,O1,KL) = ZERO
              ENDDO
            ENDDO

!  zero these terms otherwise

          ELSE
            DO I = 1, NSTREAMS
              DO O1 = 1, NSTOKES
                A_DIFFUSE_BEAM(I,O1,KL) = ZERO
              ENDDO
            ENDDO
          ENDIF

!  For bidirectional reflecting surface
!  ------------------------------------

        ELSE

!  Particular solution

          DO I = 1, NSTREAMS
            DO O1 = 1, NSTOKES
              REFL_B = ZERO
              DO J = 1, NSTREAMS
                REFL_B_S = ZERO
                DO O2 = 1, NSTOKES
                  OM = MUELLER_INDEX(O1,O2)
                  REFL_B_S = REFL_B_S +
     &                     AXBID(I,J,OM,KL) * WLOWER(J,O2,NL)
                ENDDO
                REFL_B = REFL_B + REFL_B_S
              ENDDO
              A_DIFFUSE_BEAM(I,O1,KL) = FACTOR_KERNEL * REFL_B
            ENDDO
          ENDDO

        ENDIF

!  End loop over albedo kernels

      ENDDO
 
!  Total terms
!  -----------

      DO KL = 1, N_BRDF_KERNELS
        DO I = 1, NSTREAMS
          DO O1 = 1, NSTOKES
            R2_BEAM(I,O1) = R2_BEAM(I,O1) + A_DIFFUSE_BEAM(I,O1,KL)
          ENDDO
        ENDDO
      ENDDO

!  debug

!      DO I = 1, NSTREAMS
!        DO O1 = 1, NSTOKES
!          write(*,'(2i3,1p2e14.6)')I,O1,R2_BEAM(I,O1)
!        ENDDO
!      ENDDO
!      pause'r2beam'
      
!  Finish

      RETURN
      END SUBROUTINE BVP_BRDF_SURFACE_BEAM

!

      SUBROUTINE BOA_BRDF_SOURCE
     I    ( DO_INCLUDE_SURFACE,
     I      DO_INCLUDE_SURFEMISS,
     I      DO_INCLUDE_DIRECTBEAM,
     I      FOURIER_COMPONENT, IBEAM,
     I      SURFACE_FACTOR )

!  Bottom of the atmosphere source term
!    For Kernel-based BRDF system

!  Include files
!  -------------

!  include file of dimensions and numbers

      INCLUDE 'VLIDORT.PARS'

!  Include files of input and bookkeeping variables

      INCLUDE 'VLIDORT_INPUTS.VARS'
      INCLUDE 'VLIDORT_BOOKKEEP.VARS'

!  include files of solution, setup & reflectance variables (input, output)

      INCLUDE 'VLIDORT_BRDF.VARS'
      INCLUDE 'VLIDORT_SOLUTION.VARS'
      INCLUDE 'VLIDORT_REFLECTANCE.VARS'
      INCLUDE 'VLIDORT_SETUPS.VARS'

!  Subroutine input arguments
!  --------------------------

!  local control flags

      LOGICAL          DO_INCLUDE_SURFACE
      LOGICAL          DO_INCLUDE_DIRECTBEAM
      LOGICAL          DO_INCLUDE_SURFEMISS

!  surface multiplier and Fourier index

      DOUBLE PRECISION SURFACE_FACTOR
      INTEGER          FOURIER_COMPONENT
      INTEGER          IBEAM

!  local variables
!  ---------------

      INTEGER          N, J, I, IR, IROW, UM, O1, O2, OM, KL
      INTEGER          K, KO1, K0, K1, K2
      DOUBLE PRECISION REFLEC, S_REFLEC, KMULT, STOTAL
      DOUBLE PRECISION SPAR, HOM1, HOM2, SHOM_R
      DOUBLE PRECISION SHOM_CR, HOM1CR, HOM2CR
      DOUBLE PRECISION LXR, MXR, LXR_CR, LXR_CI, MXR_CR

!  initialise boa source function

      DO UM = LOCAL_UM_START, N_USER_STREAMS
        DO O1 = 1, NSTOKES
          BOA_DIFFUSE_SOURCE(UM,O1) = ZERO
          BOA_DIRECT_SOURCE(UM,O1)  = ZERO
        ENDDO
      ENDDO

!  Last layer and offset

      N = NLAYERS
      KO1 = K_REAL(N) + 1

!  reflectance from surface
!  ------------------------

      IF ( DO_INCLUDE_SURFACE ) THEN

!  Downward Stokes vector at surface at computational angles
!    reflectance integrand  a(j).x(j).I(-j)

        DO I = 1, NSTREAMS
          IR = ( I - 1 ) * NSTOKES
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            SPAR = WLOWER(I,O1,N)
            SHOM_R = ZERO
            DO K = 1, K_REAL(N)
              LXR  = LCON(K,N)*SOLA_XPOS(I,O1,K,N)
              MXR  = MCON(K,N)*SOLB_XNEG(I,O1,K,N)
              HOM1 = LXR * T_DELT_EIGEN(K,N)
              HOM2 = MXR
              SHOM_R = SHOM_R + HOM1 + HOM2
            ENDDO
            SHOM_CR = ZERO
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K-2
              K1 = KO1 + K0
              K2 = K1  + 1
              LXR_CR =  LCON(K1,N) * SOLA_XPOS(I,O1,K1,N) -
     &                  LCON(K2,N) * SOLA_XPOS(I,O1,K2,N)
              LXR_CI =  LCON(K1,N) * SOLA_XPOS(I,O1,K2,N) +
     &                  LCON(K2,N) * SOLA_XPOS(I,O1,K1,N)
              MXR_CR =  MCON(K1,N) * SOLB_XNEG(I,O1,K1,N) -
     &                  MCON(K2,N) * SOLB_XNEG(I,O1,K2,N)
              HOM1CR = LXR_CR*T_DELT_EIGEN(K1,N)
     &                -LXR_CI*T_DELT_EIGEN(K2,N)
              HOM2CR = MXR_CR
              SHOM_CR = SHOM_CR + HOM1CR + HOM2CR
            ENDDO
            STOTAL = SPAR + SHOM_R + SHOM_CR
            STOKES_DOWNSURF(I,O1) = QUAD_STRMWTS(I) * STOTAL
          ENDDO
        ENDDO

!  reflected multiple scatter intensity at user defined-angles
!  -----------------------------------------------------------

!  start loop over BRDF kernels

        DO KL = 1, N_BRDF_KERNELS

          KMULT = SURFACE_FACTOR * BRDF_FACTORS(KL)

!  ###### Lambertian reflectance (same for all user-streams)

          IF ( LAMBERTIAN_KERNEL_FLAG(KL) ) THEN

            IF ( FOURIER_COMPONENT .EQ. 0 ) THEN
              O1 = 1
              REFLEC = ZERO
              DO J = 1, NSTREAMS
                REFLEC = REFLEC + STOKES_DOWNSURF(J,O1)
              ENDDO
              REFLEC = KMULT * REFLEC
              DO UM = LOCAL_UM_START, N_USER_STREAMS
                A_USER_DIFFUSE_SOURCE(UM,O1,KL) = REFLEC
              ENDDO
              DO O1 = 2, NSTOKES
                DO UM = LOCAL_UM_START, N_USER_STREAMS
                  A_USER_DIFFUSE_SOURCE(UM,O1,KL) = ZERO
                ENDDO
              ENDDO
            ELSE
              DO O1 = 1, NSTOKES
                DO UM = LOCAL_UM_START, N_USER_STREAMS
                  A_USER_DIFFUSE_SOURCE(UM,O1,KL) = ZERO
                ENDDO
              ENDDO
            ENDIF

!  ###### bidirectional reflectance

          ELSE

!  Code activated 7 May 2005.

            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO O1 = 1, NSTOKES
                REFLEC = ZERO
                DO J = 1, NSTREAMS
                  S_REFLEC = ZERO
                  DO O2 = 1, NSTOKES
                    OM = MUELLER_INDEX(O1,O2)
                    S_REFLEC = S_REFLEC + STOKES_DOWNSURF(J,O2) *
     &                                  USER_BIREFLEC(OM,KL,UM,J)
                  ENDDO
                  REFLEC = REFLEC + S_REFLEC
                ENDDO
                A_USER_DIFFUSE_SOURCE(UM,O1,KL) = KMULT * REFLEC
              ENDDO
            ENDDO

          ENDIF

!  End loop over albedo kernels

        ENDDO

!  Total BOA Diffuse source term

        DO KL = 1, N_BRDF_KERNELS
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO O1 = 1, NSTOKES
              BOA_DIFFUSE_SOURCE(UM,O1) = BOA_DIFFUSE_SOURCE(UM,O1) +
     &           A_USER_DIFFUSE_SOURCE(UM,O1,KL)
            ENDDO
          ENDDO
        ENDDO

!  Add direct beam if flagged

        IF ( DO_INCLUDE_DIRECTBEAM ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO O1 = 1, NSTOKES
              BOA_DIRECT_SOURCE(UM,O1) =
     &               BOA_DIRECT_SOURCE(UM,O1)     +
     &               USER_DIRECT_BEAM(UM,IBEAM,O1)
            ENDDO
          ENDDO
        ENDIF

!  End inclusion of surface terms

      ENDIF

!  Add surface emission term if flagged

      IF ( DO_INCLUDE_SURFEMISS ) THEN
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO O1 = 1, NSTOKES
            BOA_DIFFUSE_SOURCE(UM,O1) = BOA_DIFFUSE_SOURCE(UM,O1) +
     &               USER_EMISSIVITY(UM,O1)
          ENDDO
        ENDDO
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE BOA_BRDF_SOURCE

!

      SUBROUTINE VLIDORT_DBCORRECTION (FLUXMULT)

!  Prepares Exact Direct Beam reflection for the BRDF case

!  Include files
!  -------------

!  include file of dimensions and numbers

      INCLUDE 'VLIDORT.PARS'

!  Include files of input and bookkeeping variables 

      INCLUDE 'VLIDORT_INPUTS.VARS'
      INCLUDE 'VLIDORT_BOOKKEEP.VARS'

!  Include file of setup variables (Input to the present module)

      INCLUDE 'VLIDORT_SETUPS.VARS'

!  Input arguments controlling type of surface

      INCLUDE 'VLIDORT_REFLECTANCE.VARS'

!  output arguments also go in here.

      INCLUDE 'VLIDORT_BRDF.VARS'

!  Output goes in the correction file

      INCLUDE 'VLIDORT_SINGSCAT.VARS'

!  input argument

      DOUBLE PRECISION FLUXMULT

!  Local variables
!  ---------------

      INTEGER          N, NUT, NSTART, NUT_PREV, NLEVEL, KL
      INTEGER          UT, UTA, UM, UA, NC, IB, V, O1, O2, OM
      DOUBLE PRECISION FINAL_SOURCE, TR
      DOUBLE PRECISION X0_FLUX, X0_BOA, ATTN, REFL_ATTN, SUM

!  first stage
!  -----------

!  Initialize

      DO V = 1, N_GEOMETRIES
        DO O1 = 1, NSTOKES
          EXACTDB_SOURCE(V,O1) = ZERO
          DO UTA = 1, N_USER_LEVELS
            STOKES_DB(UTA,V,O1)      = ZERO
          ENDDO
        ENDDO
      ENDDO

!  return if no upwelling

      IF ( .NOT.DO_UPWELLING ) RETURN

!  Reflection of solar beam
!  ------------------------

!  Must use adjusted values here.

! Old code --------------------------------------------------
!  Code to deal with refractive geometry case
!   R. Spurr, 7 May 2005. RT Solutions Inc.
!      DO IB = 1, NBEAMS
!        IF ( DO_REFLECTED_DIRECTBEAM(IB) ) THEN
!  Solar attenuation
!          IF ( DO_REFRACTIVE_GEOMETRY ) THEN
!            X0_BOA = DCOS(SZA_LOCAL_INPUT(NLAYERS,IB)*DEG_TO_RAD)
!          ELSE
!            X0_BOA = X0(IB)
!          ENDIF 
!c          X0_FLUX        = FOUR * X0_BOA * FLUX_FACTOR
!          X0_FLUX        = FOUR * X0_BOA
!          ATTN           = X0_FLUX * SOLAR_BEAM_OPDEP(IB)
!          ATMOS_ATTN(IB) = ATTN
!  Start loop over albedo kernels
!          DO KL = 1, N_BRDF_KERNELS
!           REFL_ATTN = BRDF_FACTORS(KL) * ATTN
!           DO UM = 1, N_USER_STREAMS
!            DO UA = 1, N_USER_RELAZMS
!             V = UMOFF(IB,UM) + UA
!             DO O1 = 1, NSTOKES
!              SUM = ZERO
!              DO O2 = 1, NSTOKES
!               OM = MUELLER_INDEX(O1,O2)
!               SUM = SUM + FLUXVEC(O2)*EXACTDB_BRDFUNC(OM,KL,UM,UA,IB)
!              ENDDO
!              A_EXACTDB_SOURCE(V,O1,KL) = REFL_ATTN * SUM
!             ENDDO
!            ENDDO
!           ENDDO
!          ENDDO
!  Finish loop over beams
!
!        ENDIF
!      ENDDO

!  New Code  R. Spurr, 6 August 2007. RT Solutions Inc.
!  ====================================================

!  Start geoemtry loops

      DO UM = 1, N_USER_STREAMS
        DO IB = 1, NBEAMS
          IF ( DO_REFLECTED_DIRECTBEAM(IB) ) THEN
            DO UA = 1, N_USER_RELAZMS
              V = VZA_OFFSETS(IB,UM) + UA

!  Beam attenuation

              IF ( DO_SSCORR_OUTGOING ) THEN
                X0_BOA = DCOS(SZANGLES_ADJUST(UM,IB,UA)*DEG_TO_RAD)
                ATTN = BOA_ATTN(V)
              ELSE
                IF ( DO_REFRACTIVE_GEOMETRY ) THEN
                  X0_BOA = DCOS(SZA_LOCAL_INPUT(NLAYERS,IB)*DEG_TO_RAD)
                ELSE
                  X0_BOA = COS_SZANGLES(IB)
                ENDIF 
                ATTN = SOLAR_BEAM_OPDEP(IB)
              ENDIF
              X0_FLUX = FOUR * X0_BOA
              ATTN    = ATTN * X0_FLUX
              ATTN_DB_SAVE(V) = ATTN

!  Loop over albedo kernels, assign reflection

              DO KL = 1, N_BRDF_KERNELS
                REFL_ATTN = BRDF_FACTORS(KL) * ATTN
                DO O1 = 1, NSTOKES
                  SUM = ZERO
                  DO O2 = 1, NSTOKES
                    OM = MUELLER_INDEX(O1,O2)
                    SUM = SUM +
     &                FLUXVEC(O2)*EXACTDB_BRDFUNC(OM,KL,UM,UA,IB)
                  ENDDO
                  A_EXACTDB_SOURCE(V,O1,KL) = REFL_ATTN * SUM
                ENDDO
              ENDDO

!  Finish loops over geometries

            ENDDO
          ENDIF
        ENDDO
      ENDDO

!  Build together the kernels. Multiply result by solar Flux

      DO V = 1, N_GEOMETRIES
        DO O1 = 1, NSTOKES
          DO KL = 1, N_BRDF_KERNELS
            EXACTDB_SOURCE(V,O1) = EXACTDB_SOURCE(V,O1) +
     &                   A_EXACTDB_SOURCE(V,O1,KL)
          ENDDO
          EXACTDB_SOURCE(V,O1) = EXACTDB_SOURCE(V,O1) * FLUXMULT
        ENDDO
      ENDDO

!  Upwelling recurrence: transmittance of exact source term
!  --------------------------------------------------------

!  initialize cumulative source term
!    (Already flux-multiplied, because EXACTBD is a final result)

      NC =  0
      DO V = 1, N_GEOMETRIES
        DO O1 = 1, NSTOKES
          DB_CUMSOURCE(V,O1,NC) = EXACTDB_SOURCE(V,O1)
        ENDDO
      ENDDO

!  initialize optical depth loop

      NSTART = NLAYERS
      NUT_PREV = NSTART + 1

!  Main loop over all output optical depths

      DO UTA = N_USER_LEVELS, 1, -1

!  Layer index for given optical depth

        NLEVEL = UTAU_LEVEL_MASK_UP(UTA)
        NUT    = NLEVEL + 1

!  Cumulative layer transmittance :
!    loop over layers working upwards to level NUT

        DO N = NSTART, NUT, -1
          NC = NLAYERS + 1 - N

          DO IB = 1, NBEAMS
            DO UM = 1, N_USER_STREAMS
              DO UA = 1, N_USER_RELAZMS
                V = VZA_OFFSETS(IB,UM) + UA
                IF ( DO_SSCORR_OUTGOING ) THEN
                  TR = UP_LOSTRANS(N,V)
                ELSE
                  TR = T_DELT_USERM(N,UM)
                ENDIF
                DO O1 = 1, NSTOKES
                  DB_CUMSOURCE(V,O1,NC) = TR * DB_CUMSOURCE(V,O1,NC-1)
                ENDDO
              ENDDO
            ENDDO
          ENDDO

!  end layer loop

        ENDDO

!  Offgrid output : partial layer transmittance, then set result

        IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

          UT = PARTLAYERS_OUTINDEX(UTA)
          N  = PARTLAYERS_LAYERIDX(UT)
          DO IB = 1, NBEAMS
            DO UM = 1, N_USER_STREAMS
              DO UA = 1, N_USER_RELAZMS
                V = VZA_OFFSETS(IB,UM) + UA
                IF ( DO_SSCORR_OUTGOING ) THEN
                  TR = UP_LOSTRANS_UT(UT,V)
                ELSE
                  TR = T_UTUP_USERM(UT,UM)
                ENDIF
                DO O1 = 1, NSTOKES
                  STOKES_DB(UTA,V,O1) = TR * DB_CUMSOURCE(V,O1,NC)
                ENDDO
              ENDDO
            ENDDO
          ENDDO

!  Ongrid output : Set final cumulative source directly

        ELSE

          DO IB = 1, NBEAMS
            DO UM = 1, N_USER_STREAMS
              DO UA = 1, N_USER_RELAZMS
                V = VZA_OFFSETS(IB,UM) + UA
                DO O1 = 1, NSTOKES
                  FINAL_SOURCE = DB_CUMSOURCE(V,O1,NC)
                  STOKES_DB(UTA,V,O1) = FINAL_SOURCE
               ENDDO
              ENDDO
            ENDDO
          ENDDO

        ENDIF

!  Check for updating the recursion 

        IF ( NUT. NE. NUT_PREV ) NSTART = NUT - 1
        NUT_PREV = NUT

!  end optical depth loop

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_DBCORRECTION


end module m_vlidort_brdf_supplement
