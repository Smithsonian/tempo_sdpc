MODULE parameters_module

  
  ! -------------
  ! Default Sizes
  ! -------------
  INTEGER, PARAMETER :: maxChar = 512 ! Default maximum string length
  INTEGER, PARAMETER :: MaxProfPar = 5 ! Maximum number of parameters controlling a parameterized profile

  ! -------------
  ! Unit Scalings
  ! -------------

  ! This is to scale the photons/cm2/s/nm units to O(1)
  REAL(KIND=8), PARAMETER :: PhotonScalingUnit = 1.0d14

  ! ----------
  ! Precisions
  ! ----------
  
  ! =====================================================
  ! Define KIND variables for single and double precision
  ! =====================================================
  INTEGER, PARAMETER :: i1 = SELECTED_INT_KIND(2**1)
  INTEGER, PARAMETER :: i2 = SELECTED_INT_KIND(2**2)
  INTEGER, PARAMETER :: i4 = SELECTED_INT_KIND(2**3)
  INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND(2**4)
  INTEGER, PARAMETER :: r4 = KIND(1.0)
  INTEGER, PARAMETER :: r8 = KIND(1.0D0)

  ! ------------------
  ! Physical Constants
  ! ------------------
  REAL(KIND=8), PARAMETER :: Constants_R  = 8.3144598d0     ! Universal Gas Constant (J/K/mol)
  REAL(KIND=8), PARAMETER :: Constants_NA = 6.0221409d23    ! Avogadros constant (g/mol)
  REAL(KIND=8), PARAMETER :: Constants_c  = 2.99792458d2    ! Speed of light (m/s)
  REAL(KIND=8), PARAMETER :: Constants_kb = 1.38064852d-23  ! Boltzmann Constant (J/K)
  REAL(KIND=8), PARAMETER :: Constants_h  = 6.626070040d-34 ! Planck Constant (J.s)
  REAL(KIND=8), PARAMETER :: Constants_pi = 3.14159265358979d0
  REAL(KIND=8), PARAMETER :: deg2rad = Constants_pi / 180.d0
  REAL(KIND=8), PARAMETER :: rad2deg = 180.d0 / Constants_pi

  ! ---------------------------------------------------------------
  ! Reserved file unit numbers
  ! ---------------------------------------------------------------
  INTEGER, PARAMETER :: ctrunit = 11 ! Input control file
  INTEGER, PARAMETER :: srfunit = 12 ! Surface reading
  INTEGER, PARAMETER :: eofunit = 14
  INTEGER, PARAMETER :: lckunit = 15
  INTEGER, PARAMETER :: htrunit = 22 ! Hitran file access
  INTEGER, PARAMETER :: errunit = 16 ! Error Log
  INTEGER, PARAMETER :: pcaunit = 17 ! PCA Control files
  
  ! ---------------------------------------------------------------
  ! VLIDORT BRDF Kernel Stuff
  ! ---------------------------------------------------------------

! The new supplement in v2p8 has more kernels
  ! Short names for brdf kernels (used for output)
  CHARACTER(LEN=5), DIMENSION(16), PARAMETER ::  short_kern_name =(/  &
                'isotr', 'rthin', 'rthck', 'lsprs', 'ldens', 'hapke', &
                'rjean', 'rhman', 'cxmnk', 'gcmnk', 'gcmcr', 'bpdfv', &
                'bpdfs', 'bpdfn', 'ncmgl', 'hpke5'/)

  ! Full BRDF names checked in VLIDORT supplement
  CHARACTER (LEN=10), DIMENSION(16), PARAMETER :: BRDF_CHECK_NAMES = (/ &
    'Lambertian', 'Ross-thin ', 'Ross-thick', 'Li-sparse ', 'Li-dense  ', 'Hapke     ', &
    'Roujean   ', 'Rahman    ', 'Cox-Munk  ', 'GissCoxMnk', 'GCMcomplex', 'BPDF-Vegn ', &
    'BPDF-Soil ', 'BPDF-NDVI ', 'NewCMGlint', 'Hapke5Par '/)

  ! Number of parameters for each kernel
  INTEGER, DIMENSION(16), PARAMETER :: vl_brdf_npar = &
                                      (/0,0,0,2,2,3,0,3,2,2,2,2,2,2,3,5/)

  ! ---------------------------------------------------------------
  ! VLIDORT Phase function indices
  ! ---------------------------------------------------------------
  
  INTEGER, DIMENSION(8), PARAMETER :: greekmat_idxs = (/1, 2, 5, 6, 11, 12, 15, 16/)
  INTEGER, DIMENSION(8), PARAMETER :: phasmoms_idxs = (/1, 5, 5, 2,  3,  6,  6,  4/)

  ! ---------------------------------------------------------------
  ! Some default types
  ! ---------------------------------------------------------------
  
  TYPE DiagSpcOpt
    LOGICAL                             :: DoDiag
    LOGICAL                             :: DoQU
    INTEGER                             :: nSpc
    CHARACTER(LEN=maxChar), ALLOCATABLE :: Name(:)
    INTEGER,                ALLOCATABLE :: Idx(:)
  ENDTYPE DiagSpcOpt
  
  TYPE ParJacType
    INTEGER                             :: ProfileType
    INTEGER                             :: nPar
    CHARACTER(LEN=maxChar), ALLOCATABLE :: ParName(:)
    REAL(KIND=8),           ALLOCATABLE :: K(:,:,:)
  ENDTYPE ParJacType
  
  ! Program Error Codes
  INTEGER(KIND=2), PARAMETER :: ErrorCode_Input   = 1
  INTEGER(KIND=2), PARAMETER :: ErrorCode_FileIO  = 2
  INTEGER(KIND=2), PARAMETER :: ErrorCode_RTM     = 3
  INTEGER(KIND=2), PARAMETER :: ErrorCode_OptProp = 4
  INTEGER(KIND=2), PARAMETER :: ErrorCode_Profile = 5
  INTEGER(KIND=2), PARAMETER :: ErrorCode_Level1  = 6
  INTEGER(KIND=2), PARAMETER :: ErrorCode_SpecFit = 7


! ============================================================================
! The following is from the parameters module of the fparser code
! ============================================================================
!
! Copyright (c) 2000-2008, Roland Schmehl. All rights reserved.
!
! This software is distributable under the BSD license. See the terms of the
! BSD license in the documentation provided with this software.
INTEGER, PARAMETER :: rn = KIND(0.0d0)          ! Precision of real numbers
INTEGER, PARAMETER :: is = SELECTED_INT_KIND(1) ! Data type of bytecode


! For some things that need flexible arrays
TYPE ArrType_2D
  REAL(KIND=8), ALLOCATABLE :: Arr(:,:)
ENDTYPE ArrType_2D

TYPE ArrType_3D
  REAL(KIND=8),   ALLOCATABLE :: Arr(:,:,:)
  INTEGER(KIND=2),ALLOCATABLE :: Arr_i2(:,:,:)
ENDTYPE ArrType_3D

TYPE ArrType_4D
  REAL(KIND=8), ALLOCATABLE :: Arr(:,:,:,:)
ENDTYPE ArrType_4D

TYPE ArrSetType_2D
  INTEGER                       :: n
  INTEGER,          ALLOCATABLE :: imx(:), jmx(:)
  TYPE(ArrType_2D), ALLOCATABLE :: Data(:)
ENDTYPE ArrSetType_2D

TYPE ArrSetType_3D
  INTEGER                       :: n
  INTEGER,          ALLOCATABLE :: imx(:), jmx(:), kmx(:)
  TYPE(ArrType_3D), ALLOCATABLE :: Data(:)
ENDTYPE ArrSetType_3D

TYPE ArrSetType_4D
  INTEGER                       :: n
  INTEGER,          ALLOCATABLE :: imx(:), jmx(:), kmx(:), lmx(:)
  TYPE(ArrType_4D), ALLOCATABLE :: Data(:)
ENDTYPE ArrSetType_4D

END MODULE parameters_module
