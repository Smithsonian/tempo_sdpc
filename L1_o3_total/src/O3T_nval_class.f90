!!****************************************************************************
!!F90
!
!!Description:
!
!  MODULE O3T_nval_class
! 
! read in tables for calculation of forward model quantity: n-value
!
!!Input Parameters:
! None
!
!!Output Parameters:
! None
! 
!!Return
! None 
!
!!Revision History:
! Initial version 03/26/2002  Kai Yang/UMBC
!
!!Team-unique Header:
! This software was developed by the OMI Science Team Support
! Group for the National Aeronautics and Space Administration, Goddard
! Space Flight Center, under NASA Task 916-003-1
!
!!References and Credits
! Written by 
! Kai Yang 
! University of Maryland Baltimore County
! email: Kai.Yang-1@nasa.gov
! 
!!Design Notes
!
!!END
!!****************************************************************************
MODULE O3T_nval_class
  USE PGS_PC_class
  USE OMI_LUN_set
  USE OMI_SMF_class

  IMPLICIT NONE
  PUBLIC
  INTEGER (KIND = 4), PARAMETER, PRIVATE :: zero = 0, one = 1, two = 2
  INTEGER (KIND = 4), PARAMETER, PRIVATE :: three = 3, four = 4, five = 5
  INTEGER (KIND = 4), PARAMETER, PRIVATE :: np_nval_read = 13

  TYPE, PUBLIC :: nvalLUT_t
    CHARACTER( LEN=256 ) :: nval_fn
    LOGICAL :: nval_read, ringAdded
    INTEGER(KIND=4) :: nvza
    INTEGER(KIND=4) :: nsza
    INTEGER(KIND=4) :: nwl
    INTEGER(KIND=4) :: npres
    INTEGER(KIND=4) :: nprof
    INTEGER(KIND=4), DIMENSION(5) :: n5size
    INTEGER(KIND=4), DIMENSION(3) :: n3size
    REAL(KIND=4) :: wl_min, wl_max
    REAL(KIND=4), DIMENSION(:), POINTER :: sza, vza, pres, wl
    REAL(KIND=8), DIMENSION(:), POINTER :: presLog
    REAL(KIND=4), DIMENSION(:), POINTER :: c0, c1, c2
    REAL(KIND=4), DIMENSION(:), POINTER :: lgi0, z1, &
         z2, tr, sb, knb
  END TYPE nvalLUT_t

  TYPE(nvalLUT_t), PUBLIC :: nvELS, nvRRS

  REAL (KIND = 4), PUBLIC :: wlT_min,   wl_mix,    wl_oz, &
       wl_refl_l, wl_refl_h, wlT_max

  LOGICAL :: nval_read
  INTEGER (KIND = 4), PUBLIC :: nvzaT
  INTEGER (KIND = 4), PUBLIC :: nszaT
  INTEGER (KIND = 4), PUBLIC :: nwl_com       ! the number of wl stored here
  INTEGER (KIND = 4), PUBLIC :: nwlLUT        ! the number of wl in HDF tab     
  INTEGER (KIND = 4), PUBLIC :: npresT
  INTEGER (KIND = 4), PUBLIC :: nprofT 
  INTEGER (KIND = 4), PUBLIC :: iwl_refl_l, iwl_refl_h, iwl_ozone, iwl_mix 
  INTEGER (KIND = 4), PUBLIC :: iwl_s, iwl_e  ! refer to wl_t first then
  ! refer to wl at the end
  INTEGER (KIND = 4), DIMENSION(1), PUBLIC :: RingAdded !0 for FALSE, 1 for TRUE

  !data grids
  REAL (KIND = 4), DIMENSION(:), ALLOCATABLE:: sza, vza, pressure, wl_com
  REAL (KIND = 8), DIMENSION(:), ALLOCATABLE, PUBLIC:: log_pres

  ! c0, c1, c2 is 1-D array with size equal to nwl_com
  REAL (KIND = 4), DIMENSION(:), ALLOCATABLE, PUBLIC:: c0, c1, c2 

  !
  ! For all the standard profiles at layer 0 in the h5 files
  ! lgi0, z1, z2, tr are arrays with a shape as 
  ! (nvzaT, nszaT, nprofT, nwl_com, npresT ) 
  ! sb is an array with shape ( nprofT, nwl_com, npresT )
  !
  INTEGER (KIND = 4), DIMENSION(5) :: n5size  
  INTEGER (KIND = 4), DIMENSION(3) :: n3size  
  REAL (KIND = 4), DIMENSION(:), ALLOCATABLE:: lgi0_n, z1_n, &
       z2_n, tr_n, sb_n, knb_n

  PUBLIC  :: O3T_nval_setup
  !  PUBLIC  :: O3T_nval_setup_omps
  PUBLIC  :: O3T_nval_dispose
  PUBLIC  :: O3T_getLambdaSet
  PRIVATE :: O3T_getLambdaRange

CONTAINS  
  FUNCTION O3T_nval_setup( nval_lun, nv ) RESULT( status )
    USE HDF5
    USE O3T_LUT_fs

    CHARACTER ( LEN = 200 ) :: fn_nval

    INTEGER (KIND = 4), INTENT(IN) :: nval_lun
    TYPE(nvalLUT_t), INTENT(OUT), OPTIONAL :: nv

    INTEGER (HID_T) :: file_id
    INTEGER (KIND = 4) :: status
    INTEGER :: ip_nval, error, ierr, i_foo(1) ! ie, iprof, is,
    !INTEGER :: iwl, ipres, isza, ivza 
    !INTEGER, DIMENSION(5) :: data_dims
    TYPE (DSH5_T), DIMENSION(np_nval_read) :: nval_parameters
    INTEGER, DIMENSION(np_nval_read) :: rank_nval =  &
         (/ 1,1,1,5,5,5,5,5,3,1,1,1,1 /)
    REAL (KIND = 4), DIMENSION(:,:,:,:,:), ALLOCATABLE :: &
         lgi0_t, z1_t, z2_t, tr_t, knb_t
    REAL (KIND = 4), DIMENSION(:,:,:), ALLOCATABLE :: sb_t
    REAL (KIND = 4), DIMENSION(:), ALLOCATABLE :: wl_t, c0_t, c1_t, c2_t 
    CHARACTER ( LEN = 200 ) :: msg
    CHARACTER ( LEN = 20 ) :: FUNCTIONNAME = "O3T_nval_setup"
    !INTEGER (HID_T)  :: datatype
    INTEGER :: version

    !! get the nval table filename and path
    version = 1
    status = PGS_PC_GetReference( nval_lun, version, fn_nval )
    IF( status /= PGS_S_SUCCESS ) THEN
      WRITE( msg,'(A,I9)' ) "get LUT name failed at LUN = ", nval_lun
      ierr = OMI_SMF_setmsg( status, msg, FUNCTIONNAME, zero )
      status = OZT_E_FAILURE
      nval_read = .FALSE.
      RETURN
    ENDIF

    nval_parameters  = (/ ds_sza, ds_vza, ds_pres, ds_lgi0, &
         ds_z1,  ds_z2,  ds_tr,   ds_knb, ds_sb, &
         ds_wl,  ds_c0,  ds_c1,   ds_c2 /)

    ! Initialize FORTRAN interface.
    CALL h5open_f(error)

    CALL h5fopen_f( fn_nval, H5F_ACC_RDONLY_F, file_id, error )
    IF( error < zero ) THEN
      status = OZT_E_FAILURE
      WRITE( msg,* ) "h5fopen_f failed on file: ", fn_nval
      ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), FUNCTIONNAME, zero )
      nval_read = .FALSE.
      RETURN
    ENDIF

    RingAdded(:) = 0 ! initialized to no LRRS ring
    ds_RingAdded%datatype = H5T_NATIVE_INTEGER
    ds_RingAdded%dims(1)  = 1
    CALL h5tget_class_f(ds_RingAdded%datatype, ds_RingAdded%class, error )
    status = UTIL_lh5_selectDS( file_id, ds_RingAdded )
    status = UTIL_lh5_get( ds_RingAdded, RingAdded )
    status = UTIL_lh5_disposeDS( ds_RingAdded )

    DO ip_nval = 1, np_nval_read
      !! if no ring added, skip knr2
      IF( RingAdded(1) == 0 .AND. &
           INDEX( nval_parameters(ip_nval)%name, "knr2")>0 ) CYCLE

      nval_parameters(ip_nval)%datatype = H5T_NATIVE_REAL
      CALL h5tget_class_f(nval_parameters(ip_nval)%datatype, &
           nval_parameters(ip_nval)%class, error )
      IF( error < zero ) THEN
        status = OZT_E_FAILURE
        WRITE( msg,* ) "h5tget_class_f failed on dataset: ", &
             nval_parameters(ip_nval)%name
        ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
             FUNCTIONNAME, zero )
        nval_read = .FALSE.
        RETURN
      ENDIF
      nval_parameters(ip_nval)%rank = rank_nval(ip_nval)
      status = UTIL_lh5_selectDS( file_id, nval_parameters(ip_nval) )
      IF( status .NE. OZT_S_SUCCESS ) THEN
        status = OZT_E_FAILURE
        WRITE( msg,* ) "UTIL_lh5_selectDS failed for ", &
             nval_parameters(ip_nval)%name, "in file ", fn_nval
        ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
             FUNCTIONNAME, zero )
        nval_read = .FALSE.
        RETURN
      ENDIF
    ENDDO

    nvzaT  = int(nval_parameters(4)%dims(1), kind=4)
    nszaT  = int(nval_parameters(4)%dims(2), kind=4)
    nprofT = int(nval_parameters(4)%dims(3), kind=4)
    nwlLUT = int(nval_parameters(4)%dims(4), kind=4)
    npresT = int(nval_parameters(4)%dims(5), kind=4)

    ! storeage space for the all the layers, all profiles, all 
    ! wavelengths, all pressures, all solar and view angles
    ALLOCATE( lgi0_t( nvzaT, nszaT, nprofT, nwlLUT, npresT ), &
         z1_t(   nvzaT, nszaT, nprofT, nwlLUT, npresT ), &
         z2_t(   nvzaT, nszaT, nprofT, nwlLUT, npresT ), &
         tr_t(   nvzaT, nszaT, nprofT, nwlLUT, npresT ), &
         knb_t(  nvzaT, nszaT, nprofT, nwlLUT, npresT ), &
         sb_t(                 nprofT, nwlLUT, npresT ), &
         wl_t(nwlLUT), c0_t(nwlLUT), c1_t(nwlLUT), c2_t(nwlLUT), &  
         sza(nszaT), vza(nvzaT), &
         pressure(npresT), log_pres(npresT), &
         STAT = ierr )
    IF( ierr /= zero ) THEN
      status = OZT_E_FAILURE
      ierr = OMI_SMF_setmsg( OZT_E_MEM_ALLOC, ".. allocation failure", &
           FUNCTIONNAME, zero )
      nval_read = .FALSE.
      RETURN
    ENDIF

    !! read the data stored in the NVAL....h5 file
    status = UTIL_lh5_get( nval_parameters(1),      sza )
    status = UTIL_lh5_get( nval_parameters(2),      vza )
    status = UTIL_lh5_get( nval_parameters(3), pressure )
    status = UTIL_lh5_get( nval_parameters(4),   lgi0_t )
    status = UTIL_lh5_get( nval_parameters(5),     z1_t )
    status = UTIL_lh5_get( nval_parameters(6),     z2_t )
    status = UTIL_lh5_get( nval_parameters(7),     tr_t )
    IF( RingAdded(1) == 1 ) THEN
      status = UTIL_lh5_get( nval_parameters(8),    knb_t )
    ENDIF
    status = UTIL_lh5_get( nval_parameters(9),     sb_t )
    status = UTIL_lh5_get( nval_parameters(10),    wl_t )
    status = UTIL_lh5_get( nval_parameters(11),    c0_t )
    status = UTIL_lh5_get( nval_parameters(12),    c1_t )
    status = UTIL_lh5_get( nval_parameters(13),    c2_t )
    log_pres = LOG10( DBLE( pressure ) )
    wl_t  = 0.1*wl_t !! convert A to nm

    DO ip_nval = 1, np_nval_read
      !! if no ring added, skip knr2
      IF( RingAdded(1) == 0 .AND. &
           INDEX( nval_parameters(ip_nval)%name, "knr2" )>0 ) CYCLE
      status = UTIL_lh5_disposeDS( nval_parameters(ip_nval ) )
    ENDDO
    CALL h5fclose_f( file_id, error)
    CALL h5close_f(error)

    status = O3T_getLambdaRange()  !! this sets wlT_min & wlT_max
    IF( status /= OZT_S_SUCCESS ) THEN
      status = OZT_E_FAILURE
      ierr = OMI_SMF_setmsg( OZT_E_INPUT, "O3T_getLambdaRange failed.", &
           FUNCTIONNAME, zero )
      nval_read = .FALSE.
      RETURN
    ENDIF

    IF( wlT_min > MAXVAL( wl_t ) ) THEN
      status = OZT_E_FAILURE
      ierr = OMI_SMF_setmsg( OZT_E_INPUT, "wlT_min too large", &
           FUNCTIONNAME, zero )
      nval_read = .FALSE.
      RETURN
    ELSE IF( wlT_min < MINVAL( wl_t )) THEN
      iwl_s = 1
    ELSE
      i_foo = MAXLOC( wl_t, MASK = wl_t-wlT_min  <= 0.0 )
      iwl_s = i_foo( 1 )
    ENDIF

    IF( wlT_max < MINVAL( wl_t ) ) THEN
      status = OZT_E_FAILURE
      ierr = OMI_SMF_setmsg( OZT_E_INPUT, "wlT_max too small", &
           FUNCTIONNAME, zero )
      nval_read = .FALSE.
      RETURN
    ELSE IF( wlT_max > MAXVAL( wl_t )) THEN
      iwl_e = nwlLUT
    ELSE
      i_foo = MINLOC( wl_t, MASK = wl_t-wlT_max  >= 0.0 )
      iwl_e = i_foo( 1 )
    ENDIF
    nwl_com = iwl_e - iwl_s + 1

    WRITE( msg,'(4(A,I5))') "wavelength used in nval calculations:iwl_s =", &
         iwl_s,",iwl_e =",iwl_e, &
         ",nwl_com =", nwl_com, ",nwlLUT =", nwlLUT
    ierr = OMI_SMF_setmsg(OZT_S_SUCCESS, msg, FUNCTIONNAME, four ) 

    n5size = (/ nvzaT*nszaT*nprofT*nwl_com*npresT, &
         nvzaT*nszaT*nprofT*nwl_com      , &
         nvzaT*nszaT*nprofT          , &
         nvzaT*nszaT                , &
         nvzaT                      /)

    n3size = (/ nprofT*nwl_com*npresT, &
         nprofT*nwl_com      , &
         nprofT            /)

    ALLOCATE( wl_com(nwl_com), c0(nwl_com), c1(nwl_com), c2(nwl_com), &  
         lgi0_n( n5size(1) ), &
         z1_n(   n5size(1) ), &
         z2_n(   n5size(1) ), &
         tr_n(   n5size(1) ), &
         knb_n(  n5size(1) ), &            
         sb_n(   n3size(1) ), &            
         STAT = ierr )
    ! pick out the unpeturbed layer (layer 0) for all the profiles
    ! this will be stored in n values for radiance calculation
    ! HDF used 0 based index, Fortran use 1 based index
    lgi0_n(:) = RESHAPE( lgi0_t(:,:,:,iwl_s:iwl_e,:), (/ n5size(1) /) )
    z1_n(:)   = RESHAPE(   z1_t(:,:,:,iwl_s:iwl_e,:), (/ n5size(1) /) )
    z2_n(:)   = RESHAPE(   z2_t(:,:,:,iwl_s:iwl_e,:), (/ n5size(1) /) )
    tr_n(:)   = RESHAPE(   tr_t(:,:,:,iwl_s:iwl_e,:), (/ n5size(1) /) )
    IF( RingAdded(1) == 1 ) THEN
      knb_n(:) = RESHAPE(knb_t(:,:,:,iwl_s:iwl_e,:), (/ n5size(1) /) )
    ELSE
      knb_n(:)  = 0.0
    ENDIF
    sb_n(:)   = RESHAPE(   sb_t(    :,iwl_s:iwl_e,:), (/ n3size(1) /) )
    wl_com(:) =            wl_t(      iwl_s:iwl_e )
    c0(:)     =            c0_t(      iwl_s:iwl_e )
    c1(:)     =            c1_t(      iwl_s:iwl_e )
    c2(:)     =            c2_t(      iwl_s:iwl_e )

    !!OPEN( UNIT=1,FILE='NVAL.txt',FORM='FORMATTED',STATUS='NEW')
    !!WRITE(1,*) lgi0_t(:,:,:,(/1,4,7,8,9,10/),:)
    !!WRITE(1,*) lgi0_t(:,:,:,:,:)
    !!CLOSE(1)
    nval_read = .TRUE.

    IF( PRESENT( nv ) ) THEN
      nv%nval_fn   = fn_nval
      nv%nval_read = nval_read

      IF( RingAdded(1) == 0  ) THEN
        nv%ringAdded = .FALSE.
      ELSE IF( RingAdded(1) == 1  ) THEN
        nv%ringAdded = .TRUE.
      ELSE
        status = OZT_E_FAILURE
        ierr = OMI_SMF_setmsg( OZT_E_INPUT, "unknown RingAdded value", &
             FUNCTIONNAME, zero )
        nv%nval_read = .FALSE.
        RETURN
      ENDIF
      nv%nvza      = nvzaT
      nv%nsza      = nszaT
      nv%nwl       = nwl_com
      nv%npres     = npresT
      nv%nprof     = nprofT
      nv%n5size(:) = n5size(:)
      nv%n3size(:) = n3size(:)
      nv%wl_min    = wlT_min
      nv%wl_max    = wlT_max
      ALLOCATE( nv%wl(nv%nwl), nv%c0(nv%nwl), nv%c1(nv%nwl), nv%c2(nv%nwl), &  
           nv%sza(nv%nsza), nv%vza(nv%nvza), nv%pres(nv%npres), &
           nv%presLog( nv%npres ), nv%lgi0( nv%n5size(1) ), &
           nv%z1(   nv%n5size(1) ), nv%z2(   nv%n5size(1) ), &
           nv%tr(   nv%n5size(1) ), nv%knb(  nv%n5size(1) ), &            
           nv%sb(   nv%n3size(1) ), STAT = ierr )
      nv%presLog(:) = log_pres(:)
      nv%sza(:)     = sza(:) 
      nv%vza(:)     = vza(:)
      nv%pres(:)    = pressure(:)
      nv%lgi0(:)    = lgi0_n(:) 
      nv%z1(:)      = z1_n(:)
      nv%z2(:)      = z2_n(:)
      nv%tr(:)      = tr_n(:)
      IF( nv%ringAdded ) THEN
        nv%knb(:)  = knb_n(:)
      ELSE
        nv%knb(:)  = 0.0
      ENDIF
      nv%sb(:)      = sb_n(:)
      nv%wl(:)      = wl_com(:)
      nv%c0(:)      = c0(:)
      nv%c1(:)      = c1(:)
      nv%c2(:)      = c2(:)
    ENDIF

    DEALLOCATE( lgi0_t, z1_t, z2_t, tr_t, sb_t, knb_t, &
         wl_t, c0_t, c1_t, c2_t, STAT = ierr )

    status = OZT_S_SUCCESS
  END FUNCTION O3T_nval_setup
  !...............................................................................

!  FUNCTION O3T_nval_setup_omps( nval_lun, nv ) RESULT( status )
!    ! 
!    ! Description: 
!    ! 
!    !   create a version of O3T_nval_setup() specifically for OMPS instrument 
!    !   where I have to deal one extra dimension: lgi0, z1, z2, tr, knb, sb
!    !   and f0flux.
!    !
!    ! Author: Jason Li  
!    !
!    ! Initial Version: 2011-APR-12 
!
!    USE HDF5
!    USE O3T_LUT_fs
!
!    CHARACTER ( LEN = 200 ) :: fn_nval
!
!    INTEGER (KIND = 4), INTENT(IN) :: nval_lun
!    TYPE(nvalLUT_t), INTENT(OUT), OPTIONAL :: nv
!
!    INTEGER (HID_T) :: file_id
!    INTEGER (KIND = 4) :: status
!    INTEGER :: nXtrackPixels 
!    INTEGER :: ip_nval, is, ie, iprof, error, ierr, i_foo(1)
!    INTEGER :: iwl, ipres, isza, ivza, middlePixelIndex
!    INTEGER, DIMENSION(5) :: data_dims
!    INTEGER(HSIZE_T), dimension(6) :: start, count
!    TYPE (DSH5_T), DIMENSION(np_nval_read) :: nval_parameters
!    !     INTEGER, DIMENSION(np_nval_read) :: rank_nval =  &
!    !                                         (/ 1,1,1,5,5,5,5,5,3,1,1,1,1 /)
!    INTEGER, DIMENSION(np_nval_read) :: rank_nval =  &
!         (/ 1,1,1,6,6,6,6,6,4,1,1,1,1 /)
!    REAL (KIND = 4), DIMENSION(:,:,:,:,:), ALLOCATABLE :: &
!         lgi0_t, z1_t, z2_t, tr_t, knb_t
!    REAL (KIND = 4), DIMENSION(:,:,:), ALLOCATABLE :: sb_t
!    REAL (KIND = 4), DIMENSION(:), ALLOCATABLE :: wl_t, c0_t, c1_t, c2_t 
!    CHARACTER ( LEN = 200 ) :: msg
!    CHARACTER ( LEN = 20 ) :: FUNCTIONNAME = "O3T_nval_setup_omps"
!    INTEGER (HID_T)  :: datatype
!    INTEGER :: version
!
!    !! get the nval table filename and path
!    version = 1
!    status = PGS_PC_GetReference( nval_lun, version, fn_nval )
!    IF( status /= PGS_S_SUCCESS ) THEN
!      WRITE( msg,'(A,I9)' ) "get LUT name failed at LUN = ", nval_lun
!      ierr = OMI_SMF_setmsg( status, msg, FUNCTIONNAME, zero )
!      status = OZT_E_FAILURE
!      nval_read = .FALSE.
!      RETURN
!    ELSE
!      print *, "[O3T_nval_setup_omps]: nvalue table = ", trim(fn_nval)
!    ENDIF
!
!    ! COMMENT-JYLI: in O3T_LUT_fs.f90, ds_knb actually points to knr2 dataset 
!    ! in the nvalue lookup table, not the knb as the name may lead you to 
!    ! believe. knb and knr2 have different ranks.  Will cause out-of-bounds 
!    ! problems if we do not pay attention!
!
!    nval_parameters  = (/ ds_sza, ds_vza, ds_pres, ds_lgi0, &
!         ds_z1,  ds_z2,  ds_tr,   ds_knb, ds_sb, &
!         ds_wl,  ds_c0,  ds_c1,   ds_c2 /)
!
!    ! Initialize FORTRAN interface.
!    CALL h5open_f(error)
!
!    CALL h5fopen_f( fn_nval, H5F_ACC_RDONLY_F, file_id, error )
!    IF( error < zero ) THEN
!      status = OZT_E_FAILURE
!      WRITE( msg,* ) "h5fopen_f failed on file: ", fn_nval
!      ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), FUNCTIONNAME, zero )
!      nval_read = .FALSE.
!      RETURN
!    ENDIF
!
!    RingAdded(:) = 0 ! initialized to no LRRS ring
!    ds_RingAdded%datatype = H5T_NATIVE_INTEGER
!    ds_RingAdded%dims(1)  = 1
!    CALL h5tget_class_f(ds_RingAdded%datatype, ds_RingAdded%class, error )
!    status = UTIL_lh5_selectDS( file_id, ds_RingAdded )
!    status = UTIL_lh5_get( ds_RingAdded, RingAdded )
!    status = UTIL_lh5_disposeDS( ds_RingAdded )
!
!    DO ip_nval = 1, np_nval_read
!      !! if no ring added, skip knr2
!      IF( RingAdded(1) == 0 .AND. &
!           INDEX( nval_parameters(ip_nval)%name, "knr2")>0 ) CYCLE
!
!      nval_parameters(ip_nval)%datatype = H5T_NATIVE_REAL
!      CALL h5tget_class_f(nval_parameters(ip_nval)%datatype, &
!           nval_parameters(ip_nval)%class, error )
!      IF( error < zero ) THEN
!        status = OZT_E_FAILURE
!        WRITE( msg,* ) "h5tget_class_f failed on dataset: ", &
!             nval_parameters(ip_nval)%name
!        ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
!             FUNCTIONNAME, zero )
!        nval_read = .FALSE.
!        RETURN
!      ENDIF
!      nval_parameters(ip_nval)%rank = rank_nval(ip_nval)
!      status = UTIL_lh5_selectDS( file_id, nval_parameters(ip_nval) )
!      IF( status .NE. OZT_S_SUCCESS ) THEN
!        status = OZT_E_FAILURE
!        WRITE( msg,* ) "UTIL_lh5_selectDS failed for ", &
!             nval_parameters(ip_nval)%name, "in file ", fn_nval
!        ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
!             FUNCTIONNAME, zero )
!        nval_read = .FALSE.
!        RETURN
!      ENDIF
!    ENDDO
!
!    nvzaT         = nval_parameters(4)%dims(1)
!    nszaT         = nval_parameters(4)%dims(2)
!    nprofT        = nval_parameters(4)%dims(3)
!    nwlLUT        = nval_parameters(4)%dims(4)
!    npresT        = nval_parameters(4)%dims(5)
!    nXtrackPixels = nval_parameters(4)%dims(6)
!    middlePixelIndex = nXtrackPixels / 2 - 1
!
!    ! storeage space for the all the layers, all profiles, all 
!    ! wavelengths, all pressures, all solar and view angles
!    ALLOCATE( lgi0_t( nvzaT, nszaT, nprofT, nwlLUT, npresT ), &
!         z1_t(   nvzaT, nszaT, nprofT, nwlLUT, npresT ), &
!         z2_t(   nvzaT, nszaT, nprofT, nwlLUT, npresT ), &
!         tr_t(   nvzaT, nszaT, nprofT, nwlLUT, npresT ), &
!         knb_t(  nvzaT, nszaT, nprofT, nwlLUT, npresT ), &
!         sb_t(                 nprofT, nwlLUT, npresT ), &
!         wl_t(nwlLUT), c0_t(nwlLUT), c1_t(nwlLUT), c2_t(nwlLUT), &  
!         sza(nszaT), vza(nvzaT), &
!         pressure(npresT), log_pres(npresT), &
!         STAT = ierr )
!    IF( ierr /= zero ) THEN
!      status = OZT_E_FAILURE
!      ierr = OMI_SMF_setmsg( OZT_E_MEM_ALLOC, ".. allocation failure", &
!           FUNCTIONNAME, zero )
!      nval_read = .FALSE.
!      RETURN
!    ENDIF
!
!    !! read the data stored in the NVAL....h5 file
!    status = UTIL_lh5_get( nval_parameters(1),      sza )
!    status = UTIL_lh5_get( nval_parameters(2),      vza )
!    status = UTIL_lh5_get( nval_parameters(3), pressure )
!
!    start = (/     0,     0,      0,      0,      0, middlePixelIndex /)
!    count = (/ nvzaT, nszaT, nprofT, nwlLUT, npresT,                1 /)
!    status = UTIL_lh5_get( nval_parameters(4),   lgi0_t, start, count )
!    status = UTIL_lh5_get( nval_parameters(5),     z1_t, start, count )
!    status = UTIL_lh5_get( nval_parameters(6),     z2_t, start, count )
!    status = UTIL_lh5_get( nval_parameters(7),     tr_t, start, count )
!    IF( RingAdded(1) == 1 ) THEN
!      status = UTIL_lh5_get( nval_parameters(8), knb_t, start, count )
!    ENDIF
!
!    start = (/      0,      0,      0, middlePixelIndex, 0, 0 /)
!    count = (/ nprofT, nwlLUT, npresT,                1, 0, 0 /)      
!    status = UTIL_lh5_get( nval_parameters(9),     sb_t, start, count )
!
!    status = UTIL_lh5_get( nval_parameters(10),    wl_t )
!    status = UTIL_lh5_get( nval_parameters(11),    c0_t )
!    status = UTIL_lh5_get( nval_parameters(12),    c1_t )
!    status = UTIL_lh5_get( nval_parameters(13),    c2_t )
!    log_pres = LOG10( DBLE( pressure ) )
!    wl_t  = 0.1*wl_t !! convert A to nm
!
!    DO ip_nval = 1, np_nval_read
!      !! if no ring added, skip knr2
!      IF( RingAdded(1) == 0 .AND. &
!           INDEX( nval_parameters(ip_nval)%name, "knr2" )>0 ) CYCLE
!      status = UTIL_lh5_disposeDS( nval_parameters(ip_nval ) )
!    ENDDO
!    CALL h5fclose_f( file_id, error)
!    CALL h5close_f(error)
!
!    status = O3T_getLambdaRange()  !! this sets wlT_min & wlT_max
!    IF( status /= OZT_S_SUCCESS ) THEN
!      status = OZT_E_FAILURE
!      ierr = OMI_SMF_setmsg( OZT_E_INPUT, "O3T_getLambdaRange failed.", &
!           FUNCTIONNAME, zero )
!      nval_read = .FALSE.
!      RETURN
!    ENDIF
!
!    IF( wlT_min > MAXVAL( wl_t ) ) THEN
!      status = OZT_E_FAILURE
!      ierr = OMI_SMF_setmsg( OZT_E_INPUT, "wlT_min too large", &
!           FUNCTIONNAME, zero )
!      nval_read = .FALSE.
!      RETURN
!    ELSE IF( wlT_min < MINVAL( wl_t )) THEN
!      iwl_s = 1
!    ELSE
!      i_foo = MAXLOC( wl_t, MASK = wl_t-wlT_min  <= 0.0 )
!      iwl_s = i_foo( 1 )
!    ENDIF
!
!    IF( wlT_max < MINVAL( wl_t ) ) THEN
!      status = OZT_E_FAILURE
!      ierr = OMI_SMF_setmsg( OZT_E_INPUT, "wlT_max too small", &
!           FUNCTIONNAME, zero )
!      nval_read = .FALSE.
!      RETURN
!    ELSE IF( wlT_max > MAXVAL( wl_t )) THEN
!      iwl_e = nwlLUT
!    ELSE
!      i_foo = MINLOC( wl_t, MASK = wl_t-wlT_max  >= 0.0 )
!      iwl_e = i_foo( 1 )
!    ENDIF
!    nwl_com = iwl_e - iwl_s + 1
!
!    WRITE( msg,'(4(A,I5))') "wavelength used in nval calculations:iwl_s =", &
!         iwl_s,",iwl_e =",iwl_e, &
!         ",nwl_com =", nwl_com, ",nwlLUT =", nwlLUT
!    ierr = OMI_SMF_setmsg(OZT_S_SUCCESS, msg, FUNCTIONNAME, four ) 
!
!    n5size = (/ nvzaT*nszaT*nprofT*nwl_com*npresT, &
!         nvzaT*nszaT*nprofT*nwl_com      , &
!         nvzaT*nszaT*nprofT          , &
!         nvzaT*nszaT                , &
!         nvzaT                      /)
!
!    n3size = (/ nprofT*nwl_com*npresT, &
!         nprofT*nwl_com      , &
!         nprofT            /)
!
!    ALLOCATE( wl_com(nwl_com), c0(nwl_com), c1(nwl_com), c2(nwl_com), &  
!         lgi0_n( n5size(1) ), &
!         z1_n(   n5size(1) ), &
!         z2_n(   n5size(1) ), &
!         tr_n(   n5size(1) ), &
!         knb_n(  n5size(1) ), &            
!         sb_n(   n3size(1) ), &            
!         STAT = ierr )
!    ! pick out the unpeturbed layer (layer 0) for all the profiles
!    ! this will be stored in n values for radiance calculation
!    ! HDF used 0 based index, Fortran use 1 based index
!    lgi0_n(:) = RESHAPE( lgi0_t(:,:,:,iwl_s:iwl_e,:), (/ n5size(1) /) )
!    z1_n(:)   = RESHAPE(   z1_t(:,:,:,iwl_s:iwl_e,:), (/ n5size(1) /) )
!    z2_n(:)   = RESHAPE(   z2_t(:,:,:,iwl_s:iwl_e,:), (/ n5size(1) /) )
!    tr_n(:)   = RESHAPE(   tr_t(:,:,:,iwl_s:iwl_e,:), (/ n5size(1) /) )
!    IF( RingAdded(1) == 1 ) THEN
!      knb_n(:) = RESHAPE(knb_t(:,:,:,iwl_s:iwl_e,:), (/ n5size(1) /) )
!    ELSE
!      knb_n(:)  = 0.0
!    ENDIF
!    sb_n(:)   = RESHAPE(   sb_t(    :,iwl_s:iwl_e,:), (/ n3size(1) /) )
!    wl_com(:) =            wl_t(      iwl_s:iwl_e )
!    c0(:)     =            c0_t(      iwl_s:iwl_e )
!    c1(:)     =            c1_t(      iwl_s:iwl_e )
!    c2(:)     =            c2_t(      iwl_s:iwl_e )
!
!    !!OPEN( UNIT=1,FILE='NVAL.txt',FORM='FORMATTED',STATUS='NEW')
!    !!WRITE(1,*) lgi0_t(:,:,:,(/1,4,7,8,9,10/),:)
!    !!WRITE(1,*) lgi0_t(:,:,:,:,:)
!    !!CLOSE(1)
!    nval_read = .TRUE.
!
!    IF( PRESENT( nv ) ) THEN
!      nv%nval_fn   = fn_nval
!      nv%nval_read = nval_read
!
!      IF( RingAdded(1) == 0  ) THEN
!        nv%ringAdded = .FALSE.
!      ELSE IF( RingAdded(1) == 1  ) THEN
!        nv%ringAdded = .TRUE.
!      ELSE
!        status = OZT_E_FAILURE
!        ierr = OMI_SMF_setmsg( OZT_E_INPUT, "unknown RingAdded value", &
!             FUNCTIONNAME, zero )
!        nv%nval_read = .FALSE.
!        RETURN
!      ENDIF
!      nv%nvza      = nvzaT
!      nv%nsza      = nszaT
!      nv%nwl       = nwl_com
!      nv%npres     = npresT
!      nv%nprof     = nprofT
!      nv%n5size(:) = n5size(:)
!      nv%n3size(:) = n3size(:)
!      nv%wl_min    = wlT_min
!      nv%wl_max    = wlT_max
!      ALLOCATE( nv%wl(nv%nwl), nv%c0(nv%nwl), nv%c1(nv%nwl), nv%c2(nv%nwl), &
!           nv%sza(nv%nsza), nv%vza(nv%nvza), nv%pres(nv%npres), &
!           nv%presLog( nv%npres ), nv%lgi0( nv%n5size(1) ), &
!           nv%z1(   nv%n5size(1) ), nv%z2( nv%n5size(1) ),  &
!           nv%tr(   nv%n5size(1) ), nv%knb( nv%n5size(1) ), &  
!           nv%sb(   nv%n3size(1) ), STAT = ierr )
!      nv%presLog(:) = log_pres(:)
!      nv%sza(:)     = sza(:) 
!      nv%vza(:)     = vza(:)
!      nv%pres(:)    = pressure(:)
!      nv%lgi0(:)    = lgi0_n(:) 
!      nv%z1(:)      = z1_n(:)
!      nv%z2(:)      = z2_n(:)
!      nv%tr(:)      = tr_n(:)
!      IF( nv%ringAdded ) THEN
!        nv%knb(:)  = knb_n(:)
!      ELSE
!        nv%knb(:)  = 0.0
!      ENDIF
!      nv%sb(:)      = sb_n(:)
!      nv%wl(:)      = wl_com(:)
!      nv%c0(:)      = c0(:)
!      nv%c1(:)      = c1(:)
!      nv%c2(:)      = c2(:)
!    ENDIF
!
!    DEALLOCATE( lgi0_t, z1_t, z2_t, tr_t, sb_t, knb_t, &
!         wl_t, c0_t, c1_t, c2_t, STAT = ierr )
!
!    status = OZT_S_SUCCESS
!  END FUNCTION O3T_nval_setup_omps

  !...............................................................................


  FUNCTION O3T_getLambdaRange( ) RESULT( status )
    INTEGER (KIND = 4) :: ierr, status
    INTEGER (KIND = 4) :: ilun
    REAL (KIND = 4) :: wl_foo
    CHARACTER ( LEN = 200 ) :: msg
    INTEGER (KIND=4), DIMENSION(2) :: LUNSet = (/ WL_MIN_LUN, WL_MAX_LUN /)
    DO ilun = 1, 2
      status = PGS_PC_GetConfigData( LUNSet(ilun), msg )
      IF( status /= PGS_S_SUCCESS ) THEN
        WRITE( msg,'(A,I9)' ) "get wl failed at LUN = ", LUNSet(ilun)
        ierr = OMI_SMF_setmsg( status, msg, "O3T_getLambdaRange", zero )
        status = OZT_E_FAILURE
        RETURN
      ELSE
        READ( msg, '(F8.4)') wl_foo
        SELECT CASE( ilun )
        CASE( 1 ) 
          wlT_min = wl_foo
        CASE( 2 ) 
          wlT_max = wl_foo
        CASE DEFAULT
          ierr = OMI_SMF_setmsg( status, &
               "unknow ilun", &
               "O3T_getLambdaRange", zero )
          status = OZT_E_FAILURE
          RETURN
        END SELECT
      ENDIF
    ENDDO
    status = OZT_S_SUCCESS
    RETURN
  END FUNCTION O3T_getLambdaRange

  FUNCTION O3T_getLambdaSet( ) RESULT( status )
    INTEGER (KIND = 4) :: ierr, status, i_foo(1)
    INTEGER (KIND = 4) :: ilun
    REAL (KIND = 4) :: wl_foo
    CHARACTER ( LEN = 200 ) :: msg
    INTEGER (KIND=4), DIMENSION(4) :: LUNSet = (/ WL_MIX_LUN, &
         WL_OZ_LUN, WL_REFLLOW_LUN, WL_REFLHIGH_LUN /)

    IF( .NOT. nval_read ) THEN
      WRITE( msg,'(A,I9)' ) "NValue Table must be read in before checking."
      ierr = OMI_SMF_setmsg( status, msg, "O3T_getLambdaSet", zero )
      status = OZT_E_FAILURE
      RETURN
    ENDIF

    DO ilun = 1, 4
      status = PGS_PC_GetConfigData( LUNSet(ilun), msg )
      IF( status /= PGS_S_SUCCESS ) THEN
        WRITE( msg,'(A,I9)' ) "get wl failed at LUN = ", LUNSet(ilun)
        ierr = OMI_SMF_setmsg( status, msg, "O3T_getLambdaSet", zero )
        status = OZT_E_FAILURE
        RETURN
      ELSE
        READ( msg, '(F10.4)') wl_foo
        SELECT CASE( ilun )
        CASE( 1 ) 
          wl_mix = wl_foo
          IF( wl_mix > MAXVAL( wl_com ) .OR. &
               wl_mix < MINVAL( wl_com) ) THEN
            WRITE( msg,'(A,F6.2)') "wl_mix = ", wl_mix, &
                 " outsie LUT range." 
            ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, &
                 "O3T_getLambdaSet", zero )
            status = OZT_E_FAILURE
            RETURN
          ELSE 
            i_foo = MINLOC( ABS(wl_com-wl_mix) )
            iwl_mix = i_foo( 1 )
            WRITE( msg,'(A,F6.2,A,I3,A)') "wl_mix=", wl_com(iwl_mix), &
                 ", iwl_mix = ", iwl_mix, &
                 " used in ozone retrieval."
            ierr = OMI_SMF_setmsg( OZT_S_SUCCESS, msg,  &
                 "O3T_getLambdaSet", zero )
          ENDIF
        CASE( 2 ) 
          wl_oz = wl_foo
          IF( wl_oz > MAXVAL(wl_com) .OR. wl_oz < MINVAL(wl_com) ) THEN
            WRITE( msg,'(A,F6.2)') "wl_oz=", wl_oz, " outsie LUT range."
            ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, &
                 "O3T_getLambdaSet", zero )
            status = OZT_E_FAILURE
            RETURN
          ELSE 
            i_foo = MINLOC( ABS( wl_com-wl_oz ) )
            iwl_ozone = i_foo( 1 )
            WRITE( msg,'(A,F6.2,A,I3,A)') "wl_oz=", wl_com(iwl_ozone), &
                 ", iwl_ozone = ", iwl_ozone, &
                 " used in ozone retrieval."
            ierr = OMI_SMF_setmsg( OZT_S_SUCCESS, msg,  &
                 "O3T_getLambdaSet", zero )
          ENDIF
        CASE( 3 ) 
          wl_refl_l = wl_foo
          IF( wl_refl_l > MAXVAL(wl_com) .OR. &
               wl_refl_l < MINVAL(wl_com) ) THEN
            WRITE( msg,'(A,F6.2)') "wl_refl_l=",wl_refl_l,&
                 " outside LUT range."
            ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, &
                 "O3T_getLambdaSet", zero )
            status = OZT_E_FAILURE
            RETURN
          ELSE 
            i_foo = MINLOC( ABS( wl_com-wl_refl_l ) )
            iwl_refl_l = i_foo( 1 )
            WRITE( msg,'(A,F6.2,A,I3,A)') "wl_refl_l=", &
                 wl_com(iwl_refl_l), ", iwl_refl_l=",iwl_refl_l, &
                 " used in reflectance caclculation."
            ierr = OMI_SMF_setmsg( OZT_S_SUCCESS, msg, &
                 "O3T_getLambdaSet", zero )
          ENDIF
        CASE( 4 ) 
          wl_refl_h = wl_foo
          IF( wl_refl_h > MAXVAL(wl_com) .OR. &
               wl_refl_h < MINVAL(wl_com) ) THEN
            WRITE( msg,'(A,F6.2,A)') "wl_refl_h = ", wl_refl_h, &
                 " outsie LUT range."
            ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, &
                 "O3T_getLambdaSet", zero )
            status = OZT_E_FAILURE
            RETURN
          ELSE 
            i_foo = MINLOC( ABS( wl_com-wl_refl_h ) )
            iwl_refl_h = i_foo( 1 )
            WRITE( msg,'(A,F6.2,A,I3,A)') "wl_refl_h = ", &
                 wl_com(iwl_refl_h), ", iwl_refl_h = ", iwl_refl_h, &
                 " used in reflectance caclculation."
            ierr = OMI_SMF_setmsg( OZT_S_SUCCESS, msg, &
                 "O3T_getLambdaSet", zero )
          ENDIF
        CASE DEFAULT
          ierr = OMI_SMF_setmsg( status, &
               "unknow ilun", &
               "O3T_getLambdaSet", zero )
          status = OZT_E_FAILURE
          RETURN
        END SELECT
      ENDIF
    ENDDO
    status = OZT_S_SUCCESS
    RETURN
  END FUNCTION O3T_getLambdaSet

  FUNCTION O3T_nval_dispose( ) RESULT (status)
    INTEGER (KIND = 4) :: status
    INTEGER  :: ierr

    IF( nval_read ) THEN
      DEALLOCATE( sza, vza, pressure, log_pres, wl_com, c0, c1, c2, &
           lgi0_n, z1_n,  z2_n,  tr_n,  sb_n, knb_n, STAT = ierr )
      nval_read = .FALSE.
    ENDIF
    status = OZT_S_SUCCESS
  END FUNCTION O3T_nval_dispose

END MODULE O3T_nval_class
