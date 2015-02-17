!!****************************************************************************
!!F90
!
!!Description:
!
! MODULE O3T_dndx_class
! 
! read in tables for calculation of forward model quantities:
! dN/dX, and dN/dT.
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
MODULE O3T_dndx_class
  USE O3T_nval_class
  USE OMI_SMF_class

  IMPLICIT NONE
  PUBLIC
  INTEGER (KIND = 4), PARAMETER, PRIVATE :: np_dndx_read = 6
  INTEGER (KIND = 4), PARAMETER, PRIVATE :: zero = 0, one = 1, two = 2
  INTEGER (KIND = 4), PARAMETER, PRIVATE :: three = 3, four = 4, five = 5

  LOGICAL :: dndx_read
  INTEGER (KIND = 4), PUBLIC :: nwl_sub       ! nwl_sub = iwl_c - iwl_s + 1
  INTEGER (KIND = 4), PUBLIC :: nlyrT
  INTEGER (KIND = 4), PUBLIC :: iwl_c

  REAL (KIND = 4), DIMENSION(1), PUBLIC :: lyrDetO3! >0 actual, <= 0 assumed 10%
  !
  ! For all the profiles at all the perturbed layers in the h5 files
  ! lgi0, z1, z2, tr are arrays with a shape as 
  ! (nvzaT, nszaT, nlyrT, nwl_sub', npresT, nprofT ) 
  ! sb is an array with shape ( nlyrT, nwl_sub, npresT, nprofT )
  !
  INTEGER (KIND = 4), DIMENSION(6) :: d6size  
  INTEGER (KIND = 4), DIMENSION(4) :: d4size  
  REAL (KIND = 4), DIMENSION(:), ALLOCATABLE:: lgi0_d, z1_d, &
                                               z2_d, tr_d, sb_d
  PUBLIC :: O3T_dndx_setup
!  PUBLIC :: O3T_dndx_setup_omps
  PUBLIC :: O3T_dndx_dispose

  CONTAINS  
!!
! reads parameters from HDF5 file directly into ALLOCATABLE arrays for use 
! by the algorithm. Parameters have dimensionality of 
! (nvzaT, nszaT, nprofT, nwlLUT, npresT, nlyrT ), for the the partial 
! derivative (dndx-value) tables corresponding to (viewing zenith, 
! solar zenith, profile, wavelength, pressure, and layer). The layer
! dimension of 12 is for the unperturbed case followed by 11 umkehr 
! layers 0-9 and 10, 11, and 12 combined. The atmospheric backscatter to 
! ground terms, sb, have no angular dependence, so the dimensionality is 
! reduced.  
!!

    FUNCTION O3T_dndx_setup( wl_cutoff ) RESULT( status )
      USE HDF5
      USE O3T_LUT_fs
      USE PGS_PC_class
      USE OMI_LUN_set

      REAL(KIND=4), INTENT(IN), OPTIONAL :: wl_cutoff 
      CHARACTER ( LEN = 200 ) :: fn_dndx

      INTEGER (HID_T) :: file_id
      INTEGER (KIND = 4) :: status
      INTEGER :: ip_lut, iprof, ilyr, error, ierr, i_foo(1) ! is, ie, 
      INTEGER :: iwl, ipres, isza, ivza
      INTEGER (KIND=4) :: nvzaP, nszaP, nprofP, nwlLUTP, npresP 
      !INTEGER :: iwl_max
      !INTEGER, DIMENSION(6) :: data_dims
      TYPE (DSH5_T), DIMENSION(np_dndx_read) :: lut_parameters
      INTEGER, DIMENSION(np_dndx_read) :: rank_dndx = (/ 6,6,6,6,4,1 /) 
      INTEGER :: count6d
      INTEGER :: count4d
      REAL (KIND = 4), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: &
                                                    lgi0_t, z1_t, z2_t, tr_t
      REAL (KIND = 4), DIMENSION(:,:,:,:), ALLOCATABLE :: sb_t
      REAL (KIND = 4), DIMENSION(:), ALLOCATABLE :: wl_t
      CHARACTER ( LEN = 200 ) :: msg
      !INTEGER (HID_T)  :: datatype
      INTEGER :: version

      IF( .NOT. nval_read ) THEN
         WRITE( msg,'(A)' ) "nval table must be read in first before dndx tab"
         ierr = OMI_SMF_setmsg( status, msg, "O3T_dndx_setup", zero )
         status = OZT_E_FAILURE
         RETURN
      ENDIF

      !! get the dndx table filename and path
      version = 1
      status = PGS_PC_GetReference( OMTO3_DNDX_LUN, version, fn_dndx )
      IF( status /= PGS_S_SUCCESS ) THEN
         WRITE( msg,'(A)' ) "get LUT name failed, PGE aborting, exit code = 1."
         ierr = OMI_SMF_setmsg( status, msg, "O3T_dndx_setup", zero )
         status = OZT_E_FAILURE
         RETURN
      ENDIF

      lut_parameters  = (/ ds_lgi0, ds_z1, ds_z2, ds_tr, ds_sb, ds_wl /)

      ! Initialize FORTRAN interface. !
      CALL h5open_f(error)
 
      CALL h5fopen_f( fn_dndx, H5F_ACC_RDONLY_F, file_id, error )
      IF( error < zero ) THEN
         status = OZT_E_FAILURE
         WRITE( msg,* ) "h5fopen_f failed on file: ", fn_dndx
         ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), "O3T_dndx_setup", zero )
         dndx_read = .FALSE.
         RETURN
      ENDIF

      ds_dx%datatype = H5T_NATIVE_REAL
      ds_dx%dims(1)  = 1
      CALL h5tget_class_f(ds_dx%datatype, ds_dx%class, error )
      status = UTIL_lh5_selectDS( file_id, ds_dx )
      IF( status == OZT_S_SUCCESS ) THEN
         status = UTIL_lh5_get( ds_dx, lyrDetO3 )
         status = UTIL_lh5_disposeDS( ds_dx )
      ELSE
         lyrDetO3(:) = 0.0
      ENDIF

      DO ip_lut = 1, np_dndx_read
        lut_parameters(ip_lut)%datatype = H5T_NATIVE_REAL
        CALL h5tget_class_f(lut_parameters(ip_lut)%datatype, &
                            lut_parameters(ip_lut)%class, error )
        IF( error < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "h5tget_class_f failed on dataset: ", &
                          lut_parameters(ip_lut)%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                 "O3T_dndx_setup", zero )
           RETURN
        ENDIF     
        lut_parameters(ip_lut)%rank = rank_dndx(ip_lut)
        status = UTIL_lh5_selectDS( file_id, lut_parameters(ip_lut) )
        IF( status .NE. OZT_S_SUCCESS ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "UTIL_lh5_selectDS failed for ", &
                          lut_parameters(ip_lut)%name, "in file ", fn_dndx
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                  "O3T_dndx_setup", zero )
           dndx_read = .FALSE.
           RETURN
          ENDIF
      ENDDO

      nvzaP  = lut_parameters(1)%dims(1)
      nszaP  = lut_parameters(1)%dims(2)
      nprofP = lut_parameters(1)%dims(3)
      nwlLUTP= lut_parameters(1)%dims(4)
      npresP = lut_parameters(1)%dims(5)
      nlyrT  = lut_parameters(1)%dims(6)

      IF( nvzaP /= nvzaT .OR. nszaP /= nszaT .OR. npresP /= npresT &
          .OR. nprofP /= nprofT .OR. nwlLUTP /= nwlLUT ) THEN
         WRITE( msg,* ) "nval and dndx tables on different size grids."
         ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                "O3T_dndx_setup", zero )
         dndx_read = .FALSE.
         status = OZT_E_FAILURE
         RETURN
      ENDIF

      ! storeage space for the all the layers, all profiles, all 
      ! wavelengths, all pressures, all solar and view angles
      ALLOCATE( lgi0_t( nvzaT, nszaT, nprofT, nwlLUT, npresT, nlyrT ), &
                z1_t(   nvzaT, nszaT, nprofT, nwlLUT, npresT, nlyrT ), &
                z2_t(   nvzaT, nszaT, nprofT, nwlLUT, npresT, nlyrT ), &
                tr_t(   nvzaT, nszaT, nprofT, nwlLUT, npresT, nlyrT ), &
                sb_t(                 nprofT, nwlLUT, npresT, nlyrT ), &
                wl_t(nwlLUT), STAT = ierr )
      IF( ierr /= zero ) THEN
         status = OZT_E_FAILURE
         ierr = OMI_SMF_setmsg( OZT_E_MEM_ALLOC, ".. allocation failure", &
                               "O3T_dndx_setup", zero )
         dndx_read = .FALSE.
         RETURN
      ENDIF

      !! read the data stored in the DNDX....h5 file
      status = UTIL_lh5_get( lut_parameters(1),   lgi0_t )
      status = UTIL_lh5_get( lut_parameters(2),     z1_t )
      status = UTIL_lh5_get( lut_parameters(3),     z2_t )
      status = UTIL_lh5_get( lut_parameters(4),     tr_t )
      status = UTIL_lh5_get( lut_parameters(5),     sb_t )
      status = UTIL_lh5_get( lut_parameters(6),     wl_t )
      wl_t  = 0.1*wl_t !! convert A to nm

      DO ip_lut = 1, np_dndx_read
        status = UTIL_lh5_disposeDS( lut_parameters(ip_lut ) )
      ENDDO
      CALL h5fclose_f( file_id, error)
      CALL h5close_f(error)

      IF( PRESENT( wl_cutoff ) ) THEN
         IF( wl_cutoff > wlT_max .OR. wl_cutoff > MAXVAL( wl_t ) ) THEN
            iwl_c = iwl_e
         ELSE
            i_foo = MAXLOC( wl_t, MASK = wl_t-wl_cutoff <= 0.0 )
            iwl_c = i_foo( 1 )
         ENDIF
      ELSE
         iwl_c = iwl_e
      ENDIF

      nwl_sub = iwl_c - iwl_s + 1 
      WRITE( msg,'(4(A,I3))') "wavelength used in dndx calculations:iwl_s =", &
                             iwl_s,",iwl_c =",iwl_c, &
                             ",nwl_sub =",nwl_sub, ",nwl_com =", nwl_com 
      ierr = OMI_SMF_setmsg(OZT_S_SUCCESS, msg, "O3T_dndx_setup", four ) 

      d6size = (/ nvzaT*nszaT*nlyrT*nwl_sub*npresT*nprofT, &
                  nvzaT*nszaT*nlyrT*nwl_sub*npresT      , &
                  nvzaT*nszaT*nlyrT*nwl_sub            , &
                  nvzaT*nszaT*nlyrT                  , &
                  nvzaT*nszaT                       , &
                  nvzaT                             /)
      d4size = (/ nlyrT*nwl_sub*npresT*nprofT, &
                  nlyrT*nwl_sub*npresT      , &
                  nlyrT*nwl_sub            , &
                  nlyrT                    /)

      ALLOCATE( lgi0_d( d6size(1) ), &
                z1_d(   d6size(1) ), &
                z2_d(   d6size(1) ), &
                tr_d(   d6size(1) ), &
                sb_d(   d4size(1) ), &            
                STAT = ierr )

    ! rearrange the lgi0_t, etc... for dndx calculation
      count6d = 1
      count4d = 1
      DO iprof = 1, nprofT
      DO ipres = 1, npresT
      DO iwl   = iwl_s, iwl_c   !! iwl refers to wl_t array, NOT wl_com
        DO ilyr = 1, nlyrT
          DO isza = 1, nszaT
          DO ivza = 1, nvzaT
            lgi0_d(count6d) = lgi0_t( ivza,isza,iprof,iwl,ipres,ilyr )
              z1_d(count6d) =   z1_t( ivza,isza,iprof,iwl,ipres,ilyr )
              z2_d(count6d) =   z2_t( ivza,isza,iprof,iwl,ipres,ilyr )
              tr_d(count6d) =   tr_t( ivza,isza,iprof,iwl,ipres,ilyr )
              count6d = count6d + 1
          ENDDO
          ENDDO
          sb_d(count4d) = sb_t( iprof,iwl,ipres,ilyr )
          count4d = count4d + 1
        ENDDO
      ENDDO
      ENDDO
      ENDDO

      DEALLOCATE( lgi0_t, z1_t, z2_t, tr_t, sb_t, wl_t, STAT = ierr )
      dndx_read = .TRUE.

      status = OZT_S_SUCCESS
    END FUNCTION O3T_dndx_setup

    !...........................................................................

!    FUNCTION O3T_dndx_setup_omps( dndx_lun, wl_cutoff ) RESULT( status )
!      USE HDF5
!      USE O3T_LUT_fs
!      USE PGS_PC_class
!
!      INTEGER (KIND = 4), INTENT(IN) :: dndx_lun
!      REAL(KIND=4), INTENT(IN), OPTIONAL :: wl_cutoff 
!      CHARACTER ( LEN = 200 ) :: fn_dndx
!
!      INTEGER (HID_T) :: file_id
!      INTEGER (KIND = 4) :: status
!      INTEGER :: ip_lut, is, ie, iprof, ilyr, error, ierr, i_foo(1)
!      INTEGER :: iwl, ipres, isza, ivza
!      INTEGER (KIND=4) :: nvzaP, nszaP, nprofP, nwlLUTP, npresP 
!      INTEGER :: iwl_max
!      INTEGER, DIMENSION(6) :: data_dims
!      TYPE (DSH5_T), DIMENSION(np_dndx_read) :: lut_parameters
!      INTEGER, DIMENSION(np_dndx_read) :: rank_dndx = (/ 6,6,6,6,4,1 /) 
!      INTEGER :: count6d
!      INTEGER :: count4d
!      REAL (KIND = 4), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: &
!                                                    lgi0_t, z1_t, z2_t, tr_t
!      REAL (KIND = 4), DIMENSION(:,:,:,:), ALLOCATABLE :: sb_t
!      REAL (KIND = 4), DIMENSION(:), ALLOCATABLE :: wl_t
!      CHARACTER ( LEN = 200 ) :: msg
!      INTEGER (HID_T)  :: datatype
!      INTEGER :: version
!
!      IF( .NOT. nval_read ) THEN
!         WRITE( msg,'(A)' ) "nval table must be read in first before dndx tab"
!         ierr = OMI_SMF_setmsg( status, msg, "O3T_dndx_setup", zero )
!         status = OZT_E_FAILURE
!         RETURN
!      ENDIF
!
!      !! get the dndx table filename and path
!      version = 1
!!     status = PGS_PC_GetReference( OMTO3_DNDX_LUN, version, fn_dndx )
!      status = PGS_PC_GetReference( dndx_lun, version, fn_dndx )
!      IF( status /= PGS_S_SUCCESS ) THEN
!         WRITE( msg,'(A)' ) "get LUT name failed, PGE aborting, exit code = 1."
!         ierr = OMI_SMF_setmsg( status, msg, "O3T_dndx_setup", zero )
!         status = OZT_E_FAILURE
!         RETURN
!      ELSE
!         print *, "[O3T_dndx_setup_omps]:DNDX table =", trim(fn_dndx)
!      ENDIF
!
!      lut_parameters  = (/ ds_lgi0, ds_z1, ds_z2, ds_tr, ds_sb, ds_wl /)
!
!      ! Initialize FORTRAN interface. !
!      CALL h5open_f(error)
! 
!      CALL h5fopen_f( fn_dndx, H5F_ACC_RDONLY_F, file_id, error )
!      IF( error < zero ) THEN
!         status = OZT_E_FAILURE
!         WRITE( msg,* ) "h5fopen_f failed on file: ", fn_dndx
!         ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), "O3T_dndx_setup", zero )
!         dndx_read = .FALSE.
!         RETURN
!      ENDIF
!
!      ds_dx%datatype = H5T_NATIVE_REAL
!      ds_dx%dims(1)  = 1
!      CALL h5tget_class_f(ds_dx%datatype, ds_dx%class, error )
!      status = UTIL_lh5_selectDS( file_id, ds_dx )
!      IF( status == OZT_S_SUCCESS ) THEN
!         status = UTIL_lh5_get( ds_dx, lyrDetO3 )
!         status = UTIL_lh5_disposeDS( ds_dx )
!      ELSE
!         lyrDetO3(:) = 0.0
!      ENDIF
!
!      DO ip_lut = 1, np_dndx_read
!        lut_parameters(ip_lut)%datatype = H5T_NATIVE_REAL
!        CALL h5tget_class_f(lut_parameters(ip_lut)%datatype, &
!                            lut_parameters(ip_lut)%class, error )
!        IF( error < zero ) THEN
!           status = OZT_E_FAILURE
!           WRITE( msg,* ) "h5tget_class_f failed on dataset: ", &
!                          lut_parameters(ip_lut)%name
!           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
!                                 "O3T_dndx_setup", zero )
!           RETURN
!        ENDIF     
!        lut_parameters(ip_lut)%rank = rank_dndx(ip_lut)
!        status = UTIL_lh5_selectDS( file_id, lut_parameters(ip_lut) )
!        IF( status .NE. OZT_S_SUCCESS ) THEN
!           status = OZT_E_FAILURE
!           WRITE( msg,* ) "UTIL_lh5_selectDS failed for ", &
!                          lut_parameters(ip_lut)%name, "in file ", fn_dndx
!           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
!                                  "O3T_dndx_setup", zero )
!           dndx_read = .FALSE.
!           RETURN
!          ENDIF
!      ENDDO
!
!      nvzaP  = lut_parameters(1)%dims(1)
!      nszaP  = lut_parameters(1)%dims(2)
!      nprofP = lut_parameters(1)%dims(3)
!      nwlLUTP= lut_parameters(1)%dims(4)
!      npresP = lut_parameters(1)%dims(5)
!      nlyrT  = lut_parameters(1)%dims(6)
!
!      IF( nvzaP /= nvzaT .OR. nszaP /= nszaT .OR. npresP /= npresT &
!          .OR. nprofP /= nprofT .OR. nwlLUTP /= nwlLUT ) THEN
!         WRITE( msg,* ) "nval and dndx tables on different size grids."
!         ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
!                                "O3T_dndx_setup", zero )
!         dndx_read = .FALSE.
!         status = OZT_E_FAILURE
!         RETURN
!      ENDIF
!
!      ! storeage space for the all the layers, all profiles, all 
!      ! wavelengths, all pressures, all solar and view angles
!      ALLOCATE( lgi0_t( nvzaT, nszaT, nprofT, nwlLUT, npresT, nlyrT ), &
!                z1_t(   nvzaT, nszaT, nprofT, nwlLUT, npresT, nlyrT ), &
!                z2_t(   nvzaT, nszaT, nprofT, nwlLUT, npresT, nlyrT ), &
!                tr_t(   nvzaT, nszaT, nprofT, nwlLUT, npresT, nlyrT ), &
!                sb_t(                 nprofT, nwlLUT, npresT, nlyrT ), &
!                wl_t(nwlLUT), STAT = ierr )
!      IF( ierr /= zero ) THEN
!         status = OZT_E_FAILURE
!         ierr = OMI_SMF_setmsg( OZT_E_MEM_ALLOC, ".. allocation failure", &
!                               "O3T_dndx_setup", zero )
!         dndx_read = .FALSE.
!         RETURN
!      ENDIF
!
!      !! read the data stored in the DNDX....h5 file
!      status = UTIL_lh5_get( lut_parameters(1),   lgi0_t )
!      status = UTIL_lh5_get( lut_parameters(2),     z1_t )
!      status = UTIL_lh5_get( lut_parameters(3),     z2_t )
!      status = UTIL_lh5_get( lut_parameters(4),     tr_t )
!      status = UTIL_lh5_get( lut_parameters(5),     sb_t )
!      status = UTIL_lh5_get( lut_parameters(6),     wl_t )
!      wl_t  = 0.1*wl_t !! convert A to nm
!
!      DO ip_lut = 1, np_dndx_read
!        status = UTIL_lh5_disposeDS( lut_parameters(ip_lut ) )
!      ENDDO
!      CALL h5fclose_f( file_id, error)
!      CALL h5close_f(error)
!
!      IF( PRESENT( wl_cutoff ) ) THEN
!         IF( wl_cutoff > wlT_max .OR. wl_cutoff > MAXVAL( wl_t ) ) THEN
!            iwl_c = iwl_e
!         ELSE
!            i_foo = MAXLOC( wl_t, MASK = wl_t-wl_cutoff <= 0.0 )
!            iwl_c = i_foo( 1 )
!         ENDIF
!      ELSE
!         iwl_c = iwl_e
!      ENDIF
!
!      nwl_sub = iwl_c - iwl_s + 1 
!      WRITE( msg,'(4(A,I3))') "wavelength used in dndx calculations:iwl_s =", &
!                             iwl_s,",iwl_c =",iwl_c, &
!                             ",nwl_sub =",nwl_sub, ",nwl_com =", nwl_com 
!      ierr = OMI_SMF_setmsg(OZT_S_SUCCESS, msg, "O3T_dndx_setup", four ) 
!
!      d6size = (/ nvzaT*nszaT*nlyrT*nwl_sub*npresT*nprofT, &
!                  nvzaT*nszaT*nlyrT*nwl_sub*npresT      , &
!                  nvzaT*nszaT*nlyrT*nwl_sub            , &
!                  nvzaT*nszaT*nlyrT                  , &
!                  nvzaT*nszaT                       , &
!                  nvzaT                             /)
!      d4size = (/ nlyrT*nwl_sub*npresT*nprofT, &
!                  nlyrT*nwl_sub*npresT      , &
!                  nlyrT*nwl_sub            , &
!                  nlyrT                    /)
!
!      ALLOCATE( lgi0_d( d6size(1) ), &
!                z1_d(   d6size(1) ), &
!                z2_d(   d6size(1) ), &
!                tr_d(   d6size(1) ), &
!                sb_d(   d4size(1) ), &            
!                STAT = ierr )
!
!    ! rearrange the lgi0_t, etc... for dndx calculation
!      count6d = 1
!      count4d = 1
!      DO iprof = 1, nprofT
!      DO ipres = 1, npresT
!      DO iwl   = iwl_s, iwl_c   !! iwl refers to wl_t array, NOT wl_com
!        DO ilyr = 1, nlyrT
!          DO isza = 1, nszaT
!          DO ivza = 1, nvzaT
!            lgi0_d(count6d) = lgi0_t( ivza,isza,iprof,iwl,ipres,ilyr )
!              z1_d(count6d) =   z1_t( ivza,isza,iprof,iwl,ipres,ilyr )
!              z2_d(count6d) =   z2_t( ivza,isza,iprof,iwl,ipres,ilyr )
!              tr_d(count6d) =   tr_t( ivza,isza,iprof,iwl,ipres,ilyr )
!              count6d = count6d + 1
!          ENDDO
!          ENDDO
!          sb_d(count4d) = sb_t( iprof,iwl,ipres,ilyr )
!          count4d = count4d + 1
!        ENDDO
!      ENDDO
!      ENDDO
!      ENDDO
!
!      DEALLOCATE( lgi0_t, z1_t, z2_t, tr_t, sb_t, wl_t, STAT = ierr )
!      dndx_read = .TRUE.
!
!      status = OZT_S_SUCCESS
!    END FUNCTION O3T_dndx_setup_omps

    !...........................................................................

    FUNCTION O3T_dndx_dispose( ) RESULT (status)
      INTEGER (KIND = 4) :: status
      INTEGER  :: ierr

      IF( dndx_read ) THEN
         DEALLOCATE( lgi0_d, z1_d,  z2_d,  tr_d,  sb_d, STAT = ierr )
         dndx_read = .FALSE.
      ENDIF
      status = OZT_S_SUCCESS

    END FUNCTION O3T_dndx_dispose

END MODULE O3T_dndx_class
