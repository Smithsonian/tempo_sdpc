!JBAK 2017-07-25 : for updating OMI operational code
MODULE  m_get_omigeos5

  USE HDF5
  USE OMSAO_precision_module
  USE OMSAO_omidata_module,  ONLY : nxtrack_max, ntimes_max, &
                                    zoom_mode , zoom_p1
  USE OMSAO_variables_module, ONLY :l2_geos5_filename, coadd_uv2,nxbin,nybin,&
                                    offset_line, ncoadd
  USE OMSAO_errstat_module
  USE m_ezspline_interpolation, ONLY: interpol
  USE m_utilities, ONLY: reverse
  IMPLICIT NONE
  !-------------------------------------------------------------------- 
  !GEOS5 Data Variables 
  INTEGER, PARAMETER ::  ngl=47
  TYPE, PUBLIC :: GEOS5_BLock
    REAL (KIND=r4),DIMENSION (:,:),  POINTER :: spres, phis, lat, lon, ptrp
    REAL (KIND=r4),DIMENSION (:,:,:),POINTER :: ts
    REAL (KIND=r4),DIMENSION (:,:,:),POINTER :: ps 
  END TYPE GEOS5_Block  
  TYPE (GEOS5_Block) :: GEOS5 
  !------------------------------------------------------------------------- 
  CONTAINS
   
  SUBROUTINE read_geos5 (nt,nx,nl, pge_error_status)   
    
    IMPLICIT NONE
    !----------------------
    ! INPUT/OUTPUT variables
    !----------------------
    INTEGER, INTENT(IN)        :: nt ! origial alongtrack dimension of UV-2
    INTEGER, INTENT(IN)        :: nx ! origial crosstrack dimension of UV-2 
    INTEGER, INTENT(IN)        :: nl ! alongtrack dimension after sub-extracting and coadding
    INTEGER,INTENT(INOUT)      :: pge_error_status

    !-------------------------
    ! LOCAL variables
    !-------------------------
    ! ncep variables
    INTEGER, PARAMETER :: mrank = 3, ok_flag = 0
    INTEGER (HID_T) :: fid        ! File identifier
    INTEGER (HID_T) :: data_id    ! Dataset identifier
    INTEGER (HID_T) :: data_type
    INTEGER (HSIZE_T), DIMENSION (mrank) :: data_dims
    INTEGER         ::  rank, error ! Error flag
    CHARACTER(maxchlen) :: dname
    ! data variables
    REAL (KIND=r8) :: dtdp
    REAL (KIND=r8),DIMENSION (1:ngl) :: tmp, pmid
    REAL (KIND=r4),DIMENSION(nx, 0:nt-1)     :: lat, lon, spres, ptrp, phis
    REAL (KIND=r4),DIMENSION(1:ngl,nx, 0:nt-1) :: ts, dp
    REAL (KIND=r4),DIMENSION(0:ngl,nx, 0:nt-1) :: ps
    ! Coadding varaibles
    INTEGER :: ix, iy, iix, iiy, j, k ! used in loop
    INTEGER :: nline, nbx, nbin, sline, eline,  nx4, nx1
    INTEGER :: xoff, xoff0
    REAL (KIND=r4):: cc
    LOGICAL :: file_exist
    !----------
    INTEGER, DIMENSION(0:ngl) :: ord ! convert top-down to down-top
    LOGICAL,SAVE :: first = .true.
  
    IF (first) THEN
      ALLOCATE (geos5%spres(nxtrack_max,0:ntimes_max-1))  
      ALLOCATE (geos5%phis(nxtrack_max,0:ntimes_max-1))  
      ALLOCATE (geos5%lat(nxtrack_max,0:ntimes_max-1))  
      ALLOCATE (geos5%lon(nxtrack_max,0:ntimes_max-1))  
      ALLOCATE (geos5%ptrp(nxtrack_max,0:ntimes_max-1))  
      ALLOCATE (geos5%ts(1:ngl,nxtrack_max,0:ntimes_max-1))  
      ALLOCATE (geos5%ps(0:ngl,nxtrack_max,0:ntimes_max-1))  
      first = .false.
    ENDIF

    ! Determine if file exists or not
    INQUIRE (FILE= adjustl(trim(l2_geos5_filename)), EXIST= file_exist)
    IF (.NOT. file_exist) THEN
        WRITE(*, *) 'Warning: no l2 goes5 file!!!', adjustl(trim(l2_geos5_filename))
        pge_error_status = pge_errstat_error ; RETURN
    ENDIF

    ! OPEN 
    CALL h5open_f(pge_error_status)
    IF ( pge_error_status /= ok_flag)  RETURN
  
    CALL h5fopen_f(adjustl(trim(l2_geos5_filename)),H5F_ACC_RDONLY_F,fid, pge_error_status)
    IF ( pge_error_status /= ok_flag)  RETURN

 
    ! Read Geolocation Field
    
    rank = 2 ; data_dims(1:rank)=(/nx, nt/)
    dname = 'lat'
    data_type = H5T_NATIVE_REAL
    CALL h5dopen_f(fid,adjustl(trim(dname)), data_id, error)
    CALL h5dread_f(data_id,data_type,lat(:,0:nt-1), data_dims(1:rank), error)
    CALL h5dclose_f(data_id, error)
   
    dname = 'lon'
    data_type = H5T_NATIVE_REAL
    CALL h5dopen_f(fid,adjustl(trim(dname)), data_id, error)
    CALL h5dread_f(data_id,data_type,lon(:,0:nt-1), data_dims(1:rank), error)
    CALL h5dclose_f(data_id, error)
    
    dname = 'PS' ! pa
    data_type = H5T_NATIVE_REAL
    CALL h5dopen_f(fid,adjustl(trim(dname)), data_id, error)
    CALL h5dread_f(data_id,data_type,spres(:,0:nt-1), data_dims(1:rank), error)
    CALL h5dclose_f(data_id, error)

    dname = 'PHIS' ! m2/s2
    data_type = H5T_NATIVE_REAL
    CALL h5dopen_f(fid,adjustl(trim(dname)), data_id, error)
    CALL h5dread_f(data_id,data_type,phis(:,0:nt-1), data_dims(1:rank), error)
    CALL h5dclose_f(data_id, error)

    dname = 'TROPPB' ! pa
    data_type = H5T_NATIVE_REAL
    CALL h5dopen_f(fid,adjustl(trim(dname)), data_id, error)
    CALL h5dread_f(data_id,data_type,ptrp(:,0:nt-1), data_dims(1:rank), error)
    CALL h5dclose_f(data_id, error)

    dname = 'DELP' !pa
    rank = 3 ; data_dims(1:rank)=(/ngl, nx, nt/)
    data_type = H5T_NATIVE_REAL
    CALL h5dopen_f(fid,adjustl(trim(dname)), data_id, error)
    CALL h5dread_f(data_id,data_type,dp(1:ngl,:,0:nt-1), data_dims(1:rank), error)
    CALL h5dclose_f(data_id, error)

    dname = 'T'
    data_type = H5T_NATIVE_REAL
    CALL h5dopen_f(fid,adjustl(trim(dname)), data_id, error)
    CALL h5dread_f(data_id,data_type,ts(1:ngl,:,0:nt-1), data_dims(1:rank), error)
    CALL h5dclose_f(data_id, error)
    pge_error_status = error
     
    ! CLOSE
    CALL h5fclose_f(fid, pge_error_status)
    IF (pge_error_status /= ok_flag) RETURN 
    CALL h5close_f(pge_error_status)
    IF (pge_error_status /= ok_flag) RETURN 

    !--------------------------------------------------------------------------
    ! Co-adding process just adopting the codes in OMI cloud module
    !--------------------------------------------------------------------------
    IF (zoom_mode) THEN 
     nx1 = nx /2
    ELSE
     nx1 = nx
    ENDIF      
    sline = offset_line ; nline = nl*nybin
    eline = offset_line + nline -1
    nbin = nxbin ; IF (coadd_uv2) nbin = nbin * ncoadd
    nbx = nx1 /nbin
    IF (zoom_mode .AND. MOD(zoom_p1, 2) == 0 ) THEN 
      nbx = nbx -1 ; xoff0 = 1
    ELSE
      xoff0 = 0 
    ENDIF
    DO ix = 1, nx
      DO iy = 0, nt -1
       ps(0,ix,iy)    = 1.0*0.01
       Do k = 1, ngl 
           ps(k,ix,iy) = ps(k-1,ix,iy) + dp(k,ix,iy)*0.01
       ENDDO
      ENDDO
    ENDDO
    GEOS5%lon (1:nbx, 0:nl-1) = 0.0
    GEOS5%lat (1:nbx, 0:nl-1) = 0.0
    GEOS5%ptrp(1:nbx, 0:nl-1) = 0.0
    GEOS5%spres (1:nbx, 0:nl-1) = 0.0
    GEOS5%phis  (1:nbx, 0:nl-1) = 0.0
    GEOS5%ps    (0:ngl, 1:nbx, 0:nl-1) = 0.0
    GEOS5%ts    (1:ngl, 1:nbx, 0:nl-1) = 0.0
    DO ix = 1, nbx                      ! index after co-adding
     DO iy = 0 , nl-1
       iix  = (ix -1)*nbin + 1 + xoff0 ! origianl index befor co-adding
       iiy  = iy*nybin + sline
       cc = 0.0
       Do j = iix, iix + nbin -1
          DO k = iiy, iiy + nybin -1
          GEOS5%ptrp(ix, iy)   =GEOS5%ptrp(ix, iy)   + LOG(ptrp(j,k))
          GEOS5%spres(ix, iy)  =GEOS5%spres(ix, iy)  + LOG(spres(j,k))
          GEOS5%phis(ix, iy)   =GEOS5%phis(ix, iy)   + (phis(j,k))
          GEOS5%ps(0:ngl,ix,iy)=GEOS5%ps(0:ngl,ix,iy)+ LOG(ps(0:ngl,j,k))
          GEOS5%ts(1:ngl,ix,iy)=GEOS5%ts(1:ngl,ix,iy)+ ts(1:ngl,j,k)
          GEOS5%lat(ix, iy)    =GEOS5%lat(ix,iy) + lat(j, k)
          GEOS5%lon(ix, iy)    =GEOS5%lon(ix,iy) + lon(j, k)
          cc = cc + 1
          ENDDO
       ENDDO
       GEOS5%ptrp(ix, iy) = EXP(GEOS5%ptrp(ix, iy)/cc)*0.01 ! convert from Pa to Hpa
       GEOS5%spres(ix, iy) = EXP(GEOS5%spres(ix, iy)/cc)*0.01 ! convert from Pa to Hpa
       GEOS5%phis(ix, iy) = (GEOS5%phis(ix, iy)/cc)/9.8 * 0.001  ! geopohential height (m2/s2) /9.8(m/s2) = > surface height (m) *0.001 => km 
       GEOS5%ps(0:ngl,ix,iy)=EXP(GEOS5%ps(0:ngl, ix, iy)/cc) ! convert from pa to HPa
       GEOS5%ts(1:ngl,ix,iy)=GEOS5%ts(1:ngl,ix,iy)/cc
       GEOS5%lat(ix,iy)     =GEOS5%lat(ix,iy)/cc
       GEOS5%lon(ix,iy)     =GEOS5%lon(ix,iy)/cc       
      ENDDO
    ENDDO
    IF (zoom_mode) THEN 
      xoff = NINT((1.0*zoom_p1 - 1)/nbin)
      GEOS5%ptrp    (1+xoff:nbx+xoff,0:nl-1)=GEOS5%ptrp   (1:nbx,0:nl-1)
      GEOS5%spres   (1+xoff:nbx+xoff,0:nl-1)=GEOS5%spres   (1:nbx,0:nl-1)
      GEOS5%ps(0:ngl,1+xoff:nbx+xoff,0:nl-1)=GEOS5%ps(0:ngl,1:nbx,0:nl-1)
      GEOS5%ts(1:ngl,1+xoff:nbx+xoff,0:nl-1)=GEOS5%ts(1:ngl,1:nbx,0:nl-1)
      GEOS5%lat(1+xoff:nbx+xoff,0:nl-1) = GEOS5%lat(1+xoff:nbx+xoff,0:nl-1)
      GEOS5%lon(1+xoff:nbx+xoff,0:nl-1) = GEOS5%lon(1+xoff:nbx+xoff,0:nl-1)
    ENDIF
    RETURN   
   END SUBROUTINE read_geos5
END MODULE m_get_omigeos5
