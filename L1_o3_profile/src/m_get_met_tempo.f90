MODULE  m_get_met_tempo
  USE OMSAO_precision_module, ONLY: i4,r4,sp, dp
  USE OMSAO_parameters_module, ONLY: maxchlen
  USE OMSAO_variables_module, ONLY: l2_met_filenames, the_lon, the_lat, the_surfalt
  USE OMI_LUN_SET, ONLY:num_met_luns
  TYPE :: met_type
     LOGICAL :: do_debug=.false.
     INTEGER :: np
     REAL (kind=r4) ::ptrop, psurf, z0
     REAL (kind=dp), DIMENSION(:), allocatable :: ts, ps
  END TYPE 
  TYPE (met_type) :: thismet
  CONTAINS
  SUBROUTINE get_met_tempo (errstat)
  USE, intrinsic :: iso_c_binding, only: c_ptr, c_char,c_null_char, c_null_ptr, c_associated
  USE met_module
  IMPLICIT NONE
  !===========================
  ! input/output 
  !===========================
  INTEGER, INTENT(OUT) :: errstat
  !===========================
  ! local variables
  !===========================
  INTEGER :: np,imet, ilay, status
  INTEGER, PARAMETER :: met_flags = 7, nlecm = 30
  real(kind=dp), dimension(:), allocatable :: pprof, tprof
  real(kind=sp), dimension(:), allocatable :: tprof_r4
  real(kind=sp) :: lon_f, lat_f
  real(kind=dp) :: local_srf
  real(kind=dp), dimension(0:nlecm), parameter :: clim = & 
       (/288.2, 287.5,286.1, 284.7, 283.2, 281.8, &
        278.7, 275.5, 272.2, 268.6, 264.8, 260.8, & 
       256.6, 251.9, 246.9, 241.4,  235.4, 228.6, & 
       220.5,216.5, 216.7, 216.7,   216.7, 217.2,   &
       220.5, 223.1, 227.7, 232.4,  239.3,249.5, 257.9/)
  type (synth_met_type), SAVE :: smt
  type (c_ptr), SAVE :: met = c_null_ptr
  logical, SAVE :: have_synthetic_met_data
  logical, SAVE :: first = .true.
  
  ! Initialize dataset
  IF (first) THEN 
    if (0 /= index (l2_met_filenames(1), '.nc', .true.)) then
      have_synthetic_met_data = .true.
      call open_synth_met_data (smt, trim(l2_met_filenames(1)),errstat)
      print * , l2_met_filenames(1)
      if (errstat /= 0) then
        call tell_error (tell_runtime_error, "error opening synthetic met data", errstat)
        STOP 1
      endif
    else
      have_synthetic_met_data = .false.
      met = met_list_new (met_flags)
      if (.not. c_associated( met)) then
        call tell_error (tell_runtime_error, "met_list_new returned NULL", errstat)
        return
      endif
      do imet=1, num_met_luns
        if (0 /= index(l2_met_filenames(imet), 'grib2', .true.)) then
         status = met_list_add_file (met,trim(l2_met_filenames(imet))//c_null_char)
          if (status /= 0) then
            call met_list_free (met)
            call tell_error (tell_runtime_error, & 
            "reading:"//trim(l2_met_filenames(imet)), errstat)
            return
          endif
       endif
      enddo
    endif
    first = .false.
   endif

   
   ! define user vertical grids 
   np = nlecm 
   thismet%np = np
   if (allocated(thismet%ts)) deallocate(thismet%ts, thismet%ps)
   allocate(thismet%ts(0:thismet%np),thismet%ps(0:thismet%np))
   allocate(tprof(0:np),pprof(0:np), tprof_r4(0:np))
   pprof(0:np) = (/1013.25, 1000., 975., 950., 925., 900., 850., &
        800., 750., 700., 650., 600., 550., 500., 450., 400., 350., 300.,250.,200., 150., 100., 70., 50., &
        30., 20., 10., 7., 5., 3., 2./)    

   !Interpolation dataset into current pixel
   lon_f = real(the_lon, kind=sp) ; lat_f = real(the_lat, kind=sp)
   IF (have_synthetic_met_data) then
    CALL read_synth_met_data(smt, lon_f, lat_f, thismet%ptrop, errstat, &
                             pprof = pprof, tprof=tprof)
   ELSE
    call met_list_interp_f (met, lon_f, lat_f, errstat, &
                            thismet%psurf, ptrop=thismet%ptrop, &  
                            isobars=real(pprof, kind=sp),temp_on_isobar=tprof_r4)
   ENDIF
   thismet%ps(0:thismet%np)=pprof
   thismet%ts(0:thismet%np)=tprof
   deallocate(tprof, pprof, tprof_r4)
   IF (errstat /= 0) THEN 
      WRITE(*,*) 'Errors in get_met_tempo' ; STOP 1
   ENDIF
   
   ! heanding the pressure grids of out of range
   do ilay = 0, np
     if (thismet%ts(ilay) < 0 .or. thismet%ts(ilay) > 1000) THEN
         thismet%ts(ilay) = clim(ilay)
     endif
   enddo

   ! surface pressure 

   local_srf = REAL(the_surfalt, kind=dp)
   local_srf = 1013.0_dp * (10.0_dp ** (local_srf / 1000.0_dp / (-16.0_dp)))
   thismet%psurf = real(local_srf, kind=sp)
   thismet%z0 = real(the_surfalt, kind=sp)

   ! Debuging
   IF (thismet%do_debug) THEN 
    write(*,'(f8.1)')thismet%psurf
    write(*,'(f8.1)')thismet%ptrop
    write(*,'(f8.1)')thismet%ts
   ENDIF
   RETURN
  END SUBROUTINE get_met_tempo
END MODULE m_get_met_tempo
