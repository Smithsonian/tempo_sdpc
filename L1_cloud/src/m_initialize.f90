!>Routines to read operating settings, and print usage statement if incorrect
module m_initialize

  private ret_usage
  public initialize 

contains

  !>Read operating settings from PCF file, and optionally from command line
  !-------------------------------------------------------------------------
  !
  ! !ROUTINE:  initialize
  ! 
  ! !CALLING SEQUENCE: 
  !
  !        call initialize
  !     
  ! !INPUT PARAMETERS: none (in modules)   
  !> @param errstat error return code, non-zero indicates failure
  !
  ! !OUTPUT PARAMETERS:  
  !
  ! !SEE ALSO:  
  !
  ! !REVISION HISTORY: 
  !
  !> @author  05Jan01  Joiner      original fortran 90
  !> @author  14Mar02  Vasilkov    modified to read OMCTPo.pcf
  !>  23Mar15  O'Sullivan  updating for TEMPO
  !
  !-------------------------------------------------------------------------
  subroutine initialize(errstat)

    use m_vars
    use m_LUN_set
    use m_pgs_include
    use tell_module
    implicit none
    !-------------------------------------------------------------------------
    ! !INPUT PARAMETERS: none (in modules)   
    integer, intent(inout)         :: errstat        ! Error return code:
    !  0   all is well
    !  -1  problem
    !-------------------------------------------------------------------------
    ! Local variables
    integer :: ext, nargs!, i
    integer :: iprt=0   ! verbosity level
    integer :: iarg=0
!    integer :: argc, iargc
    character*255 ::  arg, logmsg
    character*255 ::  pcfpath 
    integer(kind=4), EXTERNAL :: pgs_pc_getnumberoffiles, pgs_pc_getreference
    integer(kind=4), EXTERNAL :: pgs_pc_getuniversalref, pgs_pc_getconfigdata
    integer(kind=4) :: returnstatus, pcf_int
    CHARACTER(LEN=200) :: buf
    !*********************************************************************

    if (errstat /= 0) return


    !===============================
    ! read command line information
    ! and do argument check
    !===============================
!    argc = iargc()
!    do i = 1, 32767
!      iarg = iarg + 1
!      if ( iarg > argc ) go to 111
!      call GetArg ( iArg, argv )
!      if(index(argv,'-h ') > 0) then
!        call ret_usage()
!      else if(index(argv,'-p ') > 0) then
!        if ( iarg+1 > argc ) call ret_usage()
!        iarg = iarg + 1
!        call GetArg ( iArg, argv )
!        read(argv,*,err=500) iprt
!      else if(index(argv,'-nc_swath ') > 0) then
!        if ( iarg+1 > argc ) call ret_usage()
!        iarg = iarg + 1
!        call GetArg ( iArg, argv )
!        read(argv,*,err=500) nc_swathname
!      else if(index(argv,'-nc_only ') > 0) then
!        read_he4 = .false.
!      else if(index(argv,'-noret ') > 0) then
!        noret = .true.
!      endif
!    enddo
!111 continue
    do
      nargs = command_argument_count()
      call get_command_argument (iarg, arg)
      if (len_trim(arg) == 0) exit
      if (trim(arg) == "-p") then
        iarg = iarg + 1
        if (iarg > nargs) then
          call tell_error(tell_usage_error, &
               "initialize: verbosity value not set", errstat)
          call ret_usage()
        else
          call get_command_argument (iarg, arg)
          read(arg,*,err=500) iprt
        endif
      else if (trim(arg) == "-h") then
        call ret_usage()
      else if (trim(arg) == "-nc_in") then
        read_nc = .false.
      else if (trim(arg) == "-nc_out") then
        write_nc = .false.
      else if (trim(arg) == "+he4_in") then
        read_he4 = .true.
      else if (trim(arg) == "+he5_out") then
        write_he5 = .true.
        read_he4 = .true.  ! he4 input required to calculate cloud mask
      else if (trim(arg) == "-tempo") then
        nc_swathname = "band_290_490_nm"
        read_he4 = .false.
        write_he5 = .false.
        have_omi_data = .false.
      else if (trim(arg) == "-noret") then
        noret = .true.
      ! FIXME - temporary flag for proof-of-concept ODL ASCII test
      else if (trim(arg) == "-wrt_odl") then
        wrt_odl = .true.
      endif
      iarg = iarg + 1
    enddo


    ! Set logging level from iprt
    call tell_set_log_level(iprt)

    !***********************************************************************
    !read OMCLDRR.pcf
    !check if PCF exists based on PGS_PC_INFO_FILE environment variable
    !---------------------------------------------------------------------
    status=1
    returnstatus=1
    call tell_log(1,'initialize: checking for pcf file')
    call get_environment_variable('PGS_PC_INFO_FILE',pcfpath)
    inquire(file=pcfpath,exist=ex)
    write(logmsg,"(A,L1)") 'initialize: file status ',ex
    call tell_log(1,logmsg)
    if (ex) then
      version = 1
      status = pgs_pc_getreference ( L1B_LUN, version, filename)
      if (status == 0) then
!        write(logmsg,"(A12,A24,I1)") 'initialize: ', &
!             'get_pc_reference status ',status
!        call tell_log(1,logmsg)
        write(logmsg,"(A13,A)") 'l1b filename ',trim(filename)
        call tell_log(1,logmsg)
      else
        call tell_error(tell_io_read_error, &
             'initialize: unable to read radiance filename from PCF', &
             errstat)
      endif

      returnstatus = pgs_pc_getconfigdata (versionid_lun, buf)
      if (returnstatus == 0) then
        read (buf,*)processing_version
        write(logmsg,"(a,i2)")'processing_version = ',processing_version
        call tell_log(1,logmsg)
      endif

      version = 1
      returnstatus = pgs_pc_getreference(IRR1B_FILE,version,irrad_filename)
      if (returnstatus == 0) then
        write(logmsg,"(A19,A)") 'l1b irrad filename ',trim(irrad_filename)
        call tell_log(1,logmsg)
      else
        call tell_error(tell_io_read_error, &
             'initialize: unable to read irradiance filename from PCF', &
             errstat)
      endif

      version = 1
      returnstatus = pgs_pc_getreference(L2_out,version,flnm_out)
      if (returnstatus == 0) then
        filename_out=trim(flnm_out)
        write(logmsg,"(A16,A)") 'output filename ',trim(filename_out)
        call tell_log(1,logmsg)
      else
        call tell_error(tell_io_read_error, &
             'initialize: unable to read output filename from PCF', &
             errstat)
      endif

      version = 1
      returnstatus = pgs_pc_getconfigdata(using_resid_LUN,buf)
      IF(returnstatus == 0 ) THEN
        read(buf,*) pcf_int
        using_resid = pcf_int == 1
        write(logmsg,"(A34,L1)") 'initialize: setting using_resid = ', &
             using_resid
        call tell_log(1,logmsg)
      endif

      returnstatus = pgs_pc_getconfigdata(write_resid_LUN,buf)
      IF(returnstatus == 0 ) THEN
        read(buf,*) pcf_int
        write_resid = pcf_int == 1
        write(logmsg,"(A34,L1)") 'initialize: setting write_resid = ', &
             write_resid
        call tell_log(1,logmsg)
      endif

      returnstatus = pgs_pc_getconfigdata(do_o3_LUN,buf)
      IF(returnstatus == 0 ) THEN
        read(buf,*) pcf_int
        do_o3 = pcf_int == 1
        write(logmsg,"(A28,L1)") 'initialize: setting do_o3 = ', &
             do_o3
        call tell_log(1,logmsg)
      endif

      returnstatus = pgs_pc_getconfigdata(write_obs_LUN,buf)
      IF(returnstatus == 0 ) THEN
        read(buf,*) pcf_int
        write_obs = pcf_int == 1
        write(logmsg,"(A32,L1)") 'initialize: setting write_obs = ', &
             write_obs
        call tell_log(1,logmsg)
      endif

      returnstatus = pgs_pc_getconfigdata(no_ret_ps_LUN,buf)
      IF(returnstatus == 0 ) THEN
        read(buf,*) pcf_int
        no_ret_ps = pcf_int == 1
        write(logmsg,"(A32,L1)") 'initialize: setting no_ret_ps = ', &
             no_ret_ps
        call tell_log(1,logmsg)
      endif

      returnstatus = pgs_pc_getconfigdata(no_ret_LUN,buf)
      IF(returnstatus == 0) THEN
        if (.NOT. noret) then ! only true if set from command line above
          read(buf,*) pcf_int
          noret = pcf_int == 1
        endif
        write(logmsg,"(A28,L1)") 'initialize: setting noret = ', &
             noret
        call tell_log(1,logmsg)
      endif

      returnstatus = pgs_pc_getconfigdata(transient_chk,buf)
      IF(returnstatus == 0 ) THEN
        read(buf,*) pcf_int
        transient_check = pcf_int == 1
        write(logmsg,"(A38,L1)") 'initialize: setting transient_check = ', &
             transient_check
        call tell_log(1,logmsg)
      endif

      returnstatus = pgs_pc_getconfigdata(test_solar_LUN,buf)
      IF(returnstatus == 0 ) THEN
        read(buf,*) pcf_int
        test_solar = pcf_int == 1
        write(logmsg,"(A33,L1)") 'initialize: setting test_solar = ', &
             test_solar
        call tell_log(1,logmsg)
      endif

      returnstatus = pgs_pc_getconfigdata(add_shift_LUN,buf)
      IF(returnstatus == 0 ) THEN
        read(buf,*) add_shift
        write(logmsg,"(A32,L1)") 'initialize: setting add_shift = ', &
             add_shift
        call tell_log(1,logmsg)
      endif

      returnstatus = pgs_pc_getconfigdata(using_spline_LUN,buf)
      IF(returnstatus == 0 ) THEN
        read(buf,*) pcf_int
        using_spline = pcf_int == 1
        write(logmsg,"(A35,L1)") 'initialize: setting using_spline = ', &
             using_spline
        call tell_log(1,logmsg)
      endif

      returnstatus = pgs_pc_getconfigdata(wmin_LUN,buf)
      IF(returnstatus == 0 ) THEN
        read(buf,*) wmin
        write(logmsg,"(A27,F5.1)") 'initialize: setting wmin = ', &
             wmin
        call tell_log(1,logmsg)
        set_wmin=.true.
      endif
      returnstatus = pgs_pc_getconfigdata(wmax_LUN,buf)
      IF(returnstatus == 0 ) THEN
        read(buf,*) wmax
        write(logmsg,"(A27,F5.1)") 'initialize: setting wmax = ', &
             wmax
        call tell_log(1,logmsg)
        set_wmax=.true.
      endif

      if (wrt_odl) then
        version = 1
        returnstatus = pgs_pc_getreference(mdlist_LUN,version,md_namelist)
        if (returnstatus == 0) then
          write(logmsg,"(A16,A)") 'metadata namelist ',trim(md_namelist)
          call tell_log(1,logmsg)
        else
          call tell_error(tell_io_read_error, &
            'initialize: unable to read metadata namelist filename from PCF', &
            errstat)
        endif
      endif

    else  !ex=.false., PCF does not exist or environment variable not set

      call tell_log(1,'initialize: PCF file not found')
      call tell_error (tell_io_error, &
           "read_cld_dimensions: failed", &
           errstat)
      return

    endif

    !sort out .nc and .he4/5 filenames
    !radiance
    ext = index(filename, '.nc')
    if (ext <= 0) then ! extension is he4
      ext = index(filename, '.he4')
    endif
    filename = trim(filename(1:ext-1))//'.he4'
    filename_in_nc = trim(filename(1:ext-1))//'.nc'
    !irradiance
    ext = index(irrad_filename, '.nc')
    if (ext <= 0) then ! extension is he4
      ext = index(irrad_filename, '.he4')
    endif
    irrad_filename = trim(irrad_filename(1:ext-1))//'.he4'
    irrad_filename_nc = trim(irrad_filename(1:ext-1))//'.nc'
    !output
    ext = index(filename_out, '.nc')
    if (ext <= 0) then ! extension is he5
      ext = index(filename_out, '.he5')
    endif
    filename_out = trim(filename_out(1:ext-1))//'.he5'
    filename_out_nc = trim(filename_out(1:ext-1))//'.nc'
    !***********************************************************************

    return
500 call ret_usage()

  end subroutine initialize

  !****************************************************************************
  !>Print usage statement to the command line
  subroutine ret_usage()

    implicit none

    print *
    print *, &
         'Usage:  L1_cloud [-p iprt] [options] '
    print *
    print *,  'where'
    print *
    print *, '-h          print this message'
    print *
    print *, '-p iprt     verbosity flag (iprt=0: no output, 1: basic output '
    print *, '             >1 detailed output, default=0)'
    print *
    print *, '-noret      do not actually perform retrieval (for testing)'
    print *
    print *, '-nc_in      do not read netCDF input'
    print *, '-nc_out     do not write netCDF output'
    print *, '-he4_in     read HDFEOS4 input'
    print *, '-he5_out    write HDFEOS5 output'
    print *
    print *, '-tempo      expect TEMPO-fortmat netCDF input'
    print *
    print *, "Last Revised: 22 September 2017      (E. O'Sullivan)  "
    print *

    stop 1

  end subroutine ret_usage

end module m_initialize
