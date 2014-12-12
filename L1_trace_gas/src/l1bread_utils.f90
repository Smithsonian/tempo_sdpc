module l1bread_utils
  use tell_module
  use tio_module
  implicit none
  private
  public lookup_swathname, read_l1_radiance_info

contains
  subroutine lookup_swathname (l1b_channel, swathname, errstat)
    implicit none
    character (len=*), intent(in) :: l1b_channel
    character (len=*), intent(inout) :: swathname
    integer, intent(inout) :: errstat

    if (errstat < 0) return

    if (l1b_channel == "UV1") then
      swathname = "band_290_490_nm"
    else if (l1b_channel == "UV2") then
      swathname = "band_540_740_nm"
    else
      call tell_error (tell_internal_error, &
                       "*** read_irradiance_data:  Unsupported value l1b_channel="//l1b_channel, &
                       errstat)
    endif

  end subroutine
  
  subroutine read_l1_radiance_info (l1bfile, l1bchannel, rpt, errstat)

    USE OMSAO_precision_module
    USE OMSAO_variables_module,  ONLY : Radiance_Paras_Type
    use tio_module

    implicit none
    character (len=*), intent (in) :: l1bfile, l1bchannel
    type(Radiance_Paras_Type), INTENT(out) :: rpt
    integer (kind=i4), intent (inout) :: errstat

    type (tiof_object_type) :: tio_l1obj

    if (errstat < 0) return

    rpt%ntimes = 0 ; rpt%nxtrack = 0 ; rpt%nwavel_ccd = 0
    rpt%l1bfilename = l1bfile
    rpt%l1bchannel = l1bchannel

    ! allow error to flow through
    call tiof_open (l1bfile, tio_l1obj, errstat)
    call lookup_swathname (l1bchannel, rpt%swathname, errstat)
    call tiof_inq_group (tio_l1obj, rpt%swathname, errstat)
    call tiof_inq_dimlen (tio_l1obj, "mirror_step", rpt%ntimes, errstat)
    call tiof_inq_dimlen (tio_l1obj, "xtrack", rpt%nxtrack, errstat)
    call tiof_inq_dimlen (tio_l1obj, "spectral_channel", rpt%nwavel_ccd, errstat)
    if (errstat < 0) return
    
    write (*,*) "DEBUG: In read_l1_radiance_info, l1bfile=", &
      trim(l1bfile), ", l1bswath=", trim(rpt%swathname)
    
    call tiof_close (tio_l1obj, errstat)

  end subroutine read_l1_radiance_info
end module l1bread_utils
