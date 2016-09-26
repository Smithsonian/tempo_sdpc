!>Routines for setting quality flags
module m_lambda_qual

  private
  public bad_rad_lambda, bad_irrad_lambda

contains

  !>Set radiance quality flags
  subroutine bad_rad_lambda(ip, iLine, errstat)

    use m_vars, ONLY: quality_flagL, qc, wmin, wmax, transient_check
    use m_cloud_pres_mod, ONLY: f1p, w1p
    use m_find

    implicit none

    integer, intent(in) :: ip, iLine, errstat
    integer :: iw                    
    integer :: iw_start, iw_end
    logical :: pxl_error, pxl_warning

    if (errstat /= 0) return

    !find starting and ending wavelengths
    !====================================
    if (count(w1p(:)<= wmax .and. w1p(:) >= wmin) > 0) then
      iw_end=maxval(find2(w1p(:) <= wmax .and. w1p(:) >= wmin, &
           count(w1p(:) <= wmax .and. w1p(:) >= wmin)))-1
    else
      iw_end=-1
    endif
    if (count(w1p(:) >= wmin) > 1) then
      iw_start=minval(find2(w1p(:) >= wmin,count(w1p(:) >= wmin)))-1
    else
      iw_start=-1
    endif

    ! check the image quality flags
    !================================
    if (iw_start >= 0 .and. iw_end > 0) then
      do iw=iw_start, iw_end
        pxl_error = (BTEST(quality_flagL(iw+1,ip+1),0) .or. & ! missing
             BTEST(quality_flagL(iw+1,ip+1),1) .or. & ! bad
             BTEST(quality_flagL(iw+1,ip+1),2) .or. & ! error
             f1p(iw+1) <= 0.)
        if (transient_check) pxl_error = pxl_error .or. BTEST(quality_flagL(iw+1,ip+1),3)
        pxl_warning = ( &
             BTEST(quality_flagL(iw+1,ip+1),4) .or. & ! RTS warning
             BTEST(quality_flagL(iw+1,ip+1),5) .or. & ! saturation warning
             BTEST(quality_flagL(iw+1,ip+1),6) .or. & ! noise warning
             BTEST(quality_flagL(iw+1,ip+1),7) .or. & ! dark current warning
             BTEST(quality_flagL(iw+1,ip+1),8) .or. & ! offset warning
             BTEST(quality_flagL(iw+1,ip+1),9) .or. & ! smear warning
             BTEST(quality_flagL(iw+1,ip+1),10)) ! stray light warning
        if (.not. transient_check) pxl_warning = pxl_warning .or. &
             BTEST(quality_flagL(iw+1,ip+1),3)

        if(pxl_error) then
          f1p(iw+1)=0.
          qc(ip,iLine)=IBSET(qc(ip,iLine),9)
        endif ! if pxl_error
        if(pxl_warning) qc(ip,iLine)=IBSET(qc(ip,iLine),10)

      enddo ! loop over 
    else
      qc(ip,iLine)=IBSET(qc(ip,iLine),9)
    endif

  end subroutine bad_rad_lambda

  !>Set irradiance quality flags
  subroutine bad_irrad_lambda(nXtrack, errstat)

    use m_vars, ONLY: irr_quality_flagL, qc, nsolwave, ws, fs, wmin, wmax, &
         transient_check 
    use m_find
    use tell_module

    implicit none

    integer, intent(in) :: nXtrack, errstat
    integer :: iw, ip
    integer :: iw_start, iw_end              !, iw_start2, iw_end2
    logical, dimension(nsolwave,nXtrack) :: pxl_error
    logical :: pxl_warning
    character (len=128) :: logmsg


    if (errstat /= 0) return

    do ip=0,nXtrack-1

      !find starting and ending wavelengths
      !====================================
      iw_end=maxval(find2(ws(:,ip) <= wmax,count(ws(:,ip) <= wmax)))-1
      iw_start=minval(find2(ws(:,ip) >= wmin,count(ws(:,ip) >= wmin)))-1

      ! check the image quality flags
      !==============================
      do iw=iw_start,iw_end
        pxl_error(iw+1,ip+1) = (BTEST(irr_quality_flagL(iw+1,ip+1),0) .or. &
             BTEST(irr_quality_flagL(iw+1,ip+1),1) .or. &
             BTEST(irr_quality_flagL(iw+1,ip+1),2) &
             .or. fs(iw,ip) < 1e12 .or. fs(iw,ip) > 1e17 &
             )
        if (transient_check) pxl_error(iw+1,ip+1) = pxl_error(iw+1,ip+1) .or. &
             BTEST(irr_quality_flagL(iw+1,ip+1),3)
        pxl_warning = ( &
             BTEST(irr_quality_flagL(iw+1,ip+1),4) .or. &
             BTEST(irr_quality_flagL(iw+1,ip+1),5) .or. &
             BTEST(irr_quality_flagL(iw+1,ip+1),6) .or. &
             BTEST(irr_quality_flagL(iw+1,ip+1),7) .or. &
             BTEST(irr_quality_flagL(iw+1,ip+1),8) .or. &
             BTEST(irr_quality_flagL(iw+1,ip+1),9) .or. &
             BTEST(irr_quality_flagL(iw+1,ip+1),10))
        if (.not. transient_check) pxl_warning = pxl_warning .or. &
             BTEST(irr_quality_flagL(iw+1,ip+1),3)

        if(pxl_error(iw+1,ip+1)) then
          qc(ip,:)=IBSET(qc(ip,:),11)
          fs(iw,ip)=0.
          write(logmsg,"(A29,I6,2X,F10.6)") 'bad irradiance scan position ', &
               ip, ws(iw,ip)
          call tell_log(1,logmsg)
        endif
        if(pxl_warning) qc(ip,:)=IBSET(qc(ip,:),12)

      enddo !iw
    enddo !ip

  end subroutine bad_irrad_lambda

end module m_lambda_qual

