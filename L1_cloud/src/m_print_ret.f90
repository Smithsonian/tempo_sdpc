module m_print_ret

contains

  subroutine print_ret()

    use m_cloud_pres_mod
    use m_vars, ONLY: iprt, squeeze, lat, lon, sza, chlcl, & 
         wave_short, wave_long, iLine, cld_frac_min, refl_clr_oc, &
         get_cloud_frac
    implicit none
    !-------------------------------------------------------------------------
    !         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
    !-------------------------------------------------------------------------
    !BOP
    !
    ! !ROUTINE:  print_ret
    ! 
    ! !DESCRIPTION: print_ret prints out Jacobian and other parameters of
    !               retrievals		
    !
    ! !CALLING SEQUENCE: 
    !
    !        call print_ret
    !     
    ! !INPUT PARAMETERS:   
    !
    ! !OUTPUT PARAMETERS:  
    !
    ! !SEE ALSO:  
    !
    ! !REVISION HISTORY: 
    !
    !  05Jan01   Joiner     original fortran 90
    !
    !EOP
    !-------------------------------------------------------------------------
    !
    !formats for write statements
    character (len=8) :: fmt100="(6f12.5)", fmt101="(6e12.3)"

    !**************************************************************************

    if (iprt >= 3) then
      print *, 'chisq, bias, std'
      write(6,fmt101) chisq, bias, std   
      if (iter == 0) then
        print *,'b_i '
        write (6,fmt101) b_i
        print *, 'lat, lon, sza, chl. clim.'
        write(6,fmt101) lat(ip+1,iLine), lon(ip+1,iLine), sza(ip,iLine), chlcl(ip)
        print *,'x_fg '
        write (6,fmt100) x_fg
      endif
      print *,'x, iter= ',iter
      write (6,fmt100) x
    endif
    if (iprt >= 6) print *,'chl ',chloro
    if (iprt >= 7 .and. iter == 1) then
      print *, iter, wave_long, wave_short
      print *, 'np0, np, nt, j, l'
      print *, np0, np, nt, j, l
      print *, 'i_np0, nc, nt_o, j_o, l_o'
      print *, i_np0, nc, nt_o, j_o, l_o
      print *, ind(0), ind(nobs-1)
      print *,'y_obs'
      write(6,fmt100) y_obs
      print *,'waves'
      write(6,fmt100) waves
      print *,'rad_cld'
      write(6,fmt100) rad_cld
      if (cld_frac > cld_frac_min .and. cld_frac < 1.0 .and. &
           get_cloud_frac .and. iprt >= 7) then
        print *,'rad_clr'
        write(6,fmt100) rad_clr
        print *,'ring_clr'
        write(6,fmt100) ring_clr
      endif
      print *,'ring_cld'
      write(6,fmt100) ring_cld
      print *, 'ix1, ix2, np, x(0,1), psurf,iter,  ip'
      print *, ix1, ix2, np, x(0,1), psurf,iter,  ip
      print *,'rad_clds 1 ',ix1
      write(6,fmt100) rad_clds(ix1,:)
      print *,'rad_clds 2 ',ix2
      write(6,fmt100) rad_clds(ix2,:)
      print *,'ring_clds 1 '
      write(6,fmt100) ring_clds(ix1,:)
      print *,'ring_clds 2 '
      write(6,fmt100) ring_clds(ix2,:)
      print *,'h(0,*)'
      write(6,fmt100) h(:,0)
      if (squeeze) then
        print *,'h(1,*)'
        write(6,fmt100) h(:,1)
      endif
      !  print *,'h(*,nst-4)'
      !  write(6,fmt100) h(:,nst-4)
      print *,'h(*,nst-3)'
      write(6,fmt100) h(:,nst-3)
      print *,'h(*,nst-2)'
      write(6,fmt100) h(:,nst-2)
      print *,'h(*,nst-1)'
      write(6,fmt100) h(:,nst-1)
      print *,'rad_tot'
      write(6,fmt100) rad_tot
      print *,'sflx'
      write(6,fmt101) sflx
      print *,'ycalc'
      write(6,fmt100) ycalc
      print *,'err_cov '
      write (6,fmt101) err_cov
    endif
    if (iprt >= 6) then
      print *,'y_calc '
      write (6,fmt100) y_calc
    endif


  end subroutine print_ret

end module m_print_ret
