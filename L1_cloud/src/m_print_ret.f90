module m_print_ret

contains

  subroutine print_ret()

    use m_cloud_pres_mod
    use m_vars, ONLY: iprt, squeeze, lat, lon, sza, chlcl, cloud_clear, &
         wave_short, wave_long, iLine, cld_frac_min, refl_clr_oc, get_cloud_frac
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
        if (.not. cloud_clear) then
          print *,'b_in '
          write (6,fmt101) b_i
        endif
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
      print *, 'np0, np, irc, ir, nt, j, l'
      print *, np0, np, irc, ir, nt, j, l
      print *, 'i_np0, nc, irco, nt_o, j_o, l_o'
      print *, i_np0, nc, irco, nt_o, j_o, l_o
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
      if (add_oc .and. ((.not. ret_chl .and. iter == 0) .or. ret_chl)) then
        if (ret_chl .and. reflec < 0.2) then
          print *,'refl_clr_oc, refl_oc, reflec ',refl_clr_oc, refl_oc, reflec
          print *,'ring_ocs 1 ',ic1
          write(6,fmt100) ring_ocs(ic1,:)
          print *,'ring_ocs 2 ',ic2
          write(6,fmt100) ring_ocs(ic2,:)
          print *,'ring_oc'
          write(6,fmt100) ring_oc
          print *,'rad_clr_oc'
          write(6,fmt100) rad_clr_oc
          if (ret_chl) then
            print *,'h_chl'
            write(6,fmt100) h(:,nst-nterms-2-nsh)
          endif
          !stop
        endif
      endif
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
      if (ret_chl) then
        print *,'h(chloro,*)'
        write(6,fmt101) h(:,nst-nterms-2-nsh)
      endif
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
