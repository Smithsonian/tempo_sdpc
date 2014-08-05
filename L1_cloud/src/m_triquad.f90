module m_triquad

  private

  public triquad, biquad

contains

  function interpol(ix,iy) result(p)

    use m_cloud_pres_mod, ONLY: temp3D
    implicit none

    !INPUT/OUTPUT variables
    real (KIND=8), intent(in) :: ix
    integer, intent(in) :: iy
    real (KIND=8), dimension(size(temp3D,dim=3)) :: p
    !local variables
    integer :: s, ii
    real (KIND=8) :: x0, x1, x2

    s=(ix) 
    s=maxval((/s,2/))
    s=minval((/s,size(temp3D(:,1,1))-1/))   
    x1=s
    x0=x1-1
    x2=x1+1
    do ii=1, size(p)
      p(ii) = temp3d(s-1,iy,ii) * (ix-x1) * (ix-x2) / ((x0-x1) * (x0-x2)) + &
           temp3d(s,iy,ii) *   (ix-x0) * (ix-x2) / ((x1-x0) * (x1-x2)) + &
           temp3d(s+1,iy,ii) * (ix-x0) * (ix-x1) / ((x2-x0) * (x2-x1))
    enddo

  end  function interpol



  function biquad( ix, jy ) result(xy)   

    use m_cloud_pres_mod, ONLY: temp3D
    implicit NONE          ! *** IDL2F9O ***

    !INPUT/OUTPUT variables
    real (KIND=8), intent(in) :: ix, jy
    real (KIND=8), dimension(size(temp3D,dim=3)) :: xy
    !local variables
    integer :: ii, s
    real (KIND=8), dimension(size(temp3D,dim=3)) :: p0, p1, p2
    integer :: y0, y1, y2


    s=(jy) 
    s=maxval((/s,2/))
    s=minval((/s,size(temp3D(1,:,1))-1/))   
    y1=s
    y0=y1-1
    y2=y1+1
    p0=interpol(ix,y0)
    p1=interpol(ix,y1)
    p2=interpol(ix,y2)
    do ii=1, size(xy)
      xy(ii) = p0(ii) * (jy-y1) * (jy-y2) / ((y0-y1) * (y0-y2)) + &
           p1(ii) * (jy-y0) * (jy-y2) / ((y1-y0) * (y1-y2)) + &
           p2(ii) * (jy-y0) * (jy-y1) / ((y2-y0) * (y2-y1))
    enddo

  end  function biquad



  function triquad( zp,xp,yp ) result(interp) 

    use m_cloud_pres_mod, ONLY: table, temp3D
    implicit NONE          

    !INPUT/OUTPUT variables
    real (KIND=8), dimension(size(table,dim=4)) :: interp
    real (KIND=8), intent(in) :: xp, yp, zp
    !local variables
    integer :: iz, z0, z1, z2, s, ii
    real (KIND=8), dimension(size(table,dim=4)) :: interp0, interp1, interp2

    !find bracketing z 
    !==============================
    s=(zp) 
    s=maxval((/s,2/))
    s=minval((/s,size(temp3D(1,:,1))-1/))   
    z1=s
    z0=z1-1
    z2=z1+1

    !interpolate in x, y first
    !==========================
    temp3D  => table(z0,:,:,:)
    interp0 = biquad(xp,yp)   
    temp3D  => table(z1,:,:,:)
    interp1 = biquad(xp,yp) 
    temp3D  => table(z2,:,:,:)
    interp2 = biquad(xp,yp) 

    !now interpolate in z
    !======================
    iz=zp
    do ii=1, size(interp)
      interp(ii) = interp0(ii) * (iz-z1) * (iz-z2) / ((z0-z1) * (z0-z2)) + &
           interp1(ii) * (iz-z0) * (iz-z2) / ((z1-z0) * (z1-z2)) + &
           interp2(ii) * (iz-z0) * (iz-z1) / ((z2-z0) * (z2-z1))
    enddo

  end function triquad

end module m_triquad
