module m_trilin

private

public trilinear, bilin

contains

function bilin( ix, jy ) result(xy)   

use m_cloud_pres_mod, ONLY: temp3D
implicit NONE          ! *** IDL2F9O ***

real (KIND=8), intent(in) :: ix, jy
real (KIND=8), dimension(size(temp3D,dim=3)) :: xy

integer :: i, j, ip, jp, ii
real (KIND=8) :: dx, dy, dx1, dy1
real (KIND=8) :: d1, d2, d3, d4

        i=int(ix, kind=4) 
        j=int(jy, kind=4)   
        ip=i+1   
        jp=j+1   
        ip=minval((/ip,size(temp3D(:,1,1))/))   
        jp=minval((/jp,size(temp3D(1,:,1))/))   
        dx=ix-i 
        dy=jy-j   
        dx1=1.-dx 
        dy1=1.-dy   
        d1=dx1*dy1
        d2=dx1*dy
        d3=dx *dy1
        d4=dx *dy 
    !print *, 'bilin ',dx,dy,dx1,dy1,i,j,ip,jp
    !print *, 'table ',table( i, j,: )
!     write(6,*) size(temp3D,1), size(temp3D,2), size(temp3D,3)
!     write(6,*) size(xy), i,j,ip,jp
!     write(6,*) dx1,dy1,dx,dy
     do ii=1, size(xy)
!       print *, ii, size(xy)
!       print *, temp3D(i,j,1)
       !stop
       xy(ii) =  temp3D( i, j,ii ) *d1       &   
         +       temp3D( i,jp,ii ) *d2       &   
         +       temp3D( ip,j,ii ) *d3       &   
         +       temp3D( ip,jp,ii )*d4   
     enddo
! ***********************************************************
end  function bilin
   
function trilinear( zp,xp,yp ) result(interp) 

use m_cloud_pres_mod, ONLY: table, temp3D
implicit NONE          

real (KIND=8), dimension(size(table,dim=4)) :: interp
real (KIND=8), intent(in) :: xp, yp, zp

integer :: iz, iz2
real (KIND=8), dimension(size(table,dim=4)) :: interp1, interp2
   
!find bracketing z 
!==============================
iz=int(zp)
   
!bounds check
!============
iz=maxval((/iz,1/))   
! JJ bug fix
!OK for v1+ of OMCLDRR
!iz=minval((/iz,size(table(1,1,:,1))-1/))   
iz=minval((/iz,size(table(:,1,1,1))-1/))   
iz2=iz+1   
   
!interpolate in x, y first
!==========================
!interp1 = bilin(table(:,:,iz,:),xp,yp)   
!interp2 = bilin(table(:,:,iz2,:),xp,yp) 
temp3D  => table(iz,:,:,:)
interp1 = bilin(xp,yp)   
temp3D  => table(iz2,:,:,:)
interp2 = bilin(xp,yp) 

!now interpolate in z
!======================
interp =(zp-iz)*(interp2-interp1)+interp1
   
! ***********************************************************
end function trilinear
   
end module m_trilin   
