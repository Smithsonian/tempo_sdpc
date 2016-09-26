!>Interpolation routines for pressure, radiance
module m_interp_pres

contains

  !-------------------------------------------------------------------------
  !
  ! !ROUTINE:  interp_pres
  ! 
  ! !DESCRIPTION: 
  !> Interpolate 3rd param at chosen wavelength to match
  !> pressure, and update jacobian 
  !
  ! !INPUT PARAMETERS:   
  !> @param[in]   ix1  cloud pressure index
  !> @param[in]   ix2  cloud pressure index+1
  !> @param[in]   pres_int  current value of pressure
  !> @param[in]   pres  vector of pressures to interpolate against
  !
  ! !OUTPUT PARAMETERS:  
  !> @param[out] rad  total reflected radiance at top of atmosphere (I think)
  !> @param[out] jacob  jacobian
  !
  ! !SEE ALSO:  
  !
  ! !REVISION HISTORY: 
  !
  !> @author  05Jan01   Joiner     original fortran 90
  !> @author  12Aug14  O'Sullivan  added documentation, 
  !>                               some guesswork involved
  !
  !-------------------------------------------------------------------------
  subroutine interp_pres(ix1, ix2, rad, pres, pres_int, jacob)

    use m_cloud_pres_mod, ONLY: temp2D
    implicit none

    ! input/output variables
    integer, intent(in) :: ix1, ix2
    real (KIND=8), intent(in) :: pres_int
    real (KIND=8), dimension(:), intent(in) :: pres
    real (KIND=8), dimension(:), intent(out) :: rad, jacob
    ! local variables
    real (KIND=8) :: temp1, temp2
    integer :: i

    !**************************************************************************

    temp1=(pres_int-pres(ix1))
    temp2=(pres(ix2)-pres(ix1))

    do i=1,size(rad)
      jacob(i)=(temp2D(2,i)-temp2D(1,i))/temp2 
    enddo
    if (temp1 /= 0) then
      do i=1,size(rad)
        rad(i)=temp2D(1,i) + (jacob(i))*temp1 
      enddo
    else
      rad=temp2D(1,:)
    endif

  end subroutine interp_pres




  !-------------------------------------------------------------------------
  !
  ! !ROUTINE:  interp_rads
  ! 
  ! !DESCRIPTION: 
  !> Interpolates radiance parameters to match the
  !>               current pressure value (or so it appears)
  !
  ! !INPUT PARAMETERS:   
  !> @param[in]    ix1 cloud pressure index
  !> @param[in]    ix2 cloud pressure index+1
  !> @param[in]   pres vector of pressure values
  !> @param[in]   pres_int current pressure value
  !> @param[in]   i0_1 backscattered intensity for ix1
  !> @param[in]   i0_2 backscattered intensity for ix2
  !> @param[in]   tr_1 transmittance factor for ix1 ?
  !> @param[in]   tr_2 transmittance factor for ix2 ?
  !> @param[in]   sb_1 surface light lost to scattering for ix1 ?
  !> @param[in]   sb_2 surface light lost to scattering for ix2 ?
  !
  ! !OUTPUT PARAMETERS:  
  !> @param[out]   i0  backscattered intensity in pixel
  !> @param[out]   tr  transmittance factor in pixel?
  !> @param[out]   sb  surface light lost to scattering in pixel?
  !
  ! !SEE ALSO:  
  !   m_read_tables.f90
  !
  ! !REVISION HISTORY: 
  !
  !> @author  05Jan01   Joiner     original fortran 90
  !> @author  12Aug14  O'Sullivan  added documentation, some guesswork involved
  !
  !-------------------------------------------------------------------------
  subroutine interp_rads(ix1, ix2, pres, pres_int, i0_1, i0_2, sb_1, sb_2, &
       tr_1, tr_2, i0, sb, tr )

    implicit none

    !input/ouput variables
    integer, intent(in)  :: ix1, ix2
    real (KIND=8), intent(in)  :: pres_int, i0_1, i0_2, sb_1, sb_2, tr_1, tr_2
    real (KIND=8), dimension(:), intent(in) :: pres
    real (KIND=8), intent(out) :: i0, sb, tr 
    !local variables
    real (KIND=8) :: temp

    !**************************************************************************

    temp=(pres(ix2)-pres_int)/(pres(ix2)-pres(ix1))
    i0=i0_2 - (i0_2-i0_1)*temp 
    sb=sb_2 - (sb_2-sb_1)*temp 
    tr=tr_2 - (tr_2-tr_1)*temp 

  end subroutine interp_rads

end module m_interp_pres
