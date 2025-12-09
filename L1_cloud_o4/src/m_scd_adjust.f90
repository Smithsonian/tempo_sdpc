!*************
module m_scd_adjust
!*************
! addition for O4 scd temperature dependence correction

contains

   !-----------------------------------------

   subroutine scd_adjust_gmi(pp,tt,cpp,scdm,scdadj,t8p)

   use m_vars, only: gmi_np, a263, b263, a293, b293, TrefO4

   real(kind=4), dimension(gmi_np+1), intent(in) :: pp
   real(kind=4), dimension(gmi_np), intent(in) :: tt
   real(kind=4), intent(in) :: cpp ! pressure to output t8p
   real, intent(in) :: scdm ! reference scd
   real, intent(out) ::scdadj, t8p ! adjusted SCD, T at cpp

   !local variables
   integer(kind=4):: kk, kkfound
   real(kind=4) :: ptop, pbot, wtop, wbot, ttop, tbot, ppmax
   real(kind=4) :: y1, y2, w1,w2

   t8p = -999.
   ptop = -999.
   pbot = -999.
   kkfound = -9

   ! initialize 
   scdadj = scdm

   !scdm<0 should have been skipped,just a safeguard
   if (scdm .lt. 0.) then
       return !in this case t8p=-999.
   endif

   ppmax = pp(gmi_np)

   ! get temperatre at cpp
   ! P-T profiles are from TOA to BOA
   do kk = 1, gmi_np
      if ((cpp .ge. pp(kk)) .and. (cpp .lt. pp(kk+1))) then
          kkfound = kk
          ptop = pp(kk)
          wtop = cpp - ptop
          ttop = tt(kk)
          pbot = pp(kk+1)
          wbot = pbot - cpp
          exit
      endif
   enddo

   if (kkfound .lt. 0) then ! did not find a layer
      if (cpp .lt. pp(1)) then
          t8p = tt(1)
      else if (cpp .gt. ppmax) then !assume same T near BOA
          t8p = tt(gmi_np)
      end if
   else ! found a layer
      if (cpp .ge. ppmax) then ! between BOA and pp(gmi_np)
         tbot = tt(gmi_np)
      else ! between pp(kkfound) and pp(kkfound+1)
         tbot = tt(kkfound+1)
      end if
      ! linear interpolation
      t8p = (tbot*wtop + ttop*wbot) / (wtop + wbot)
   end if

   ! ensure t8p is within Finkenzeller T range
   if (t8p .lt. 223.) t8p = 223.
   if (t8p .gt. 293.) t8p = 293.

   ! original L2 SCD are for TrefO4=223K
   ! T correction within range, clip on both ends
   ! scdm and scdadj need to be normalized by 1.e43 to use the a & b coeffs
   if (t8p .le. 223.) then
       scdadj = scdm 
   else if ((t8p .gt. 223.) .and. (t8p .le. 263.)) then
       y1 = scdm 
       y2 = a263 * scdm + b263
       w1 = t8p - 223.
       w2 = 263. - t8p
       scdadj = (y1*w2 + y2*w1)/(w1+w2)
   else if ((t8p .gt. 263.) .and. (t8p .le. 293.)) then 
       y1 = a263 * scdm + b263
       y2 = a293 * scdm + b293
       w1 = t8p - 263.
       w2 = 293. - t8p
       scdadj = (y1*w2 + y2*w1)/(w1+w2)
   else if (t8p .gt. 293.) then
       scdadj = a293 * scdm + b293
   endif

   !ensure scdadj is positive
   if (scdadj .lt. 0.) then
      scdadj = scdm
      t8p = TrefO4 
   endif

   end subroutine scd_adjust_gmi
   !-----------------------------------------

   subroutine scd_adjust_geos(pp,tt,cpp,scdm,scdadj,t8p)

   use m_vars, only: geos_np, &
                     a263, b263, a293, b293, TrefO4

   real(kind=4), dimension(:), intent(in) :: pp ! pp(geos_np+1)
   real(kind=4), dimension(:), intent(in) :: tt ! tt(geos_np)
   real, intent(in) :: cpp ! pressure to output t8p
   real, intent(in) :: scdm ! reference scd
   real, intent(out) ::scdadj, t8p ! adjusted SCD, T at cpp

   !local variables
   integer(kind=4):: kk, kkfound
   real :: ptop, pbot, wtop, wbot, ttop, tbot, ppmax
   real :: y1, y2, w1,w2

   t8p = -999.
   ptop = -999.
   pbot = -999.
   kkfound = -9

   ! initialize
   scdadj = scdm
   if (scdm .lt. 0.) then !safeguard, should have been skipped
       return ! t8p = -999. when this happens
   endif

   ppmax = pp(geos_np)

   ! get temperatre at cpp
   ! P-T profiles are from TOA to BOA
   do kk = 1, geos_np
      if ((cpp .ge. pp(kk)) .and. (cpp .lt. pp(kk+1))) then
          kkfound = kk
          ptop = pp(kk)
          wtop = cpp - ptop
          ttop = tt(kk)
          pbot = pp(kk+1)
          wbot = pbot - cpp
          exit
      endif
   enddo

   if (kkfound .lt. 0) then ! did not find a layer
      if (cpp .lt. pp(1)) then
          t8p = tt(1)
      else if (cpp .gt. ppmax) then !assume same T near BOA
          t8p = tt(geos_np)
      end if
   else ! found a layer
      if (cpp .ge. ppmax) then ! between BOA and pp(geos_np)
         tbot = tt(geos_np)
      else ! between pp(kkfound) and pp(kkfound+1)
         tbot = tt(kkfound+1)
      end if
      ! linear interpolation
      t8p = (tbot*wtop + ttop*wbot) / (wtop + wbot)
   end if

   ! ensure t8p is within O4 T range(223., 293.)K
   ! thus, t8p is for cpp only when within range
   if (t8p .lt. 223.) t8p = 223.
   if (t8p .gt. 293.) t8p = 293.

   ! original L2 SCD are for TrefO4=223K
   ! T correction within Finkenzeller T range, clip on both ends
   ! scdm and scdadj need to be normalized by 1.e43 to use the a & b coeffs
   if (t8p .le. 223.) then ! actually only t8p=223 will happen 
       scdadj = scdm
   else if ((t8p .gt. 223.) .and. (t8p .le. 263.)) then
       y1 = scdm 
       y2 = a263 * scdm + b263
       w1 = t8p - 223.
       w2 = 263. - t8p
       scdadj = (y1*w2 + y2*w1)/(w1+w2)
   else if ((t8p .gt. 263.) .and. (t8p .le. 293.)) then
       y1 = a263 * scdm + b263
       y2 = a293 * scdm + b293
       w1 = t8p - 263.
       w2 = 293. - t8p
       scdadj = (y1*w2 + y2 * w1)/(w1 + w2)
   else if (t8p .gt. 293.) then !should not happen, safeguard
       scdadj = a293 * scdm + b293
   endif

   ! ensure scdadj is positive
   if (scdadj .lt. 0.) then
      scdadj = scdm !scdm>0. otherwise should have returned
      t8p = TrefO4
   endif

   end subroutine scd_adjust_geos
!------------

!************
end module m_scd_adjust
!************

