!*************
module m_scd_adjust
!*************
!hqw addition for O4 scd temperature dependence correction

contains

   !-----------------------------------------

   subroutine scd_adjust_gmi(pp,tt,cpp,scdm,scdadj,t8p)

   use m_vars, only: gmi_np, a203, b203, a233, b233, &
                      a253, b253, a293, b293

   real(kind=4), dimension(gmi_np+1), intent(in) :: pp
   real(kind=4), dimension(gmi_np), intent(in) :: tt
   real(kind=4), intent(in) :: cpp ! pressure to output t8p
   real, intent(in) :: scdm ! 273K reference scd
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
   if (scdm .lt. 0.) then
       return
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

   ! ensure t8p is within 10K of (203., 303.)K
   if (t8p .lt. 193.) t8p = 193.
   if (t8p .gt. 303.) t8p = 303.

   ! original L2 SCD are for 273K
   ! T correction within [203, 293]K, clip on both ends
   ! scdm and scdadj need to be normalized by 1.e43 to use the a & b coeffs
   if (t8p .le. 203.) then
       scdadj = a203 * scdm + b203
   else if ((t8p .gt. 203.) .and. (t8p .le. 233.)) then
       y1 = a203 * scdm + b203
       y2 = a233 * scdm + b233
       w1 = t8p - 203.
       w2 = 233. - t8p
       scdadj = (y1*w2 + y2*w1)/(w1+w2)
   else if ((t8p .gt. 233.) .and. (t8p .le. 253.)) then
       y1 = a233 * scdm + b233
       y2 = a253 * scdm + b253
       w1 = t8p - 233.
       w2 = 253. - t8p
       scdadj = (y1*w2 + y2 * w1)/(w1 + w2)
   else if ((t8p .gt. 253.) .and. (t8p .le. 273.)) then
       y1 = a253 * scdm + b253
       y2 = scdm
       w1 = t8p - 253.
       w2 = 273. - t8p
       scdadj = (y1*w2 + y2*w1) / (w1 + w2)
   else if ((t8p .gt. 273.) .and. (t8p .le. 293.)) then
       y1 = scdm
       y2 = a293 * scdm + b293
       w1 = t8p - 273.
       w2 = 293. - t8p
       scdadj = (y1*w2+y2*w1)/(w1 + w2)
   else if (t8p .gt. 293.) then
       scdadj = a293 * scdm + b293
   endif

   !ensure scdadj is positive
   if (scdadj .lt. 0.) then
      scdadj = scdm
      t8p =  273.
   endif

   end subroutine scd_adjust_gmi
   !-----------------------------------------

   subroutine scd_adjust_geos(pp,tt,cpp,scdm,scdadj,t8p)

   use m_vars, only: geos_np, a203, b203, a233, b233, &
                      a253, b253, a293, b293

   real(kind=4), dimension(geos_np+1), intent(in) :: pp
   real(kind=4), dimension(geos_np), intent(in) :: tt
   real, intent(in) :: cpp ! pressure to output t8p
   real, intent(in) :: scdm ! 273K reference scd
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
   if (scdm .lt. 0.) then
       return
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

   ! ensure t8p is within 10K of (203., 303.)K
   if (t8p .lt. 193.) t8p = 193.
   if (t8p .gt. 303.) t8p = 303.

   ! original L2 SCD are for 273K
   ! T correction within [203, 293]K, clip on both ends
   ! scdm and scdadj need to be normalized by 1.e43 to use the a & b coeffs
   if (t8p .le. 203.) then
       scdadj = a203 * scdm + b203
   else if ((t8p .gt. 203.) .and. (t8p .le. 233.)) then
       y1 = a203 * scdm + b203
       y2 = a233 * scdm + b233
       w1 = t8p - 203.
       w2 = 233. - t8p
       scdadj = (y1*w2 + y2*w1)/(w1+w2)
   else if ((t8p .gt. 233.) .and. (t8p .le. 253.)) then
       y1 = a233 * scdm + b233
       y2 = a253 * scdm + b253
       w1 = t8p - 233.
       w2 = 253. - t8p
       scdadj = (y1*w2 + y2 * w1)/(w1 + w2)
   else if ((t8p .gt. 253.) .and. (t8p .le. 273.)) then
       y1 = a253 * scdm + b253
       y2 = scdm
       w1 = t8p - 253.
       w2 = 273. - t8p
       scdadj = (y1*w2 + y2*w1) / (w1 + w2)
   else if ((t8p .gt. 273.) .and. (t8p .le. 293.)) then
       y1 = scdm
       y2 = a293 * scdm + b293
       w1 = t8p - 273.
       w2 = 293. - t8p
       scdadj = (y1*w2+y2*w1)/(w1 + w2)
   else if (t8p .gt. 293.) then
       scdadj = a293 * scdm + b293
   endif

   !ensure scdadj is positive
   if (scdadj .lt. 0.) then
      scdadj = scdm
      t8p =  273.
   endif

   end subroutine scd_adjust_geos
!------------

!************
end module m_scd_adjust
!************

