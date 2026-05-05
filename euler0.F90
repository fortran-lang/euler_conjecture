program euler0
   ! 27^5 + 84^5 + 110^5 + 133^5 = 144^5 = 61917364224
   implicit none

   integer, parameter :: ik = selected_int_kind(18)  ! extended precision integer kind.
   integer, parameter :: imax = 6207                 ! search for solutions i=1:imax
   integer, parameter :: findmax = 1                 ! max number of solutions.

   real(selected_real_kind(14)) :: time0, time1
   integer :: nfound
   integer(ik) :: i, j, k, l, m
   character(*), parameter :: tfmt='(/a,1x,f0.4)'
   character(*), parameter :: fmtf = '(i0,": ",3(i0,"^5 + "),3(i0:"^5 = "))'

   call cpu_time( time0 )
   nfound = 0
   outer: do i = 1, imax
      do j = 1, i
         do k = 1, i
            do l = 1, i
               do m = 1, i
                  if ( j**5 + k**5 + l**5 + m**5 == i**5 ) then
                     nfound = nfound + 1
                     print fmtf, nfound, j, k, l, m, i, i**5
                     if ( nfound >= findmax ) exit outer
                  endif                  
               enddo
            enddo
         enddo
      enddo
   enddo outer
   call cpu_time( time1 )
   print tfmt, 'cpu time:', (time1-time0)
end program euler0
