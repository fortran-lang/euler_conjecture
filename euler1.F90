program euler1
   ! 133^5 + 110^5 + 84^5 + 27^5 = 144^5 = 61917364224
   implicit none

   integer, parameter :: ik = selected_int_kind(18)  ! extended precision integer kind.
   integer, parameter :: imax = 6207                 ! search for solutions i=1:imax
   integer, parameter :: findmax = 5                 ! max number of solutions.

   real(selected_real_kind(14)) :: time0, time1
   integer :: i, j, k, l, m, nfound
   integer(ik) :: i5, sj, sk, sl, sm, diff
   integer(ik) :: n5(imax)
   character(*), parameter :: tfmt='(/a,1x,f0.4)'
   character(*), parameter :: fmtf = '(i0,": ",3(i0,"^5 + "),3(i0:"^5 = "))'

   call cpu_time( time0 )
   nfound = 0
   outer: do i = 1, imax
      i5    = int(i,ik)**5
      n5(i) = i5
      do j = index_min( i5/4 ), i
         sj = n5(j)
         diff = i5 - sj
         if ( diff < 3 ) exit
         do k = index_min( diff/3 ), j
            sk = sj + n5(k)
            diff = i5 - sk
            if ( diff < 2 ) exit
            do l = index_min( diff/2 ), k
               sl = sk + n5(l)
               diff = i5 - sl
               if ( diff < 1 ) exit
               
               mloop: do m = index_min( diff ), l    ! although written as a loop, m only cycles 1 or 2 times.
                  sm = sl + n5(m)
                  diff = i5 - sm
                  select case (diff)
                  case (:-1)  ! i5 < sm.
                     exit mloop
                  case (0)   ! sm == i5.
                     nfound = nfound + 1
                     print fmtf, nfound, j, k, l, m, i, i5
                     if ( nfound >= findmax ) then
                        exit outer
                     else
                        exit mloop  ! exit m loop and continue searching.
                     endif
                  case default    ! (1:)  sm < i5
                     cycle mloop
                  end select
               enddo mloop
               
            enddo
         enddo
      enddo
   enddo outer
   call cpu_time( time1 )
   print tfmt, 'cpu time:', (time1-time0)

contains

   pure function index_min( v ) result( r )
      ! return an integer value r such that r**5 <= v.
      ! it need not be the largest value, but the closer to v the better.
      implicit none
      integer                 :: r
      integer(ik), intent(in) :: v

      r = int( real(v)**0.2 )       ! the floating point result can be too
      do while (int(r,ik)**5 > v)   ! large due to rounding. correct it here.
         r = r - 1
      enddo
      r = max( 1, r )               ! enforce lower bound of 1.

      return
   end function index_min

end program euler1
