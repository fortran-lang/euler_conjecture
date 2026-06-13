module ring_mod

   ! this module isolates all of the modulo ring arithmetic associated
   ! with the Euler conjecture search code.

   implicit none

   private
   public :: ik, ring_modulo_t, jloop, kloop, lloop

   ! ik is the kind value for an extended precision integer
   ! to compute i**5 values without overflow.

   integer, parameter :: ik = selected_int_kind(18)

   integer, parameter :: jloop=1, kloop=2, lloop=3   ! j, k, and l loop level enums.
   character, parameter :: cloop(jloop:lloop) = ['j', 'k', 'l']

   type level_t
      integer              :: nval    ! number of allowed values.
      integer, allocatable :: val(:)  ! (1:nval) list of allowed values.
      integer, allocatable :: map(:)  ! (0:N-1) map from mod(i,N) to contiguous index, ix=1:nval.
   end type level_t

   type allow_t   ! info associated with each modulo N value.
      integer                    :: N          ! modulus value.
      integer                    :: i5modNx    ! mapped mod(i**5,N) values
      type(level_t)              :: i          ! target i info.
      type(level_t), allocatable :: jkl(:,:)   ! (1:nval,jloop:lloop) j, k, and l loop level info.
   contains
      procedure :: print => print_allow
      procedure :: fill  => fill_allow
   end type allow_t

   type ring_modulo_t  ! overall ring modulo info.
      integer(ik)                :: kntt(jloop:lloop)  ! total count of modulo tests for each loop level.
      integer(ik)                :: kntf(jloop:lloop)  ! count of false modulo tests.
      type(allow_t), allocatable :: allow(:)           ! zero or more modulo tests.
   contains
      procedure :: init   => init_ring_modulo
      procedure :: test   => modulo_test
      procedure :: set_i5 => set_i5modN
      procedure :: print  => print_stats
   end type ring_modulo_t

   character(*), parameter :: cfmt = '(*(g0.3,1x))'

contains

   impure elemental subroutine fill_allow( allow, N )

      ! allocate and fill in the allow derived type based on the modulus N.

      ! the goal is to limit the j, k, and l loop indices when possible. this is
      ! achieved by using modulo arithmetic to limit the possibile
      ! matches of j**5 + k**5 + l**5 + m**5 with the target i**5.

      implicit none
      class(allow_t), intent(out) :: allow
      integer, intent(in)         :: N      ! modulus value with range up to about 1000.

      integer :: i, i5x, lx, kx, nval
      integer :: val(1:N)          ! temp work arrays.
      logical :: mask(0:N-1), jmask(0:N-1), kmask(0:N-1)

      allow%N       = N
      allow%i5modNx = -1  ! this is filled in later.

      ! compute mask(:) for the possible mod(i**5,N) values.

      mask(:) = .false.
      do i = 1, N
         mask( mod(int(i,ik)**5,int(N,ik)) ) = .true.
      enddo

      jmask(:) = mask(:)  ! these are also the possible j values.

      allocate( allow%i%map(0:N-1) )
      call set_val_map( mask, allow%i%map, val, nval )
      allow%i%nval = nval
      allow%i%val  = val(1:nval)  ! allocate on assignment.

      ! compute the kmask(:) array.
      ! the .true. entries are the only sk values that can be generated within the j loop.

      kmask(:) = .false.
      do i5x = 1, allow%i%nval
         call update_mask( N, allow%i%val(i5x), allow%i%val, kmask )
      enddo

      ! allocate and fill in the level l, k, and j info.

      allocate( allow%jkl(1:nval,jloop:lloop) )  ! j=>1, k=>2, l=>3.

      do i5x = 1, allow%i%nval
         ! fill in the level l info for this i5x.
         mask(:) = .false.
         call update_mask( N, allow%i%val(i5x), allow%i%val, mask )
         allocate( allow%jkl(i5x,lloop)%map(0:N-1) )
         call set_val_map( mask, allow%jkl(i5x,lloop)%map, val, nval )
         allow%jkl(i5x,lloop)%nval = nval
         allow%jkl(i5x,lloop)%val  = val(1:nval)

         ! fill in the level k info for this i5x.
         mask(:) = .false.
         do lx = 1, allow%jkl(i5x,lloop)%nval
            call update_mask( N, allow%jkl(i5x,lloop)%val(lx), allow%i%val, mask )
         enddo
         mask = mask .and. kmask  ! mask out any unreachable values.
         allocate( allow%jkl(i5x,kloop)%map(0:N-1) )
         call set_val_map( mask, allow%jkl(i5x,kloop)%map, val, nval )
         allow%jkl(i5x,kloop)%nval = nval
         allow%jkl(i5x,kloop)%val  = val(1:nval)

         ! fill in the level j info for this i5x.
         mask(:) = .false.
         do kx = 1, allow%jkl(i5x,kloop)%nval
            call update_mask( N, allow%jkl(i5x,kloop)%val(kx), allow%i%val(:), mask )
         enddo
         mask = mask .and. jmask  ! mask out any unreachable values.
         allocate( allow%jkl(i5x,jloop)%map(0:N-1) )
         call set_val_map( mask, allow%jkl(i5x,jloop)%map, val, nval )
         allow%jkl(i5x,jloop)%nval = nval
         allow%jkl(i5x,jloop)%val  = val(1:nval)
      enddo

      return
   end subroutine fill_allow

   subroutine print_stats( ring_modulo )
      ! print the modulo ring loop statistics.
      implicit none
      class (ring_modulo_t), intent(in) :: ring_modulo

      integer     :: loop
      integer(ik) :: total, rejected

      print cfmt
      print cfmt, 'modulo ring loop statistics with modvals(:)=', ring_modulo%allow(:)%N
      do loop = jloop, lloop
         total    = ring_modulo%kntt(loop)           ! total number of index tests.
         rejected = total - ring_modulo%kntf(loop)   ! number of rejected index values.
         print cfmt, cloop(loop), 'loop:', &
              & 'total=', total, 'rejected=', rejected, &
              & 'ratio=', real(rejected)/real(total)
      enddo

      return
   end subroutine print_stats

   impure elemental subroutine print_allow( allow )
      ! print the contents of a scalar or array allow structures with annotations.
      implicit none
      class(allow_t), intent(in) :: allow

      integer :: i5x, loop

      print cfmt, '    ==== Allowed index values for modulus N=', allow%N, '===='

      print cfmt
      print cfmt, 'Loop level i, nval=', allow%i%nval, &
           & 'nonzero ratio=', real(allow%i%nval)/allow%N
      print cfmt, 'Allowed modulo(i^5,N) values, val(:)=', allow%i%val

      do loop = jloop, lloop
         print cfmt
         print cfmt, 'Loop level', cloop(loop), 'info:'
         do i5x = 1, allow%i%nval
            print cfmt, 'For modulo(i^5,N)=', allow%i%val(i5x), 'nval=', allow%jkl(i5x,loop)%nval, &
              & 'nonzero ratio=', real(allow%jkl(i5x,loop)%nval)/real(allow%N)
            print cfmt, 'Allowed modulo(sum,N) values, val(:)=', allow%jkl(i5x,loop)%val
         enddo
      enddo
      print cfmt

      return
   end subroutine print_allow

   subroutine set_val_map( mask, map, val, nval )
      ! search mask(0:) for .true. entries, set map(0:), and save the indices in val(:).
      ! the arrays map(0:) and val(:) are inverse mappings.
      implicit none
      logical, intent(in)  :: mask(0:) ! (0:N-1).
      integer, intent(out) :: map(0:)  ! (0:N-1)
      integer, intent(out) :: val(:)   ! allow(1:nval) are the allowed values.
      integer, intent(out) :: nval     ! number of allowed values.

      integer :: p

      map(:) = 0
      nval   = 0
      do p = 0, ubound(mask,dim=1)
         if ( mask(p) ) then
            nval      = nval + 1
            val(nval) = p
            map(p)    = nval
         endif
      enddo

      return
   end subroutine set_val_map

   subroutine update_mask( N, m, val, mask )
      ! update mask(:) using m and val(:).
      ! this is the modulo N ring addition of m with the additive inverse of the elements of val(:).
      ! mask(k)==.true. if k is in the output set.
      implicit none
      integer, intent(in)    :: N        ! the current modulus.
      integer, intent(in)    :: m        ! the target value.
      integer, intent(in)    :: val(:)   ! set of input values.
      logical, intent(inout) :: mask(0:) ! updated ring value indexes.

      integer :: p

      do p = 1, size(val)
         mask( modulo( m - val(p), N ) ) = .true.
      enddo

      return
   end subroutine update_mask

   subroutine init_ring_modulo( ring_modulo, modvals )
      ! initialize, allocate, and fill the ring_modulo derived type.
      implicit none
      class(ring_modulo_t), intent(out) :: ring_modulo
      integer, intent(in)               :: modvals(:)

      ring_modulo%kntt(:) = 0   ! ring modulo total tests.
      ring_modulo%kntf(:) = 0   ! ring modulo false tests.

      allocate( ring_modulo%allow(size(modvals)) )

      call ring_modulo%allow%fill( modvals )  ! (:)

      return
   end subroutine init_ring_modulo

   subroutine set_i5modN( ring_modulo, i5 )
      ! set all of the mapped modulo(i**5,N) values.
      implicit none
      class(ring_modulo_t), intent(inout) :: ring_modulo
      integer(ik), intent(in)             :: i5       ! i**5 value.

      integer     :: p
      integer(ik) :: N

      do p = 1, size(ring_modulo%allow)
         N = ring_modulo%allow(p)%N
         ring_modulo%allow(p)%i5modNx = ring_modulo%allow(p)%i%map( modulo(i5,N) )
      enddo

      return
   end subroutine set_i5modN

   logical function modulo_test( ring_modulo, loop, sum )
      ! perform modulo tests on the loop index.
      implicit none
      class(ring_modulo_t), intent(inout) :: ring_modulo
      integer, intent(in)                 :: loop   ! loop level: j=>1, k=>2, l=>3.
      integer(ik), intent(in)             :: sum    ! j**5 + k**5 +... partial sum.

      integer     :: p, i5modNx
      integer(ik) :: N

      ring_modulo%kntt(loop) = ring_modulo%kntt(loop) + 1  ! total number of tests.

      ! loop in order.

      modulo_test = .false.
      do p = 1, size(ring_modulo%allow)
         N           = ring_modulo%allow(p)%N
         i5modNx     = ring_modulo%allow(p)%i5modNx
         modulo_test = ring_modulo%allow(p)%jkl(i5modNx,loop)%map(mod(sum,N)) .eq. 0
         if ( modulo_test ) return   ! return ASAP.
      enddo

      ring_modulo%kntf(loop) = ring_modulo%kntf(loop) + 1  ! number of .false. tests.

      return
   end function modulo_test

end module ring_mod

program euler2

   ! 133^5 + 110^5 + 84^5 + 27^5 = 144^5 = 61917364224

   use ring_mod, only: ik, ring_modulo_t, jloop, kloop, lloop

   implicit none

   integer, parameter :: imax = 6207          ! integer range to search i=1:imax.
   integer, parameter :: findmax = 5          ! max number of solutions to find.
   logical, parameter :: printring = .false.  ! print the modulo ring derived type info.
   integer, parameter :: modvals(*)=[775,275] ! modulo values. some good values are 11, 31, 100, 275, 775.

   real(selected_real_kind(14)) :: time0, time1
   integer :: i, j, k, l, m, nfound
   integer(ik) :: i5, sj, sk, sl, sm, diff
   integer(ik) :: n5(imax)
   character(*), parameter :: tfmt='(/a,1x,f0.4)'
   character(*), parameter :: fmtf = '(i0,": ",3(i0,"^5 + "),3(i0:"^5 = "))'

   type(ring_modulo_t) :: ring_modulo

   call cpu_time( time0 )

   call ring_modulo%init( modvals )
   if ( printring ) call ring_modulo%allow%print()

   nfound = 0
   outer: do i = 1, imax
      i5 = int(i,ik)**5
      n5(i) = i5
      call ring_modulo%set_i5( i5 )
      do j = index_min( i5/4 ), i
         sj = n5(j)
         diff = i5 - sj
         if ( diff < 3 ) exit
         if ( ring_modulo%test( jloop, sj ) ) cycle
         do k = index_min( diff/3 ), j
            sk = sj + n5(k)
            diff = i5 - sk
            if ( diff < 2 ) exit
            if ( n5(k) < diff/2 ) cycle  ! no solution in l=1:k.
            if ( ring_modulo%test( kloop, sk ) ) cycle
            do l = index_min( diff/2 ), k
               sl = sk + n5(l)
               diff = i5 - sl    ! for a solution, m**5 == diff is required.
               if ( diff < 1 ) exit
               if ( n5(l) < diff ) cycle  ! no solution in m=1:l.
               if ( ring_modulo%test( lloop, sl ) ) cycle

               mloop: do m = index_min( diff ), l   ! although written as a loop, m only cycles 1 or 2 times.
                  sm = sl + n5(m)
                  diff = i5 - sm
                  select case (diff)
                  case (:-1)   ! i5 < sm
                     exit mloop
                  case (0)     ! i5 == sm
                     nfound = nfound + 1
                     print fmtf, nfound, j, k, l, m, i, i5
                     if ( nfound >= findmax ) then
                        exit outer
                     else
                        exit mloop  ! exit m loop and continue searching.
                     endif
                  case default !(1:) sm < i5
                     cycle mloop
                  end select
               enddo mloop

            enddo
         enddo
      enddo
   enddo outer
   call ring_modulo%print()
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

end program euler2
