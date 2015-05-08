subroutine matgen(imat, a, m, n, d, info)
  implicit none

  real(kind(1.d0)), allocatable :: a(:,:)
  integer                       :: imat, m, n, info

  real(kind(1.d0)), allocatable :: work(:), d(:)
  integer, allocatable          :: aux(:)
  character(len=3)              :: path
  character                     :: type, dist, pack
  integer                       :: kl, ku, mode, i
  real(kind(1.d0))              :: anorm, cndnum, rcondc
  integer                       :: iseed(4)
  real(kind(1.d0)), parameter   :: done=1.d0, dmone=-1.d0

  iseed = (/1988, 1989, 1990, 1991/)
  info  = 0
  path(1:1) = 'Double precision'
  path(2:3) = 'SY'

  select case(imat)
  case(1)
     dist   = 'U'  ! uniform distribution
     type   = 'N'  ! symmetric, positive-definite
     mode   = 0    ! use D
     cndnum = n    ! the condition number
     anorm  = n    ! the matrix norm
     kl     = m-1  ! lower bandwidth
     ku     = n-1  ! upper bandwidth
     pack   = 'N'  ! no pack
     write(*,'("=======================")')
     write(*,'("Dist  : ",a1)')dist
     write(*,'("Type  : ",a1)')type
     write(*,'("Mode  : ",i1)')mode
     write(*,'("Anorm : ",es10.3)')anorm
     write(*,'("Cond  : ",es10.3)')cndnum
     write(*,'(" ")')
      
     allocate(a(m,n), work(3*max(m,n)), d(n), stat=info)
     if(info .ne. 0) return
    
     do i=1,n
        d(i) = n-i+1
     end do

     call dlatms( m, n, dist, iseed, type, &
          & d, mode, cndnum, anorm, kl, ku, &
          & pack, a, m, work, info )
     
     deallocate(work, stat=info)

     if(info .ne. 0) return

  case(2)
     dist   = 'N'  ! uniform distribution
     type   = 'N'  ! symmetric, positive-definite
     mode   = 3    ! D(I)=COND**(-(I-1)/(N-1))
     cndnum = 1e4 ! the condition number
     anorm  = 50.d0 ! the matrix norm
     ku     = n-1  ! upper bandwidth
     kl     = m-1  ! lower bandwidth
     pack   = 'N'  ! no pack
     write(*,'("=======================")')
     write(*,'("Dist  : ",a1)')dist
     write(*,'("Type  : ",a1)')type
     write(*,'("Mode  : ",i1)')mode
     write(*,'("Anorm : ",es10.3)')anorm
     write(*,'("Cond  : ",es10.3)')cndnum
     write(*,'(" ")')
      
     allocate(a(m,n), work(3*max(m,n)), d(n), stat=info)
     if(info .ne. 0) return
    
     call dlatms( m, n, dist, iseed, type, &
          & d, mode, cndnum, anorm, kl, ku, &
          & pack, a, m, work, info )
     if(info .ne. 0) return
     
     deallocate(work, stat=info)

  case(3)
     dist   = 'N'  ! uniform distribution
     type   = 'N'  ! symmetric, positive-definite
     mode   = 3    ! D(I)=COND**(-(I-1)/(N-1))
     cndnum = 1e5  ! the condition number
     anorm  = 10.d0 ! the matrix norm
     ku     = n-1  ! upper bandwidth
     kl     = m-1  ! lower bandwidth
     pack   = 'N'  ! no pack
     write(*,'("=======================")')
     write(*,'("Dist  : ",a1)')dist
     write(*,'("Type  : ",a1)')type
     write(*,'("Mode  : ",i1)')mode
     write(*,'("Anorm : ",es10.3)')anorm
     write(*,'("Cond  : ",es10.3)')cndnum
     write(*,'(" ")')
      
     allocate(a(m,n), work(3*max(m,n)), d(n), stat=info)
     if(info .ne. 0) return
    
     call dlatms( m, n, dist, iseed, type, &
          & d, mode, cndnum, anorm, kl, ku, &
          & pack, a, m, work, info )
     if(info .ne. 0) return
     
     deallocate(work, stat=info)

  case(4)
     dist   = 'N'  ! uniform distribution
     type   = 'N'  ! symmetric, positive-definite
     mode   = 5    ! D(i)=1 - (i-1)/(N-1)*(1 - 1/COND)
     cndnum = 1e6  ! the condition number
     anorm  = real(max(m,n)) ! the matrix norm
     ku     = n-1  ! upper bandwidth
     kl     = m-1  ! lower bandwidth
     pack   = 'N'  ! no pack
     write(*,'("=======================")')
     write(*,'("Dist  : ",a1)')dist
     write(*,'("Type  : ",a1)')type
     write(*,'("Mode  : ",i1)')mode
     write(*,'("Anorm : ",es10.3)')anorm
     write(*,'("Cond  : ",es10.3)')cndnum
     write(*,'(" ")')
      
     allocate(a(m,n), work(3*max(m,n)), d(n), stat=info)
     if(info .ne. 0) return
    
     call dlatms( m, n, dist, iseed, type, &
          & d, mode, cndnum, anorm, kl, ku, &
          & pack, a, m, work, info )
     if(info .ne. 0) return
     
     deallocate(work, stat=info)

     allocate(aux(n+2))

     call mergesortd(n, d, aux, -1)
     call mergeswapd(n, aux, d)

  case default
     write(*,'("Wrong matrix id")')
     info = 1
  end select

  write(*,'("Matrix generated")')

  return
end subroutine matgen
   


subroutine mergeswapd(n, l, d1)

  integer     :: i, lp, iswap, n
  integer     :: l(0:n+1)
  real(kind(1.d0)) :: d1(n), dswap

  lp = l(0)
  i  = 1
  do 
     if ((lp .eq. 0).or.(i>n)) exit
     do 
        if (lp >= i) exit
        lp = l(lp)
     end do
     dswap  = d1(lp)
     d1(lp) = d1(i)
     d1(i)  = dswap

     iswap  = l(lp)
     l(lp)  = l(i)
     l(i)   = lp
     lp     = iswap 
     i      = i + 1
  enddo

  return

end subroutine mergeswapd



subroutine mergesortd(n, k, l, order)

  !      Plain implementation of the merge-sort algorithm
  !      as described in:

  !      D. E. Knuth "The Art of Computer Programming,"
  !      vol.3: Sorting and Searching, Addison-Wesley, 1973
  implicit none
  integer                :: n
  integer, intent(inout) :: l(0:n+1)
  real(kind(1.d0))       :: k(n)
  integer                :: order

  integer    :: p, q, s, t, i, iord

  iord = order

  if((iord .ne. 1) .and. (iord .ne. -1)) then
     write(*,'("Wrong input in mergesort")')
     return
  end if

100 continue
  l(0) = 1
  t = n + 1
  do  p = 1,n - 1
     if (iord*k(p) .le. iord*k(p+1)) then
        l(p) = p + 1
     else
        l(t) = - (p+1)
        t = p
     end if
  end do
  l(t) = 0
  l(n) = 0

  if (l(n+1) .eq. 0) then
     return 
  else
     l(n+1) = iabs(l(n+1))
  end if

200 continue
  s = 0
  t = n+1
  p = l(s)
  q = l(t)

  if(q .eq. 0) return

300 continue
  if(iord*k(p) .gt. iord*k(q)) goto 600 

400 continue
  l(s) = sign(p,l(s))
  s = p
  p = l(p)
  if (p .gt. 0) goto 300

500 continue
  l(s) = q
  s = t
  do
     t = q
     q = l(q)
     if (q .le. 0) exit
  end do
  goto 800

600 continue
  l(s) = sign(q, l(s))
  s = q
  q = l(q)
  if (q .gt. 0) goto 300

700 continue
  l(s) = p
  s = t
  do
     t = p
     p = l(p)
     if (p .le. 0) exit
  end do

800 continue
  p = -p
  q = -q
  if(q.eq.0) then
     l(s) = sign(p, l(s))
     l(t) = 0
     goto 200
  end if

  goto 300

  return

end subroutine mergesortd
