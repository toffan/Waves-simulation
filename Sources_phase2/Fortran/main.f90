! This file is provided as part of the "projet long" for the Algebre Lineaire Numerique course
! at ENSEEIHT
! Date: 04/04/2012
! Authors: P. Amestoy, P. Berger, A. Buttari, Y. Diouane, S. Gratton, F.H. Rouet, E. Simon
!
! This file contains the driver for testing the developed methods
! 
program main
  ! use m_subspace_iter
  implicit none
  integer                       :: ierr, m, n, l, p, maxit, lwork, i, j, it, n_sv, disp, version, imat
  real(kind(1.d0)), allocatable :: a(:,:), v(:,:), u(:,:), work(:), b(:,:), s_vals(:)
  real(kind(1.d0)), allocatable :: s(:), res_sv(:)
  integer, allocatable          :: it_sv(:)
  integer                       :: t_start, t_end, t_rate, len
  real(kind(1.d0))              :: time, eps, res=0d0, percentage, percentage_ev, dummy
  real(kind(1.d0)), parameter   :: done=1.d0, dzero=0.d0

  interface
     subroutine matgen(imat, a, m, n, d, info)
       real(kind(1.d0)), allocatable :: a(:,:)
       integer :: imat, m, n, info
       real(kind(1.d0)), allocatable :: d(:)
     end subroutine matgen
  end interface


  !! parameter initializations 
  ierr=0

  maxit=1000
  eps=1d-13
  percentage=10d-3

  !! specify the dimensions n, p and the number of the dominante eigen vectors needed
  write(*,'("Version (0 :symmetric or 1 :rectangular?")')
  read(*,*) version
  write(*,'("p (block-product size)?")')
  read(*,*) p
  write(*,'("m (number of rows of the matrix a)?")')
  read(*,*) m
  write(*,'("n (number of columns of the matrix a)?")')
  read(*,*) n
  write(*,'("imat (matrix generator: 1, 2, 3 or 4)?")')
  read(*,*) imat
  write(*,'("l ( dimension of the invariant subspace )? ")')
  read(*,*) l
  write(*,'("Percentage (the percentage of the trace needed ) ")')
  read(*,*) percentage
  write(*,'("Do you want to display the singular values?(1 :yes or 0 : no) ")')
  read(*,*) disp

  if(n.gt.m) then
     write(*,'("m should be >= n. Aborting")')
     stop
  end if
  
  ! Generate the matrix
  call matgen(imat, a, m, n, s_vals, ierr)

  ! call random_number(a)
  
  if(disp.eq.1) then
     write(*,'("Here are the first ",i2," singular values: ")', advance='no')min(n,10)
     do i=1, min(n,10)
        write(*,'(f10.3,x)',advance='no')s_vals(i)
     end do
     write(*,'(" ")')
  end if

  write(*,'(" ")')
  write(*,'(" ")')
  select case(version)
  case(0)
     write(*,'("Calling the code for symmetric matrices")')
     allocate(b(m,m), s(l), res_sv(l), it_sv(l), u(m,l), stat=ierr)
     if(ierr .ne.0) then
        write(*,'("Can not allocate data, Input matrix is too big for the symmetric version")')
        stop
     endif

     !! initialize the subspace v0
     call random_number(u)
     s = 0.d0

     call system_clock(t_start)
     percentage_ev = 2*percentage -percentage*percentage
     ! Build B=A*A'
     call dgemm('n', 't', m, m, n, done, a, m, a, m, dzero, b, m )
     call subspace_iter_ev(b, u, s, m, l, p, percentage_ev, maxit, eps, res_sv, it_sv, it, n_sv, ierr)
     call system_clock(t_end, t_rate)
     time = real(t_end-t_start)/real(t_rate)
     deallocate(b, stat=ierr)
     if(ierr.eq.0) then
        do i=1, n_sv
           s(i) = sqrt(s(i))
        end do
     end if
  case(1)
     write(*,'("Calling the code for unsymmetric matrices")')
     allocate(s(l), res_sv(l), it_sv(l), u(m,l), stat=ierr)
     if(ierr .ne.0) then
        write(*,'("Can not allocate data.")')
        stop
     endif

     !! initialize the subspace v0
     call random_number(u)
     s = 0.d0
     
     call system_clock(t_start)
     call subspace_iter_sv(a, u, s, m, n, l, p, percentage, maxit, eps, res_sv, it_sv, it, n_sv, ierr)
     call system_clock(t_end, t_rate)
     time = real(t_end-t_start)/real(t_rate)
     
  case(2)
     write(*,'("Calling LAPACK DGESVD")')
     ! LAPACK, for testing
     lwork=2*max(1,3*min(m,n)+max(m,n),5*min(m,n))
     allocate(s(n), work(lwork), u(m,n), res_sv(n), it_sv(n), stat=ierr)
     if(ierr .ne.0) then
        write(*,'("Can not allocate data, Input matrix is too big for the LAPACK version")')
        stop
     endif
     call system_clock(t_start)
     call dgesvd('S', 'N', m, n, a, m, s, u, m, dummy, m, work, lwork, ierr)
     call system_clock(t_end, t_rate)
     time = real(t_end-t_start)/real(t_rate)
     n_sv   = n
     it_sv  = 1
     res_sv = 0.d0
  case default
     write(*,'("Version is not valid")')
     
  end select
  
  write(*,'(" ")')
  write(*,'(" ")')
  if(ierr.eq.0) then
     write(*,'("Done in ",f5.2,"  seconds")')time
     if(disp.eq.1) then
        do i=1, min(n_sv,10)
           write(*,'("Singular value ",i4,":",f10.3,"  with a residual:",es8.1," which converged at it:",i4,".    %",f4.0)')&
                & i,s(i),res_sv(i),it_sv(i),100.d0*s(i)/s(1)
        end do
        write(*,'(" ")')
     end if
  else
     write(*,'("An error was found. IERR=",i4)')ierr
  end if
     
  deallocate(s, u, res_sv, it_sv, stat=ierr)
  if(allocated(work)) deallocate(work, stat=ierr)

  stop

end program main
