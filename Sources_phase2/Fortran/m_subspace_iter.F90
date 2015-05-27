! This file is provided as part of the "projet long" for the Algebre Lineaire Numerique course
! at ENSEEIHT
! Date: 04/2015
! Authors: P. Amestoy, P. Berger, A. Buttari, Y. Diouane, S. Gratton, F.H. Rouet, E. Simon
!
! This file contains the implementation of
! - one routine for computing the eigenvalues of a symmetric matrix
! - one routine for computing the singular values of an unsymmetric matrix
!
! both are based on the subspace iteration method
!
! module m_subspace_iter
  ! implicit none
  ! private
  ! public :: subspace_iter_ev, subspace_iter_sv
! contains




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Subspace iteration with Rayleigh-Ritz.
! In this case the convergence can be checked for each eigenvector separately
! and the method can be stopped when the convereged eigenvectors capture
! as much as the trace of A as requested by the "percentage" argument
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine subspace_iter_ev(a, v, d, n, l, p, percentage, maxit, eps,  &
       res_ev, it_ev, it, n_ev, ierr)
    implicit none
    !! the matrix and subspace dimensions
    integer,          intent(in)                            :: n, l
    !! the number of products per iteration
    integer,          intent(in)                            :: p
    !! the traget matrix
    double precision, dimension(n,n), intent(in)            :: a
    !! maximum # of iteration
    integer,          intent(in)                            :: maxit
    !! the number of the dominant eigenvectors to compute
    double precision, intent(in)                            :: percentage
    !! the tolerance for the stopping criterion
    double precision, intent(in)                            :: eps
    !! the starting subspace. The computed eigenvectors will be
    !! returned in this array
    double precision, dimension(n,l), intent(inout)         :: v
    !! the residuals for each each eigenvector
    double precision, dimension(l),   intent(out)           :: res_ev
    !! the n_ev dominant eigenvalues
    double precision, dimension(l),   intent(out)           :: d
    !! the number of iteration to converge for each eigenvector
    integer,          dimension(l),   intent(out)           :: it_ev
    !! the global number iteration to converge
    integer,          intent(out)                           :: it
    !! the number of converged eigenvectors
    integer,          intent(out)                           :: n_ev
    !! a flag for signaling errors
    integer,          intent(out)                           :: ierr

    !! local variables
    integer                                                 :: i, j, conv, lwork
    double precision, allocatable, dimension(:,:)           :: y
    double precision, allocatable, dimension(:,:)           :: h
    double precision, allocatable, dimension(:)             :: w
    double precision                                        :: lambda, res, condition, tmp
    integer, parameter                                      :: ione=1

    ! external functions
    double precision, external                              :: dlange, dnrm2

#if defined(mex)
    character(len=80) :: string
    integer*4, external :: mexPrintf
    integer*4 :: k
#endif

    ierr = 0
    !! initialization process ...
    if((percentage.gt.1d0)  .or. (percentage.lt.0d0)) then
       ierr=1
       return
    end if
    if(l.gt.n) then
       ierr=1
       return
    end if

    allocate(y(n,l), h(l,l), w(l), stat=ierr)
    if(ierr .ne. 0) return

    lwork = n*l
    
    ! a natural choice of lambda is an estimate of some norm of a.
    ! use y as a workspace
    lambda=dlange('f', n, n, a, n, y)

    call orth_basis(v, n, l, y)

    n_ev      = 0
    v         = y
    res_ev    = 0d0
    d         = 0d0
    it_ev     = 0
    condition = 0
    
    it=0
    do while( (n_ev.lt.l) .and. (it.lt.maxit) )

       it = it+1
       !!compute  y=a*v
       do i=1, p
          call dgemm('n', 'n', n, l, n, 1.d0, a, n, v, n, 0.d0, y, n)
          v=y
       end do

       !! orthogonalization using Gramm Schmidt procedure
       !! compute  v=orth(y)
       call orth_basis(y, n, l, v)

       !!compute  y=a*v
       call dgemm('n','n', n, l, n, 1d0, a, n, v, n, 0d0, y, n)
       !!compute  h=v'*y
       call dgemm('t', 'n', l, l, n, 1d0, v, n, y, n, 0d0, h, l)
       ! Compute eig-decomposition of h. Use y as a workspace
       call dsyev('v', 'u', l, h, l, w, y, lwork, ierr)
       if( ierr .ne.0 )then
          goto 9999
       end if
       !! Sort in decreasing order (in-place)
       !!(we suppose that all the eigen values are positive)
       do i=1,l/2
          tmp        = w(i)
          w(i)       = w(l-i +1)
          w(l-i+1)   = tmp
          ! use the first column of y as a workspace
          y(1:l,1)   = h(:,i)
          h(:,i)     = h(:,l-i +1)
          h(:,l-i+1) = y(1:l,1)
       end do

       !! v=v*h
       y = v
       call dgemm('n', 'n', n, l, l, 1d0, y, n, h, l, 0d0, v, n)
       conv = 0
       i    = n_ev +1
       !! the larger eigenvalue will converge more swiftly than
       !! those corresponding to the smaller eigenvalue.
       !! for this reason, we test the convergence in the order
       !! i=1,2,.. and stop with the first one to fail the test
       do
          if( i .gt. l) exit
          !! compute res=norm(a*v(:,i) - v(:,i)*t(i,i),2)/lambda;
          !! use the first column of y as a workspace
          y(:,1) = v(:,i)
          call dgemv('n', n, n, 1d0, a, n, v(1,i), 1, -w(i), y, 1)
          res = dnrm2(n, y, ione)/lambda
          if(res.gt.eps) exit
          conv         = conv+1
          d(i)         = w(i)
          res_ev(i)    = res
          it_ev(i)     = it
          i            = i+1
          condition    = 1.d0-d(n_ev+conv)/d(1)
       end do

       n_ev = n_ev + conv

#if defined (mex)
       write(string,'(" IT:",i5," -- Found ",i4," eigenvalues (",f8.5,"% achieved)")')it,n_ev,condition
       k = mexprintf(string//achar(10))
#else
       write(*,'(" IT:",i5," -- Found ",i4," eigenvalues (",f8.5,"% achieved)",a)',advance='no')it,n_ev,condition,char(13)
       write(*,'(" ")')
#endif
       if( condition .gt. percentage) exit
       if( n_ev .ge. l) exit

    end do


    if(it .ge. maxit)then
       ierr=4
       goto 9999
    end if

9999 continue
    deallocate(y, h, w)

    return

  end subroutine subspace_iter_ev

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine subspace_iter_sv(a, u, s, m, n, l, p, percentage, maxit, eps,   &
       res_sv, it_sv, it, n_sv, ierr)
    implicit none
    !! the subspace dimensions
    integer,          intent(in)                            :: m, n, l
    !! the number of products per iteration
    integer,          intent(in)                            :: p
    !! the target matrix
    double precision, dimension(m,n), intent(in)            :: a
    !! maximum # of iteration
    integer,          intent(in)                            :: maxit
    !! the number of the dominant eigenvectors to compute
    double precision, intent(in)                            :: percentage
    !! the tolerance for the stopping criterion
    double precision, intent(in)                            :: eps
    !! the starting subspace. The computed eigenvectors will be
    !! returned in this array
    double precision, dimension(m,l), intent(inout)         :: u
    !! the residuals for each each eigenvector
    double precision, dimension(l),   intent(out)           :: res_sv
    !! the n_sv dominant eigenvalues
    double precision, dimension(l),   intent(out)           :: s
    !! the number of iteration to converge for each eigenvector
    integer,          dimension(l),   intent(out)           :: it_sv
    !! the global number iteration to converge
    integer,          intent(out)                           :: it
    !! the number of converged eigenvectors
    integer,          intent(out)                           :: n_sv
    !! a flag for signaling errors
    integer,          intent(out)                           :: ierr

    !! local variables
    integer                                                 :: i, j, conv, lwork
    double precision, allocatable, dimension(:,:)           :: y, z
    double precision, allocatable, dimension(:,:)           :: h
    double precision, allocatable, dimension(:)             :: w
    double precision                                        :: lambda, res, condition, tmp
    integer, parameter                                      :: ione=1

    ! external functions
    double precision, external                              :: dlange, dnrm2

#if defined(mex)
    character(len=80) :: string
    integer*4, external :: mexPrintf
    integer*4 :: k
#endif

   
    ierr = 0
    !! initialization process ...
    if((percentage.gt.1d0)  .or. (percentage.lt.0d0)) then
       ierr=1
       return
    end if
    if(l.gt.n) then
       ierr=1
       return
    end if

    allocate(y(m,l), z(n,l), h(l,l), w(l), stat=ierr)
    if(ierr .ne. 0) return
 
    lwork = n*l
    
    ! a natural choice of lambda is an estimate of some norm of a.
    ! use y as a workspace
    lambda=dlange('f', m, n, a, m, y)

    call orth_basis(u, n, l, y)

    res_sv    = 0d0
    it_sv     = 0d0
    s         = 0d0
    it        = 0
    n_sv      = 0
    
    it=0
    do while( (n_sv.lt.l) .and. (it.lt.maxit) )

       it = it+1
       !!compute  y=a*v
       do i=1, p
          call dgemm('t', 'n', n, l, m, 1.d0, a, m, u, m, 0.d0, z, n)
          call dgemm('n', 'n', m, l, n, 1.d0, a, m, z, n, 0.d0, y, m)
          u=y
       end do

       !! orthogonalization using Gramm Schmidt procedure
      !! compute  v=orth(y)
       call orth_basis(y, m, l, u)

       !!compute  y=a'*v
       call dgemm('t','n', n, l, m, 1d0, a, m, u, m, 0d0, z, n)
       !!compute  h=y'*y
       call dgemm('t', 'n', l, l, n, 1d0, z, n, z, n, 0d0, h, l)
       ! Compute eig-decomposition of h. Use y as a workspace
       call dsyev('v', 'u', l, h, l, w, y, lwork, ierr)
       if( ierr .ne.0 )then
           goto 9999
       end if
       !! Sort in decreasing order (in-place)
       !!(we suppose that all the eigen values are positive)
       do i=1,l/2
          tmp        = w(i)
          w(i)       = w(l-i +1)
          w(l-i+1)   = tmp
          ! use the first column of y as a workspace
         y(1:l,1)   = h(:,i)
          h(:,i)     = h(:,l-i +1)
          h(:,l-i+1) = y(1:l,1)
       end do

       !! v=v*h
       y = u
       call dgemm('n', 'n', m, l, l, 1d0, y, m, h, l, 0d0, u, m)
       conv = 0
       i    = n_sv +1
       !! the larger eigenvalue will converge more swiftly than
       !! those corresponding to the smaller eigenvalue.
       !! for this reason, we test the convergence in the order
       !! i=1,2,.. and stop with the first one to fail the test
       do
          if( i .gt. l) exit
          !! compute res=norm(a*v(:,i) - v(:,i)*t(i,i),2)/lambda;
          !! use the first column of y as a workspace
          y(:,1) = u(:,i)
          call dgemv('t', m, n, 1d0, a, m, u(1,i), 1, 0d0, y(1,2), 1)
          call dgemv('n', m, n, 1d0, a, m, y(1,2), 1, -w(i), y, 1)
          res = dnrm2(m, y, ione)/(lambda**2)
          print*, res, "!=", eps
          if(res.gt.eps) exit
          conv         = conv+1
          s(i)         = sqrt(w(i))
          res_sv(i)    = res
          it_sv(i)     = it
          i            = i+1
          condition    = 1.d0-s(n_sv+conv)/s(1)
          if(condition.gt.percentage) exit
       end do

       n_sv = n_sv + conv

#if defined (mex)
       write(string,'(" IT:",i5," -- Found ",i4," eigenvalues (",f8.5,"% achieved)")')it,n_sv,condition
       k = mexprintf(string//achar(10))
#else
       write(*,'(" IT:",i5," -- Found ",i4," eigenvalues (",f8.5,"% achieved)",a)',advance='no')it,n_sv,condition,char(13)
       write(*,'(" ")')
#endif
       if( condition .gt. percentage) exit
       if( n_sv .ge. l) exit

    end do


    if(it .ge. maxit)then
       ierr=4
       goto 9999
    end if

9999 continue
    deallocate(y, h, w)

    return


  end subroutine subspace_iter_sv



  !!========================================================================
  !! orthogonalization using gram shmidt procedure
  !!========================================================================
  subroutine orth_basis(u_in,n,m,u_out)
    implicit none
    integer          ,intent(in)                            :: n,m
    double precision,dimension(n,m),intent(in)              :: u_in
    double precision,dimension(n,m),intent(out)              :: u_out
    !! local variables
    integer ::i,j
    do i=1,m
       u_out(:,i)=u_in(:,i)
    end do
    do j=1,m
       if( j.gt. 1)then
          do i=1,j-1
             u_out(:,j)= u_out(:,j) - dot_product(u_out(:,j),              &
                  u_out(:,i))*u_out(:,i)
          end do
       end if
       u_out(:,j)=u_out(:,j)/sqrt(dot_product(u_out(:,j),u_out(:,j)))
    end do
  end subroutine orth_basis


! end module m_subspace_iter
