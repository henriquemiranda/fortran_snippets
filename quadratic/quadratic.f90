!
! 20/03/2019
! Author: Henrique Miranda
! Create a linear system of equations to fit a quadratic function in 3D
! the function has the general form:
!
!  f = a +
!      ax*x + ay*y + az*z +
!      axy*x*y + axz*x*z + ayz*y*z +
!      ax2*x*x + ayy*y*y + azz*z*z     (1)
!
! Use this as an interpolating function inside the tetrahedron
! and apply the hybrid tetrahedron method described in
!
! A.H. MacDonald, S.H. Vosko, and P.T. Coleridge,
! Journal of Physics C: Solid State Physics 12, 2991 (1979).
!
module m_defs
 integer, parameter :: dp=8
 real(dp), parameter :: zero = 0.0_dp, half = 0.5_dp, one = 1.0_dp, two=2.0_dp
 real(dp), parameter :: tol8  = 0.00000001_dp
 real(dp), parameter :: tol14 = 0.000000000000001_dp
 real(dp), parameter :: pi=3.141592653589793238462643383279502884197_dp
 real(dp), parameter :: two_pi = two*pi
end module m_defs

module m_quadratic
 use m_defs
 implicit none
 integer,parameter  :: lwork=200
 real(dp),parameter :: rcond=tol8
 real(dp) :: work(lwork)

 contains
 real(dp) function f(kpt,parm)
   ! Compute the quadratic function on kpt
   real(dp),intent(in) :: kpt(3)
   real(dp),intent(in) :: parm(10)
   real(dp) :: x,y,z
   x = kpt(1);y = kpt(2);z = kpt(3)
   !f = parm(1) + &
   !    parm(2)*x   + parm(3)*y   + parm(4)*z + &
   !    parm(5)*x*y + parm(6)*x*z + parm(7)*y*z + &
   !    parm(8)*x*x + parm(9)*y*y + parm(10)*z*z
   f = parm(1)+x*(parm(2)+y*parm(5)+z*parm(6)+x*parm(8))+&
               y*(parm(3)+z*parm(7)          +y*parm(9))+&
               z*(parm(4)                    +z*parm(10))
 end function f

 function df(kpt,parm)
   ! Compute the derivative of the quadratic function with respect to kx,ky,kz
   real(dp),intent(in) :: kpt(3)
   real(dp),intent(in) :: parm(10)
   real(dp) :: df(3)
   real(dp) :: x,y,z
   x = kpt(1);y = kpt(2);z = kpt(3)
   df(1) = parm(2) + &
           parm(5)*y + parm(6)*z + &
           2.0_dp*parm(8)*x
   df(2) = parm(3) + &
           parm(5)*x + parm(7)*z + &
           2.0_dp*parm(9)*y
   df(3) = parm(4) + &
           parm(6)*x + parm(7)*y + &
           2.0_dp*parm(10)*z
 end function df

 subroutine fit4(kpts,eigs,parm,info)
  ! Fit a linear function using the values of the function evaluated in 4 points
  ! This is done solving a linear system of equations of the form Ax=b
  real(dp),intent(in)  :: kpts(3,10)
  real(dp),intent(in)  :: eigs(4)
  real(dp),intent(out) :: parm(10)
  integer,intent(out)  :: info

  integer :: ii,irank
  real(dp) :: x,y,z
  real(dp) :: sgval(4)
  real(dp) :: a(4,4), b(4)

  ! Fill in the a matrix A
  do ii=1,4
    x=kpts(1,ii);y=kpts(2,ii);z=kpts(3,ii)
    a(ii,:) = [1.0_dp,x,y,z]
  end do

  ! Fill in the vector b
  b = eigs

  !solve a system of linear equations
  call dgelss(4,4,1,a,4,b,4,sgval,rcond,irank,work,lwork,info)
  parm(:4) = b(:)
  parm(5:) = zero
 end subroutine fit4

 subroutine fit10(kpts,eigs,parm,info)
  ! Fit a quadratic function using the values fo the function evaluated in 10 points
  ! This is done solving a linear system of equations of the form Ax=b
  real(dp),intent(in) :: kpts(3,10)
  real(dp),intent(in) :: eigs(10)
  real(dp),intent(out) :: parm(10)
  integer,intent(out) :: info

  integer :: ii,irank
  real(dp) :: x,y,z
  real(dp) :: sgval(10)
  real(dp) :: a(10,10), b(10)

  ! Fill in the a matrix A
  do ii=1,10
    x=kpts(1,ii);y=kpts(2,ii);z=kpts(3,ii)
    a(ii,:) = [1.0_dp,x,y,z,x*y,x*z,y*z,x*x,y*y,z*z]
  end do

  ! Fill in the vector b
  b = eigs

  !solve a system of linear equations
  call dgelss(10,10,1,a,10,b,10,sgval,rcond,irank,work,lwork,info)
  parm = b
 end subroutine fit10

 subroutine fit10_affine(eigs,parm)
  ! Fit a quadratic function using the values fo the function evaluated in 10 points
  ! This is done by inverting A of a linear system of equations of the form Ax=b
  ! A affine transformation of coordinates (x,y,z)->(a,b,c) is used following
  ! CMES: Computer Modeling in Engineering & Sciences, Vol. 94, No. 4, pp. 279-300, 2013
  ! (x1,y1,z1) -> (0,0,0)
  ! (x2,y2,z2) -> (1,0,0)
  ! (x3,y3,z3) -> (0,1,0)
  ! (x4,y4,z4) -> (0,0,1)
  !
  ! In this coordinates, the 10 vertices of the tetrahedra are:
  ! (  0,   0,   0)
  ! (  1,   0,   0)
  ! (  0,   1,   0)
  ! (  0,   0,   1)
  ! (1/2,   0,   0)
  ! (  0, 1/2,   0)
  ! (  0,   0  1/2)
  ! (1/2, 1/2,   0)
  ! (  0, 1/2, 1/2)
  ! (1/2,   0, 1/2)
  !
  ! We want to describe a quadratic function (1)
  ! whose coefficients sym are:
  ! (a,  ax,  ay,  az, axy, axz, ayz, ax2, ayy, azz)
  !
  ! A parm = eigs
  ! The matrix A is d f /d sym
  ! [1,   0,   0,   0,   0,   0,   0,   0,   0,   0]
  ! [1,   1,   0,   0,   0,   0,   0,   1,   0,   0]
  ! [1,   0,   1,   0,   0,   0,   0,   0,   1,   0]
  ! [1,   0,   0,   1,   0,   0,   0,   0,   0,   1]
  ! [1, 1/2,   0,   0,   0,   0,   0, 1/4,   0,   0]
  ! [1,   0, 1/2,   0,   0,   0,   0,   0, 1/4,   0]
  ! [1,   0,   0, 1/2,   0,   0,   0,   0,   0, 1/4]
  ! [1, 1/2, 1/2,   0, 1/4,   0,   0, 1/4, 1/4,   0]
  ! [1,   0, 1/2, 1/2,   0,   0, 1/4,   0, 1/4, 1/4]
  ! [1, 1/2,   0, 1/2,   0, 1/4,   0, 1/4,   0, 1/4]
  !
  real(dp),intent(in) :: eigs(10)
  real(dp),intent(out) :: parm(10)

  !real(dp) :: inva(10,10)

  ! Fill in the inverse of the A matrix
  !inva(:,1)  = [1,-3,-3,-3, 4, 4, 4, 2, 2, 2]
  !inva(:,2)  = [0,-1, 0, 0, 0, 0, 0, 2, 0, 0]
  !inva(:,3)  = [0, 0,-1, 0, 0, 0, 0, 0, 2, 0]
  !inva(:,4)  = [0, 0, 0,-1, 0, 0, 0, 0, 0, 2]
  !inva(:,5)  = [0, 4, 0, 0,-4,-4, 0,-4, 0, 0]
  !inva(:,6)  = [0, 0, 4, 0,-4, 0,-4, 0,-4, 0]
  !inva(:,7)  = [0, 0, 0, 4, 0,-4,-4, 0, 0,-4]
  !inva(:,8)  = [0, 0, 0, 0, 4, 0, 0, 0, 0, 0]
  !inva(:,9)  = [0, 0, 0, 0, 0, 0, 4, 0, 0, 0]
  !inva(:,10) = [0, 0, 0, 0, 0, 4, 0, 0, 0, 0]
  !call dgemv('N',10,10,one,inva,10,eigs,1,zero,parm,1)

  parm( 1) =  eigs(1)
  parm( 2) = -3.0_dp*eigs(1)-eigs(2)+4.0_dp*eigs(5)
  parm( 3) = -3.0_dp*eigs(1)-eigs(3)+4.0_dp*eigs(6)
  parm( 4) = -3.0_dp*eigs(1)-eigs(4)+4.0_dp*eigs(7)
  parm( 5) = +4.0_dp*(eigs(1)-eigs(5)-eigs(6)+eigs(8))
  parm( 6) = +4.0_dp*(eigs(1)-eigs(5)-eigs(7)+eigs(10))
  parm( 7) = +4.0_dp*(eigs(1)-eigs(6)-eigs(7)+eigs(9))
  parm( 8) = +2.0_dp*(eigs(1)+eigs(2)-2.0_dp*eigs(5))
  parm( 9) = +2.0_dp*(eigs(1)+eigs(3)-2.0_dp*eigs(6))
  parm(10) = +2.0_dp*(eigs(1)+eigs(4)-2.0_dp*eigs(7))

 end subroutine fit10_affine

 subroutine fit4d(kpts,eigs,vels,parm,info)
  ! Fit a quadratic function using the values and derivatives of the function evaluated in 4 points
  ! This is done solving a linear system of equations of the form Ax=b
  real(dp),intent(in) :: kpts(3,4)
  real(dp),intent(in) :: eigs(4)
  real(dp),intent(in) :: vels(3,4)
  real(dp),intent(out) :: parm(10)
  integer,intent(out) :: info

  integer :: ii
  integer :: irank
  real(dp) :: x,y,z
  real(dp) :: sgval(10)
  real(dp) :: a(16,10),b(16)

  ! Fill in the a matrix for eigenvalues
  do ii=1,4
    x=kpts(1,ii);y=kpts(2,ii);z=kpts(3,ii)
    a(ii,:) = [1.0_dp,x,y,z,x*y,x*z,y*z,x*x,y*y,z*z]
  end do

  ! Fill in the matrix for velocities (explicit just because)
  x=kpts(1,1);y=kpts(2,1);z=kpts(3,1)
  a(5,:)  = [0.0_dp,1.0_dp,0.0_dp,0.0_dp,y,z,0.0_dp,2.0_dp*x,0.0_dp,0.0_dp] !dx
  a(6,:)  = [0.0_dp,0.0_dp,1.0_dp,0.0_dp,x,0.0_dp,z,0.0_dp,2.0_dp*y,0.0_dp] !dy
  a(7,:)  = [0.0_dp,0.0_dp,0.0_dp,1.0_dp,0.0_dp,x,y,0.0_dp,0.0_dp,2.0_dp*z] !dz
  x=kpts(1,2);y=kpts(2,2);z=kpts(3,2)
  a(8,:)  = [0.0_dp,1.0_dp,0.0_dp,0.0_dp,y,z,0.0_dp,2.0_dp*x,0.0_dp,0.0_dp] !dx
  a(9,:)  = [0.0_dp,0.0_dp,1.0_dp,0.0_dp,x,0.0_dp,z,0.0_dp,2.0_dp*y,0.0_dp] !dy
  a(10,:) = [0.0_dp,0.0_dp,0.0_dp,1.0_dp,0.0_dp,x,y,0.0_dp,0.0_dp,2.0_dp*z] !dz
  x=kpts(1,3);y=kpts(2,3);z=kpts(3,3)
  a(11,:) = [0.0_dp,1.0_dp,0.0_dp,0.0_dp,y,z,0.0_dp,2.0_dp*x,0.0_dp,0.0_dp] !dx
  a(12,:) = [0.0_dp,0.0_dp,1.0_dp,0.0_dp,x,0.0_dp,z,0.0_dp,2.0_dp*y,0.0_dp] !dy
  a(13,:) = [0.0_dp,0.0_dp,0.0_dp,1.0_dp,0.0_dp,x,y,0.0_dp,0.0_dp,2.0_dp*z] !dz
  x=kpts(1,4);y=kpts(2,4);z=kpts(3,4)
  a(14,:) = [0.0_dp,1.0_dp,0.0_dp,0.0_dp,y,z,0.0_dp,2.0_dp*x,0.0_dp,0.0_dp] !dx
  a(15,:) = [0.0_dp,0.0_dp,1.0_dp,0.0_dp,x,0.0_dp,z,0.0_dp,2.0_dp*y,0.0_dp] !dy
  a(16,:) = [0.0_dp,0.0_dp,0.0_dp,1.0_dp,0.0_dp,x,y,0.0_dp,0.0_dp,2.0_dp*z] !dz

  ! Fill in the vector b
  b(:4)    = eigs
  b(5:7)   = vels(:,1)
  b(8:10)  = vels(:,2)
  b(11:13) = vels(:,3)
  b(14:16) = vels(:,4)

  !solve a system of linear equations
  call dgelss(16,10,1,a,16,b,16,sgval,rcond,irank,work,lwork,info)
  parm = b(:10)
 end subroutine fit4d

end module m_quadratic

module m_tetralite

  use m_defs
  implicit none
  integer :: tetra_6lshifts(3,4,6)

  private :: f
  contains
  subroutine tetralite_init()
    tetra_6lshifts(:,1,1) = [0,0,0]
    tetra_6lshifts(:,2,1) = [1,0,0]
    tetra_6lshifts(:,3,1) = [0,1,0]
    tetra_6lshifts(:,4,1) = [1,0,1]
    tetra_6lshifts(:,1,2) = [1,0,0]
    tetra_6lshifts(:,2,2) = [1,1,0]
    tetra_6lshifts(:,3,2) = [0,1,0]
    tetra_6lshifts(:,4,2) = [1,0,1]
    tetra_6lshifts(:,1,3) = [0,1,0]
    tetra_6lshifts(:,2,3) = [1,1,0]
    tetra_6lshifts(:,3,3) = [1,0,1]
    tetra_6lshifts(:,4,3) = [1,1,1]
    tetra_6lshifts(:,1,4) = [0,0,0]
    tetra_6lshifts(:,2,4) = [0,1,0]
    tetra_6lshifts(:,3,4) = [0,0,1]
    tetra_6lshifts(:,4,4) = [1,0,1]
    tetra_6lshifts(:,1,5) = [0,0,1]
    tetra_6lshifts(:,2,5) = [1,0,1]
    tetra_6lshifts(:,3,5) = [0,1,0]
    tetra_6lshifts(:,4,5) = [0,1,1]
    tetra_6lshifts(:,1,6) = [0,1,0]
    tetra_6lshifts(:,2,6) = [1,0,1]
    tetra_6lshifts(:,3,6) = [0,1,1]
    tetra_6lshifts(:,4,6) = [1,1,1]
  end subroutine

  pure function grid(divs,gp) result(ik)
    integer,intent(in) :: divs(3),gp(3)
    integer :: ik
    integer :: ix,iy,iz
    ix=gp(1); iy=gp(2); iz=gp(3)
    ix = mod(ix-1,divs(1))+1
    iy = mod(iy-1,divs(2))+1
    iz = mod(iz-1,divs(3))+1
    ik = (ix-1)*divs(2)*divs(3)+(iy-1)*divs(3)+iz
  end function

  pure function gridp(divs,gp) result(ik)
    integer,intent(in) :: divs(3),gp(3)
    integer :: ik
    integer :: ix,iy,iz
    ix=gp(1); iy=gp(2); iz=gp(3)
    ik = (ix-1)*divs(2)*divs(3)+(iy-1)*divs(3)+iz
  end function

  pure function inv_grid(divs,ik) result(gp)
    integer,intent(in) :: divs(3),ik
    integer :: gp(3)
    integer :: ix,iy,iz,lda,ldb
    lda = divs(2)*divs(3)
    ldb = divs(3)
    ix=(ik-1)/lda
    iy=(ik-ix*lda-1)/ldb
    iz=(ik-ix*lda-iy*ldb-1)
    gp = [ix+1,iy+1,iz+1]
  end function

  pure function grid_kpoint(divs,gp) result(kpt)
    ! Get the kpoint coordinates from gridpoint
    integer,intent(in) :: divs(3),gp(3)
    real(dp) :: kpt(3)
    kpt = one*(gp-1)/divs-half
  end function

  pure real(dp) function f(n,m,w,eig)
    integer,intent(in) :: n,m
    real(dp),intent(in) :: w,eig(4)
    f = ((w - eig(m)) / (eig(n) - eig(m)))
  end function

  pure real(dp) function g1(w,eig)
    real(dp),intent(in) :: w,eig(4)
    g1 = (3.0_dp *&
          f(2, 1, w, eig) *&
          f(3, 1, w, eig) /&
          (eig(4) - eig(1)))
  end function

  pure real(dp) function g2(w,eig)
    real(dp),intent(in) :: w,eig(4)
    g2 = (3.0_dp /&
          (eig(4) - eig(1)) *&
          (f(2, 3, w, eig) *&
           f(3, 1, w, eig) +&
           f(3, 2, w, eig) *&
           f(2, 4, w, eig)))
  end function

  pure real(dp) function g3(w,eig)
    real(dp),intent(in) :: w,eig(4)
    g3 = (3.0_dp *&
            f(2, 4, w, eig) *&
            f(3, 4, w, eig) /&
            (eig(4) - eig(1)))
  end function

  pure real(dp) function I10(w,eig)
    real(dp),intent(in) :: w,eig(4)
    I10 = (f(1, 2, w, eig) +&
           f(1, 3, w, eig) +&
           f(1, 4, w, eig))
  end function

  pure real(dp) function I11(w,eig)
    real(dp),intent(in) :: w,eig(4)
    I11 = f(2, 1, w, eig)
  end function

  pure real(dp) function I12(w,eig)
    real(dp),intent(in) :: w,eig(4)
    I12 = f(3, 1, w, eig)
  end function

  pure real(dp) function I13(w,eig)
    real(dp),intent(in) :: w,eig(4)
    I13 = f(4, 1, w, eig)
  end function

  pure real(dp) function I20(w,eig)
    real(dp),intent(in) :: w,eig(4)
    I20 = (f(1, 4, w, eig) +&
           f(1, 3, w, eig) *&
           f(3, 1, w, eig) *&
           f(2, 3, w, eig) /&
          (f(2, 3, w, eig) *&
           f(3, 1, w, eig) +&
           f(3, 2, w, eig) *&
           f(2, 4, w, eig)))
  end function

  pure real(dp) function I21(w,eig)
    real(dp),intent(in) :: w,eig(4)
    I21 = (f(2, 3, w, eig) +&
           f(2, 4, w, eig) *&
           f(2, 4, w, eig) *&
           f(3, 2, w, eig) /&
          (f(2, 3, w, eig) *&
           f(3, 1, w, eig) +&
           f(3, 2, w, eig) *&
           f(2, 4, w, eig)))
  end function

  pure real(dp) function I22(w,eig)
    real(dp),intent(in) :: w,eig(4)
    I22 = (f(3, 2, w, eig) +&
           f(3, 1, w, eig) *&
           f(3, 1, w, eig) *&
           f(2, 3, w, eig) /&
          (f(2, 3, w, eig) *&
           f(3, 1, w, eig) +&
           f(3, 2, w, eig) *&
           f(2, 4, w, eig)))
  end function

  pure real(dp) function I23(w,eig)
    real(dp),intent(in) :: w,eig(4)
    I23 = (f(4, 1, w, eig) +&
           f(4, 2, w, eig) *&
           f(2, 4, w, eig) *&
           f(3, 2, w, eig) /&
          (f(2, 3, w, eig) *&
           f(3, 1, w, eig) +&
           f(3, 2, w, eig) *&
           f(2, 4, w, eig)))
  end function

  pure real(dp) function I30(w,eig)
    real(dp),intent(in) :: w,eig(4)
    I30 = f(1, 4, w, eig)
  end function

  pure real(dp) function I31(w,eig)
    real(dp),intent(in) :: w,eig(4)
    I31 = f(2, 4, w, eig)
  end function

  pure real(dp) function I32(w,eig)
    real(dp),intent(in) :: w,eig(4)
    I32 = f(3, 4, w, eig)
  end function

  pure real(dp) function I33(w,eig)
    real(dp),intent(in) :: w,eig(4)
    I33 = (f(4, 1, w, eig) +&
           f(4, 2, w, eig) +&
           f(4, 3, w, eig))
  end function

  subroutine tetralite_delta(eig,wvals,nw,alpha,integral)
    ! function to evaluate delta
    ! theare are the expressions from
    ! A.H. MacDonald, S.H. Vosko, and P.T. Coleridge,
    ! Journal of Physics C: Solid State Physics 12, 2991 (1979).
    ! as implemented in spglib implemented by Atsushi Togo translated into fortran
    integer,intent(in) :: nw
    real(dp),intent(in) :: eig(4),alpha(4)
    real(dp),intent(in) :: wvals(nw)
    real(dp),intent(out) :: integral(nw)
    integer :: iw
    real(dp) :: w,g

    do iw=1,nw
      w = wvals(iw)

      ! w < e1 nothing to do
      if (w < eig(1)) cycle

      ! e1 < w < e2
      if (w < eig(2)) then
        g = g1(w,eig)/3.0_dp
        integral(iw) = integral(iw)+g*I10(w,eig)*alpha(1)
        integral(iw) = integral(iw)+g*I11(w,eig)*alpha(2)
        integral(iw) = integral(iw)+g*I12(w,eig)*alpha(3)
        integral(iw) = integral(iw)+g*I13(w,eig)*alpha(4)
        cycle
      endif

      ! e2 < eps < e3
      if (w < eig(3)) then
        g = g2(w,eig)/3.0_dp
        integral(iw) = integral(iw)+g*I20(w,eig)*alpha(1)
        integral(iw) = integral(iw)+g*I21(w,eig)*alpha(2)
        integral(iw) = integral(iw)+g*I22(w,eig)*alpha(3)
        integral(iw) = integral(iw)+g*I23(w,eig)*alpha(4)
        cycle
      endif

      ! e3 < eps
      if (w < eig(4)) then
        g = g3(w,eig)/3.0_dp
        integral(iw) = integral(iw)+g*I30(w,eig)*alpha(1)
        integral(iw) = integral(iw)+g*I31(w,eig)*alpha(2)
        integral(iw) = integral(iw)+g*I32(w,eig)*alpha(3)
        integral(iw) = integral(iw)+g*I33(w,eig)*alpha(4)
        cycle
      endif

      ! e4 < eps
      if (eig(4) < w) exit

    end do

  end subroutine tetralite_delta

  pure subroutine sort_4tetra(list,perm)

   integer,  intent(inout) :: perm(4)
   real(dp), intent(inout) :: list(4)

   integer :: ia,ib,ic,id
   integer :: ilow1,ilow2,ihigh1,ihigh2
   integer :: ilowest,ihighest
   integer :: imiddle1,imiddle2
   real(dp) :: va,vb,vc,vd
   real(dp) :: vlow1,vlow2,vhigh1,vhigh2
   real(dp) :: vlowest,vhighest
   real(dp) :: vmiddle1,vmiddle2

   va = list(1); ia = perm(1)
   vb = list(2); ib = perm(2)
   vc = list(3); ic = perm(3)
   vd = list(4); id = perm(4)

   if (va < vb) then
       vlow1 = va; vhigh1 = vb
       ilow1 = ia; ihigh1 = ib
   else
       vlow1 = vb; vhigh1 = va
       ilow1 = ib; ihigh1 = ia
   endif

   if (vc < vd) then
       vlow2 = vc; vhigh2 = vd
       ilow2 = ic; ihigh2 = id
   else
       vlow2 = vd; vhigh2 = vc
       ilow2 = id; ihigh2 = ic
   endif

   if (vlow1 < vlow2) then
       vlowest  = vlow1; vmiddle1 = vlow2
       ilowest  = ilow1; imiddle1 = ilow2
   else
       vlowest  = vlow2; vmiddle1 = vlow1
       ilowest  = ilow2; imiddle1 = ilow1
   endif

   if (vhigh1 > vhigh2) then
       vhighest = vhigh1; vmiddle2 = vhigh2
       ihighest = ihigh1; imiddle2 = ihigh2
   else
       vhighest = vhigh2; vmiddle2 = vhigh1
       ihighest = ihigh2; imiddle2 = ihigh1
   endif

   if (vmiddle1 < vmiddle2) then
       list = [vlowest,vmiddle1,vmiddle2,vhighest]
       perm = [ilowest,imiddle1,imiddle2,ihighest]
   else
       list = [vlowest,vmiddle2,vmiddle1,vhighest]
       perm = [ilowest,imiddle2,imiddle1,ihighest]
   endif

  end subroutine sort_4tetra

end module

module m_dividetetra
  ! This is to implement an quadratic hybrid tetrahedra integration method as deribed in:
  ! A.H. MacDonald, S.H. Vosko, and P.T. Coleridge,
  ! Journal of Physics C: Solid State Physics 12, 2991 (1979).
  use m_quadratic
  use m_tetralite
  implicit none
  integer :: ieigs(2,10)
  integer :: tetra8(4,8)
  integer :: tetra_6qshifts(3,10,6)

  contains
  subroutine dividetetra_init()
    ! Initialization routine
    integer :: itetra

    ! Indexes of the middlepoints w.r.t initial tetrahedron
    ieigs(:,1) = [1,2] !5
    ieigs(:,2) = [1,3] !6
    ieigs(:,3) = [1,4] !7
    ieigs(:,4) = [2,3] !8
    ieigs(:,5) = [3,4] !9
    ieigs(:,6) = [2,4] !10

    ! 8 tetrahedra from the paper
    tetra8(:,1) = [1,5, 6, 7]
    tetra8(:,2) = [2,5, 8,10]
    tetra8(:,3) = [3,8, 9, 6]
    tetra8(:,4) = [4,9, 7,10]

    tetra8(:,5) = [7,8,10, 9]
    tetra8(:,6) = [7,8,10, 5]
    tetra8(:,7) = [7,8, 6, 9]
    tetra8(:,8) = [7,8, 6, 5]

    ! Now I generate shifts for quadratic tetrahedron
    ! this is done by scaling the original shifts by 2
    ! and computing the 10 midpoints for each of the 6 tetrahedra
    do itetra=1,6
      ! the first 4 shifts are just a copy
      tetra_6qshifts(:,:4,itetra) = 2*tetra_6lshifts(:,:,itetra)
      ! the additional 6 are the midpoints
      tetra_6qshifts(:,5:,itetra) = nint(get_6midkpts(two*tetra_6lshifts(:,:,itetra)))
    end do

  end subroutine dividetetra_init

  function get_6midkpts(kpts) result(mid)
    real(dp),intent(in) :: kpts(3,4)
    real(dp) :: mid(3,6)
    integer :: ii,ik1,ik2
    ! get the 6 midpoints
    do ii=1,6
     ik1 = ieigs(1,ii)
     ik2 = ieigs(2,ii)
     mid(:,ii) = half*(kpts(:,ik1)+kpts(:,ik2))
    end do
  end function

  function get_values(kpts,n,parm) result(val)
    ! Get the values of the quadratic function on a list of kpoints
    integer,intent(in) :: n
    real(dp),intent(in) :: kpts(3,n)
    real(dp),intent(in) :: parm(10)
    real(dp) :: val(n)
    integer :: ii
    ! evaluate a functions on a list of points
    do ii=1,n
      val(ii) = f(kpts(:,ii),parm)
    end do
  end function

  subroutine get_8tetra(eig10,mat10,wvals,nw,integral)
    ! sum the contribution of 8 tetrahedra generated by dividing one tetrahedron
    integer,intent(in) :: nw
    real(dp),intent(in) :: eig10(10), mat10(10)
    real(dp),intent(in) :: wvals(nw)
    real(dp),intent(out) :: integral(nw)

    integer :: itetra
    integer :: ind(4)
    real(dp) :: eig4(4),mat4(4)

    ! Evaluate and add the contribution of 8 tetrahedra
    integral = 0.0_dp
    do itetra=1,8
      eig4(:) = eig10(tetra8(:,itetra))
      mat4(:) = mat10(tetra8(:,itetra))
      ind = [1,2,3,4]
      call sort_4tetra(eig4,ind)
      call tetralite_delta(eig4,wvals,nw,mat4(ind),integral)
    end do
  end subroutine

  recursive subroutine hybridtetra(kpt4,eig4,mat4,parm_eig,parm_mat,wvals,nw,integral,n,depth)
    ! subdivide tetrahedra n times and evaluate contibutions
    ! number of frequencies, desired depth
    integer,intent(in) :: nw,n
    ! quadratic coeficients for interpolation
    real(dp),intent(in) :: parm_eig(10),parm_mat(10)
    real(dp),intent(in) :: kpt4(3,4),eig4(4),mat4(4)
    real(dp),intent(in) :: wvals(nw)
    real(dp),intent(inout) :: integral(nw)
    integer,optional,intent(in) :: depth
    integer  :: itetra,idepth
    integer  :: ind(4)
    real(dp) :: eig(4),mat(4)
    real(dp) :: mid10(3,10),eig10(10),mat10(10)

    ! how deep are we?
    idepth = 1; if (present(depth)) idepth = depth
    if (idepth==1) integral = 0.0_dp

    ! if in the desired depth sum the contribution of this tetrahedron
    if (idepth==n) then
      eig = eig4
      mat = 1.0_dp/8.0_dp**(idepth-1)
      call sort_4tetra(eig,ind)
      call tetralite_delta(eig,wvals,nw,mat,integral)
    ! otherwise go deeper for each tetrahedron
    else
      ! copy 4 vertices data
      mid10(:,:4) = kpt4
      eig10(:4)   = eig4
      mat10(:4)   = mat4
      ! get the additional 6 middle points
      mid10(:,5:) = get_6midkpts(kpt4)
      eig10(5:)   = get_values(mid10(:,5:),6,parm_eig)
      mat10(5:)   = get_values(mid10(:,5:),6,parm_mat)
      ! Loop over the 8 sub tetrahedra
      do itetra=1,8
        ! get coordinates of k-points of this tetrahedron
        ind(:) = tetra8(:,itetra)
        ! subdivide and evaluate again
        call hybridtetra(mid10(:,ind),eig10(ind),mat10(ind),parm_eig,parm_mat,wvals,nw,integral,n,idepth+1)
      end do
    end if
  end subroutine
end module

program fit
  use m_quadratic
  use m_dividetetra
  use m_tetralite
  implicit none

  integer :: ii,info
  real(dp) :: kpts(3,10), eigs(10), vels(3,10)
  real(dp) :: parm_mat(10), parm_eig(10)

  ! Define a function to parametrize eigenvalues
  parm_eig = [ 1.0_dp,&
               1.0_dp, 1.0_dp, 2.0_dp,&
               2.0_dp, 1.2_dp, 2.0_dp,&
               1.0_dp, 8.0_dp, 1.0_dp]

  parm_mat = [ 1.0_dp,&
               0.0_dp, 0.0_dp, 0.0_dp,&
               0.0_dp, 0.0_dp, 0.0_dp,&
               0.0_dp, 0.0_dp, 0.0_dp]

  ! Generate 10 random points
  call random_number(kpts)

  ! Test fitting a quadratic function using 10 points or 4 points + derivatives
  fit_stuff : block
    ! Set eigenvalues and velocities
    do ii=1,10
      eigs(ii)   = f(kpts(:,ii),parm_eig)
      vels(:,ii) = df(kpts(:,ii),parm_eig)
    end do
    write(*,'(10f8.3)') (parm_eig(ii),ii=1,10)

    ! Fit using 10 points
    call fit10(kpts,eigs,parm_eig,info)
    write(*,'(10f8.3)') (parm_eig(ii),ii=1,10)

    ! Fit using 4 points and derivatives
    call fit4d(kpts,eigs,vels,parm_eig,info)
    write(*,'(10f8.3)') (parm_eig(ii),ii=1,10)
  end block fit_stuff

  ! now its a bit of tetrahedron stuff
  ! subdivide a tetrahedron in 8 and get the value of the integral inside it
  tetradivide_stuff : block
    integer :: nw
    real(dp) :: mid6(3,6)
    real(dp) :: eig4(4), mat4(4)
    real(dp) :: eig10(10), mat10(10)
    real(dp),allocatable :: wvals(:)
    real(dp),allocatable :: integral(:)

    ! Initialize
    parm_eig = parm_eig
    call tetralite_init()
    call dividetetra_init()

    nw = 10
    allocate(wvals(nw))
    allocate(integral(nw))
    wvals = [(ii,ii=1,nw)]

    ! 4 of the points come from the mother
    eig10(:4) = eigs(:4)
    ! Get the additional 6 middle points and their eigenvalues
    mid6 = get_6midkpts(kpts)
    eig10(5:) = get_values(mid6,6,parm_eig)
    ! Here I set the matrix elements to constants.
    ! Could interpolate them as well (linearly or quadratically)
    mat10 = 1.0_dp
    call get_8tetra(eig10,mat10,wvals,nw,integral)

    ! now test the hybrid tetrahedron routine
    eig4 = get_values(kpts,4,parm_eig)
    mat4 = get_values(kpts,4,parm_mat)
    call  hybridtetra(kpts(:,:4),eig4,mat4,parm_eig,parm_mat,wvals,nw,integral,1)
    deallocate(wvals)
    deallocate(integral)

  end block tetradivide_stuff

  tetralite_stuff: block
    integer :: ix,iy,iz
    integer :: ik,nk,nw
    integer :: itetra,isummit,idepth
    integer :: jtetra,jsummit,idx10
    integer :: band
    integer :: divs(3),ind(4),gp(3)
    character(len=100) :: fname
    real(dp) :: emin,emax,step
    real(dp) :: start_time, stop_time
    real(dp) :: kpt(3),kpt4(3,4),eig4(4),mat4(4),vel4(3,4)
    real(dp) :: kpt10(3,10),eig10(10),mat10(10)
    real(dp) :: parm_eig(10)
    real(dp),allocatable :: eig(:),vel(:,:)
    real(dp),allocatable :: wvals(:)
    real(dp),allocatable :: dweight(:,:)
    real(dp),allocatable :: integral(:), integral_tmp(:)

    ! choose TB band or parabola
    band = 1

    ! choose how many times to apply recursion
    idepth = 1
    divs = [2,2,2]
    divs = [4,4,4]
    divs = [6,6,6]
    divs = [8,8,8]
    divs = [10,10,10]
    !divs = [20,20,20]
    !divs = [30,30,30]
    !divs = [40,40,40]
    !divs = [50,50,50]
    divs = [100,100,100]
    nw = 500

    write(*,*) divs
    write(*,*) idepth

    nk = divs(1)*divs(2)*divs(3)
    allocate(eig(nk))
    allocate(vel(3,nk))
    allocate(wvals(nw))
    allocate(dweight(4,nw))
    allocate(integral(nw))
    allocate(integral_tmp(nw))

    ! Evaluate function on each kpoint
    ! We want to precompute the values of the function in the BZ and then use then to compute DOS
    ! Note that we assume that the functions are periodic, if this is not the case
    ! then the precomputed velocities will not be correct.
    ! this is the case for the parabola_band
    ! For testing purposes we can recompute the velocities for the quadratic tetrahedra bellow
    do ix=1,divs(1)
      do iy=1,divs(2)
        do iz=1,divs(3)
          kpt = grid_kpoint(divs,[ix,iy,iz])
          ik = grid(divs,[ix,iy,iz])
          select case(band)
          case(1)
           call tb_band(kpt,eig(ik),vel(:,ik))
          case(2)
           call parabola_band(kpt,eig(ik),vel(:,ik))
          end select
        end do
      end do
    end do
    ! For the moment I set the matrix elements to zero
    mat4  = one
    mat10 = one


    ! Determine DOS energy range
    emin = minval(eig)-0.2_dp
    emax = maxval(eig)+0.2_dp
    step = (emax-emin)/(nw-1)
    wvals = [(emin+ii*step,ii=0,nw-1)]

    ! Compute DOS with linear tetrahedron
    call cpu_time(start_time)
    integral = zero
    mat4 = one/nk
    do ix=1,divs(1)
      do iy=1,divs(2)
        do iz=1,divs(3)
          ! generate 6 tetrahedra for this point
          do itetra=1,6
            do isummit=1,4
              gp = [ix,iy,iz]+tetra_6lshifts(:,isummit,itetra)
              ! Get eigenvalues and velocities
              ik = grid(divs,gp)
              eig4(isummit)   = eig(ik)
              vel4(:,isummit) = vel(:,ik)
              ! Compute eigenvalues and velocities
              !kpt = grid_kpoint(divs,gp)
              !call tb_band(kpt,eig4(isummit),vel4(:,isummit))
              !call parabola_band(kpt,eig4(isummit),vel4(:,isummit))
            end do
            call sort_4tetra(eig4,ind)
            call tetralite_delta(eig4,wvals,nw,mat4,integral)
          end do
        end do
      end do
    end do
    call cpu_time(stop_time)
    write(*,'(a30,f12.6,a)') 'linear dos took:',(stop_time-start_time), ' [s]'
    write(fname,'(a,i0,a)') 'dosl',divs(1),'.dat'
    call write_file(fname,nw,wvals,integral)

    ! Compute DOS with quadratic tetrahedron (fitting only scalars)
    call cpu_time(start_time)
    integral = zero
    mat4=one
    ! loop over points
    do ix=1,divs(1),2
      do iy=1,divs(2),2
        do iz=1,divs(3),2
          ! generate 6 tetrahedra for this point
          do itetra=1,6
            ! each tetrahdra has 10 points
            do isummit=1,10
              gp = [ix,iy,iz]+tetra_6qshifts(:,isummit,itetra)
              ! Get eigenvalues
              ik = grid(divs,gp)
              eig10(isummit)   = eig(ik)
              kpt10(:,isummit) = grid_kpoint(divs,gp)
            end do
            ! fit energies
            call fit10(kpt10,eig10,parm_eig,info)
            ! now apply hibrid tetrahedron to each of the 8 tetrahedra
            do jtetra=1,8
              do jsummit=1,4
                idx10 = tetra8(jsummit,jtetra)
                eig4(jsummit)   = eig10(idx10)
                kpt4(:,jsummit) = kpt10(:,idx10)
              end do
              ! get hybrid weights
              call hybridtetra(kpt4,eig4,mat4,parm_eig,parm_mat,wvals,nw,integral_tmp,idepth)
              integral = integral + integral_tmp/nk
            end do
          end do
        end do
      end do
    end do
    call cpu_time(stop_time)
    write(*,'(a30,f12.6,a)') 'quadratic scalar dos took:',(stop_time-start_time), ' [s]'
    write(fname,'(a,i0,a,i0,a)') 'dosqs',divs(1),'d',idepth,'.dat'
    call write_file(fname,nw,wvals,integral)

    ! Compute DOS with quadratic tetrahedron (fitting only scalars and using affine transformation)
    call cpu_time(start_time)
    integral = zero
    ! loop over points
    kpt10(:,1)  = [0.0_dp,0.0_dp,0.0_dp]
    kpt10(:,2)  = [1.0_dp,0.0_dp,0.0_dp]
    kpt10(:,3)  = [0.0_dp,1.0_dp,0.0_dp]
    kpt10(:,4)  = [0.0_dp,0.0_dp,1.0_dp]
    kpt10(:,5:) = get_6midkpts(kpt10(:,:4))
    do ix=1,divs(1),2
      do iy=1,divs(2),2
        do iz=1,divs(3),2
          ! generate 6 tetrahedra for this point
          do itetra=1,6
            ! each tetrahdra has 10 points
            do isummit=1,10
              gp = [ix,iy,iz]+tetra_6qshifts(:,isummit,itetra)
              ! Get eigenvalues
              ik = grid(divs,gp)
              eig10(isummit) = eig(ik)
            end do
            ! fit energies
            call fit10_affine(eig10,parm_eig)
            ! now apply hibrid tetrahedron to each of the 8 tetrahedra
            do jtetra=1,8
              do jsummit=1,4
                idx10 = tetra8(jsummit,jtetra)
                eig4(jsummit) = eig10(idx10)
                kpt4(:,jsummit) = kpt10(:,idx10)
              end do
              ! get hybrid weights
              call hybridtetra(kpt4,eig4,mat4,parm_eig,parm_mat,wvals,nw,integral_tmp,idepth)
              integral = integral + integral_tmp/nk
            end do
          end do
        end do
      end do
    end do
    call cpu_time(stop_time)
    write(*,'(a30,f12.6,a)') 'quadratic affine dos took:',(stop_time-start_time), ' [s]'
    write(fname,'(a,i0,a,i0,a)') 'dosqsa',divs(1),'d',idepth,'.dat'
    call write_file(fname,nw,wvals,integral)



    ! Compute DOS with quadratic tetrahedron (fitting the velocities)
    call cpu_time(start_time)
    integral = zero
    ! loop over points
    do ix=1,divs(1)
      do iy=1,divs(2)
        do iz=1,divs(3)
          ! generate 6 tetrahedra for this point
          do itetra=1,6
            do isummit=1,4
              gp = [ix,iy,iz]+tetra_6lshifts(:,isummit,itetra)
              ! Get eigenvalues and velocities
              ik = grid(divs,gp)
              eig4(isummit)   = eig(ik)
              vel4(:,isummit) = vel(:,ik)
              kpt4(:,isummit) = grid_kpoint(divs,gp)
              ! Activate this for the parabola band
              ! Compute eigenvalues and velocities
              !kpt4(:,isummit) = grid_kpoint(divs,gp)
              !call tb_band(kpt4(:,isummit),eig4(isummit),vel4(:,isummit))
              !call parabola_band(kpt4(:,isummit),eig4(isummit),vel4(:,isummit))
            end do
            ! fit energies
            call fit4d(kpt4,eig4,vel4,parm_eig,info)
            ! get hybrid weights
            call hybridtetra(kpt4,eig4,mat4,parm_eig,parm_mat,wvals,nw,integral_tmp,idepth)
            integral = integral + integral_tmp/nk
          end do
        end do
      end do
    end do
    call cpu_time(stop_time)
    write(*,'(a30,f12.6,a)') 'quadratic velocity dos took:',(stop_time-start_time), ' [s]'
    write(fname,'(a,i0,a,i0,a)') 'dosq',divs(1),'d',idepth,'.dat'
    call write_file(fname,nw,wvals,integral)

    deallocate(eig)
    deallocate(vel)
    deallocate(wvals)
    deallocate(dweight)
    deallocate(integral)
    deallocate(integral_tmp)

  end block tetralite_stuff

  contains
  pure subroutine parabola_band(kpt,eig,vel)
    ! Simple parabola
    real(dp),intent(in) :: kpt(3)
    real(dp),intent(out) :: eig, vel(3)
    real(dp) :: x,y,z
    x = kpt(1)*two_pi
    y = kpt(2)*two_pi
    z = kpt(3)*two_pi
    eig = x*x+y*y+z*z
    vel(1) = 2*two_pi*x
    vel(2) = 2*two_pi*y
    vel(3) = 2*two_pi*z
  end subroutine

  pure subroutine tb_band(kpt,eig,vel)
    ! TB band
    real(dp),intent(in) :: kpt(3)
    real(dp),intent(out) :: eig, vel(3)
    real(dp) :: x,y,z
    x = kpt(1)*two_pi
    y = kpt(2)*two_pi
    z = kpt(3)*two_pi
    eig = -(cos(x)*cos(y) + cos(x)*cos(z) + cos(y)*cos(z))
    vel(1) = two_pi*(cos(y)*sin(x) + cos(z)*sin(x))
    vel(2) = two_pi*(cos(x)*sin(y) + cos(z)*sin(y))
    vel(3) = two_pi*(cos(x)*sin(z) + cos(y)*sin(z))
  end subroutine

  subroutine write_file(fname,nw,wvals,integral)
    integer :: iw,nw
    real(dp) :: wvals(nw), integral(nw)
    character(len=100) :: fname
    open(unit=1,file=fname)
    do iw=1,nw
      write(1,*) wvals(iw), integral(iw)
    end do
    close(1)
  end subroutine
end program fit

