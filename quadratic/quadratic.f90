!
! 20/03/2019
! Create a linear system of equations to fit a quadratic function in 3D
! the function has the general form:
!
!  f = a +
!      ax*x + ay*y + az*z +
!      axy*x*y + axz*x*z + ayz*y*z +
!      ax2*x*x + ayy*y*y + azz*z*z
!
module m_quadratic

 implicit none
 integer, parameter :: dp=8
 real(dp), parameter :: half = 0.5_dp

 contains
 real(dp) function f(kpt,parm)
   real(dp),intent(in) :: kpt(3)
   real(dp),intent(in) :: parm(10)
   real(dp) :: x,y,z
   x = kpt(1)
   y = kpt(2)
   z = kpt(3)
   f = parm(1) + &
       parm(2)*x   + parm(3)*y   + parm(4)*z + &
       parm(5)*x*y + parm(6)*x*z + parm(7)*y*z + &
       parm(8)*x*x + parm(9)*y*y + parm(10)*z*z
 end function

 function df(kpt,parm)
   ! Compute the derivative of the quadratic function with respect to x,y,z
   real(dp),intent(in) :: kpt(3)
   real(dp),intent(in) :: parm(10)
   real(dp) :: df(3)
   real(dp) :: x,y,z
   x = kpt(1)
   y = kpt(2)
   z = kpt(3)
   df(1) = parm(2) + &
           parm(5)*y + parm(6)*z + &
           2.0_dp*parm(8)*x
   df(2) = parm(3) + &
           parm(5)*x + parm(7)*z + &
           2.0_dp*parm(9)*y
   df(3) = parm(4) + &
           parm(6)*x + parm(7)*y + &
           2.0_dp*parm(10)*z
 end function

 subroutine fit10(kpts,eigs,parm,info)
  ! Fit a quadratic function using the values fo the function evaluated in 10 points
  ! This is done solving a linear system of equations of the form Ax=b
  real(dp),intent(in) :: kpts(3,10)
  real(dp),intent(in) :: eigs(10)
  real(dp),intent(out) :: parm(10)
  integer,intent(out) :: info

  integer :: ii,irank
  real(dp) :: x,y,z
  real(dp) :: work(100),sgval(10)
  real(dp) :: a(10,10), b(10)

  ! Fill in the a matrix A
  do ii=1,10
    x=kpts(1,ii)
    y=kpts(2,ii)
    z=kpts(3,ii)
    a(ii,:) = [1.0_dp,x,y,z,x*y,x*z,y*z,x*x,y*y,z*z]
  end do

  ! Fill in the vector b
  b = eigs

  !solve a system of linear equations
  call dgelss(10,10,1,a,10,b,10,sgval,1.0d-8,irank,work,100,info)
  parm = b(:10)
 end subroutine

 subroutine fit4(kpts,eigs,vels,parm,info)
  ! Fit a quadratic function using the values and derivatives of the function evaluated in 4 points
  ! This is done solving a linear system of equations of the form Ax=b
  real(dp),intent(in) :: kpts(3,4)
  real(dp),intent(in) :: eigs(4)
  real(dp),intent(in) :: vels(3,4)
  real(dp),intent(out) :: parm(10)
  integer,intent(out) :: info

  integer :: ii,ieq
  integer :: irank
  real(dp) :: x,y,z
  real(dp) :: work(100),sgval(10)
  real(dp) :: a(16,10),b(16)

  ! Fill in the a matrix for eigenvalues
  do ii=1,4
    x=kpts(1,ii)
    y=kpts(2,ii)
    z=kpts(3,ii)
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
  b(:4)  = eigs
  b(5:7) = vels(:,1)
  b(8:10) = vels(:,2)
  b(11:13) = vels(:,3)
  b(14:16) = vels(:,4)

  !solve a system of linear equations
  call dgelss(16,10,1,a,16,b,16,sgval,1.0d-8,irank,work,100,info)
  parm = b(:10)
 end subroutine

end module

program fit

  use m_quadratic
  implicit none

  integer :: ii,info
  real(dp) :: mid(3), kpts(3,10), eigs(10), vels(3,10), parm(10)

  ! Define a function to parametrize
  parm = [1.0_dp,&
          1.0_dp,1.0_dp,1.0_dp,&
          2.0_dp,1.0_dp,1.0_dp,&
          1.0_dp,8.0_dp,1.0_dp]

  ! Generate 10 random points
  call random_number(kpts)

  ! Set eigenvalues and velocities
  do ii=1,10
    eigs(ii)   = f(kpts(:,ii),parm)
    vels(:,ii) = df(kpts(:,ii),parm)
  end do

  ! Fit using 10 points
  call fit10(kpts,eigs,parm,info)
  write(*,'(10f8.3)') (parm(ii),ii=1,10)

  ! Fit using 4 points and derivatives
  call fit4(kpts,eigs,vels,parm,info)
  write(*,'(10f8.3)') (parm(ii),ii=1,10)

  ! The first 4 energies come from above
  ! Get remaining 6 energies from midlepoints
  mid = half*(kpts(:,1)+kpts(:,2)); eigs(5)  = f(mid,parm)
  mid = half*(kpts(:,1)+kpts(:,3)); eigs(6)  = f(mid,parm)
  mid = half*(kpts(:,1)+kpts(:,4)); eigs(7)  = f(mid,parm)
  mid = half*(kpts(:,2)+kpts(:,3)); eigs(8)  = f(mid,parm)
  mid = half*(kpts(:,2)+kpts(:,4)); eigs(9)  = f(mid,parm)
  mid = half*(kpts(:,3)+kpts(:,4)); eigs(10) = f(mid,parm)

end program fit

