module phys_parameters
  implicit none
  real*8, parameter :: pi=acos(-1.D0)                !pi
  integer*8, parameter :: N_b=( 127 )**2             !number of disk
  real*8, parameter :: eta= 0.708              !density (volume fraction)
  real*8, parameter,dimension(1:2) :: Box=(/1.D0,1.D0/)    !Box's size
  real*8, parameter :: r=sqrt(eta*Box(1)*Box(2)/(N_b*pi))    !radius
  real*8, dimension(1:2,1:N_b) :: X
end module phys_parameters
