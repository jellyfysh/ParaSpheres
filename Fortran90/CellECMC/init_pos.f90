subroutine init_pos
!******************************************************************!
!
! Creates an initial condition for a squared number of disks
!
! Works with square or rectangular box
!
!******************************************************************!
  use phys_parameters
  use others
  implicit none
  integer*8 :: i,j,k,n
  real*8, dimension(1:2) ::dx=(/1.D0,0.D0/),dy=(/0.5D0,1.D0/)
  logical :: ex
  inquire(file='init.dat', exist=ex)
  if (ex) then
     open(1,file='init.dat')
     write(*, *) 'Existing initial state'
     do i=1,N_b
        read(1,*) X(:,i)
     enddo
  else
     n=ceiling(N_b**.5)
     dx=dx*Box/n
     dy=dy*Box/n
     do i=0,n-1
        do j=0,n-1
           if (j*n+i+1>N_b) cycle
           X(:,j*n+i+1)=mod(i*dx+j*dy,Box)
           do k=1,2
              if (X(k,j*n+i+1)<=-Box(k)/2) X(k,j*n+i+1)=X(k,j*n+i+1)+Box(k)
              if (X(k,j*n+i+1)>=Box(k)/2)  X(k,j*n+i+1)=X(k,j*n+i+1)-Box(k)
           enddo
        enddo
     enddo
  endif
end subroutine init_pos
