subroutine write_pos
!******************************************************************!
!
! Output the configuration center in (0,0) in pos.dat
!
!******************************************************************!
  use phys_parameters
  use cell_parameters
  use others
  implicit none
  integer*8 :: j,k
  open(1,file='pos.dat',position='append')
  write(1,'(i17)')
  write(1,'(i17,e10.3)') T,time
  do j=1,N_cell(1)*N_cell(2)
     do k=1,Cell_ocp(j)
        write(1,'(f9.6,f10.6)') Cell(1,k,j)-box(1)/2+(modulo(j-1,N_cell(1))+0.5)*l_cell(1),&
             & Cell(2,k,j)-box(2)/2+((j-1)/N_cell(1)+0.5)*l_cell(2)
     enddo
  enddo
  close(1)
end subroutine write_pos
