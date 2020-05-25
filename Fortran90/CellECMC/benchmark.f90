module benchmark
  implicit none 
  logical :: do_benchmark
  logical :: create_init
  real*8 :: precision=1d-14
end module benchmark

subroutine init_benchmark
!******************************************************************!
!
!argument parser
!
!******************************************************************!
  use benchmark
  implicit none
  integer*8 :: number_arguments
  character(100) :: argument1
  character(100) :: argument2

  number_arguments = command_argument_count()
  create_init=.false.
  if (number_arguments==0) then
  	do_benchmark = .false.
  elseif (number_arguments>2) then
  	write(*,*) 'I dont know what you are talking about.'
  	do_benchmark = .false.
  elseif (number_arguments==1) then
  	call get_command_argument(1, argument1)
  	if (argument1=='benchmark') then
  	  write(*,*) 'Benchmark mode'
  	  write(*,'("Precision for three-body-collision detection  ",f18.16, "  vs  ", f18.16)') 1d0, 1d0-precision
  	  do_benchmark = .true.
  	elseif (argument1=='create_init') then
  	  write(*,*) 'will create new initial configuration'
  	  create_init=.true.
  	else
  	  write(*,*) 'I dont know what you are talking about.'
  	  do_benchmark = .false.
  	endif
  elseif (number_arguments==2) then
  	call get_command_argument(1, argument1)
  	call get_command_argument(2, argument2)
  	if (argument1=='benchmark') then
  	  write(*,*) 'Benchmark mode'
  	  write(*,'("Precision for three-body-collision detection",f18.16, f18.16)') 1d0, 1d0-precision
  	  do_benchmark = .true.
  	endif
  	if (argument2=='create_init') then
  	  write(*,*) 'will create new initial configuration'
  	  create_init=.true.
  	endif
  endif
end subroutine init_benchmark

subroutine precision_check(a,b)
!******************************************************************!
!
!Detect three body collision. Call in expl_x_cell, expl_y_cell and the main program
!
!******************************************************************!
  use benchmark
  implicit none
  real*8 :: a,b
  if (abs(a-b)<precision) then
  	write(*,*) 'The simulation is at a crossroads.'
        !call EXIT(1)
  endif
 end subroutine precision_check

subroutine write_init
!******************************************************************!
!
! Output the configuration center in (0,0) in init.dat
!
!******************************************************************!
  use phys_parameters
  use cell_parameters
  use others
  implicit none
  integer*8 :: j,k
  open(1,file='init.dat')
  do j=1,N_cell(1)*N_cell(2)
     do k=1,Cell_ocp(j)
        write(1,'(f18.10,f18.10)') Cell(1,k,j)-box(1)/2+(modulo(j-1,N_cell(1))+0.5)*l_cell(1),&
             & Cell(2,k,j)-box(2)/2+((j-1)/N_cell(1)+0.5)*l_cell(2)
     enddo
  enddo
  close(1)

  open(1,file='info.dat')
  write(1, '(A1, I12)') "N", N_b
  write(1, '(A5, f22.18)') "sigma", r   
  close(1)
end subroutine write_init

subroutine write_hist_cells
!******************************************************************!
!
! Output the number of particles in each cell in hist_cells.dat
!
!******************************************************************!
  use phys_parameters
  use cell_parameters
  use others
  implicit none
  integer*8 :: j,k
  open(1,file='hist_cells.dat')
  do j=1,N_cell(1)*N_cell(2)
     write(1,'(i2)') Cell_ocp(j)
  enddo
  close(1)
end subroutine write_hist_cells
