module cell_parameters
  use phys_parameters
  integer*8, parameter,dimension(1:2) :: N_cell=(/1,1/)*nint(sqrt(real(N_b))*7.D0/8.D0)    !Cells number (1024 ~ 28*28 etc)
  real*8, parameter,dimension(1:2) :: l_cell=Box/N_cell        !cells size
  integer*8, parameter :: n_cell_max=floor((l_cell(1)/r+2.D0)*(l_cell(2)/r+2.D0)/(2.D0*sqrt(3.D0)))!>max # of disks/cell
  real*8, dimension(1:2,1:n_cell_max,1:N_cell(1)*N_cell(2)) :: Cell  !positions of the disks/ center of their cell
  integer*2, dimension(1:N_cell(1)*N_cell(2)) :: Cell_ocp=0          !# of disks for each cell
  integer*4,dimension(9,N_cell(1)*N_cell(2)) :: Cell_neighb          !Cell neighbor's table
  integer*8,dimension(1:2) :: i_cell
end module cell_parameters
