! JeLLyFysh/ParaSpheres - Multithreaded event-chain Monte Carlo with local times
! - https://github.com/jellyfysh/paraspheres
! Copyright (C) 2020 The JeLLyFysh organization
! (see the AUTHORS file on jellyfysh/paraspheres for the full list of authors)
!
! This file is part of JeLLyFysh/ParaSpheres.
!
! JeLLyFysh/ParaSpheres is free software: you can redistribute it and/or modify
! it under the terms of the GNU General
! Public License as published by the Free Software Foundation, either > version
! 3 of the License, or (at your option)
! any later version.
!
! JeLLyFysh/ParaSpheres is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even 
! the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
! See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along with
! JeLLyFysh/ParaSpheres in the LICENSE file. If not, see
! <https://www.gnu.org/licenses/>.
!
! If you use JeLLyFysh/ParaSpheres in published work, please cite the following
! reference
! (see [Li2020] in References.bib):
! Botao Li, Synge Todo, A. C. Maggs, Werner Krauth
! Multithreaded event-chain Monte Carlo with local times,
! arXiv e-prints: 2004.11040 (2020), https://arxiv.org/abs/2004.11040
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
