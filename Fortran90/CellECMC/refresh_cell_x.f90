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

subroutine refresh_cell_x(dx,cell_cur,k_cur)
!******************************************************************!
!
! Compute the new cell of k_cur
!
!******************************************************************!
  use phys_parameters
  use cell_parameters
  implicit none
  real*8 :: dx
  integer*8:: cell_new,cell_cur,k_cur
  if (dx>l_cell(1)/2) then                                      !disk k_cur is now in cell_new=Cell_neighb(6,cell_cur)
     dx=dx-l_cell(1)                                            ! position in the new cell
     cell_new=Cell_neighb(6,cell_cur)                           !new cell
     Cell_ocp(cell_new)=Cell_ocp(cell_new)+1                    ! refresh # of disk in the new cell
     Cell(1,Cell_ocp(cell_new),cell_new)=dx                     ! x position in the new cell
     Cell(2,Cell_ocp(cell_new),cell_new)=Cell(2,k_cur,cell_cur) ! y position in the new cell
     Cell(:,k_cur,cell_cur)=Cell(:,Cell_ocp(cell_cur),cell_cur) !swap the last disk in the old cell to take k_cur's vacant place
     Cell_ocp(cell_cur)=Cell_ocp(cell_cur)-1                    ! refresh # of disks in the old cell
     cell_cur=cell_new                                          ! The cell of the "current" disk is now cell_new....
     k_cur=Cell_ocp(cell_new)                                   ! and its place is the last one
  else                                                          !disk k_cur stays in cell_cur
     Cell(1,k_cur,cell_cur)=dx
  endif
end subroutine refresh_cell_x
