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

subroutine pick_move(cell_cur,k_cur)
!******************************************************************!
!
!Picks one disk in the cells. equi-probability for each disk
!
!******************************************************************!
  use alg_parameters
  use phys_parameters
  use cell_parameters
  use benchmark
  implicit none
  integer*8 :: k_cur,cell_cur
  real*8 ::g
  if (do_benchmark) then
    if (cell_cur/=-1) then
      goto 123
    else
      k_cur=1
      cell_cur=1
      do
        if (cell_ocp(cell_cur)<1) then
          cell_cur=cell_cur+1
        else
          exit
        endif
      enddo
      write(*,'("Initial cell=",i2,"  initial active particle=(",d10.4,", ",d10.4,")")') cell_cur, &
          Cell(1,k_cur,cell_cur), Cell(2,k_cur,cell_cur)
      write(*,'("Global position = (",d12.6,", ",d12.6,")")') Cell(1,k_cur,cell_cur)-box(1)/2+&
          (modulo(cell_cur-1,N_cell(1))+0.5)*l_cell(1),&
          Cell(2,k_cur,cell_cur)-box(2)/2+((cell_cur-1)/N_cell(1)+0.5)*l_cell(2)
      goto 123
    endif
  endif
  do
     if (do_benchmark) write(*,*) 'choose a random active particle'
     call random_number(g)
     cell_cur=floor(g*N_cell(1)*N_cell(2))+1
     call random_number(g)
     k_cur=floor(g*n_cell_max)+1
     if (k_cur>cell_ocp(cell_cur)) cycle
     exit
  enddo
123 continue
end subroutine pick_move
