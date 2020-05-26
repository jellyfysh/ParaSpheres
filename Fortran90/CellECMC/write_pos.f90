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
