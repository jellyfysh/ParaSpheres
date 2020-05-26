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
subroutine init_cell
!*********************************************************************!
!
! Establish the cell-scheme and put the inital configuration in the cells
! Restart from last valid conf of 'pos.dat'
!
!*********************************************************************!
  use cell_parameters
  use phys_parameters
  use alg_parameters
  use others
  use benchmark
  implicit none
  integer*8 :: i,i_2,j_2,j,c_1,c_2,n_1,n_2
  logical :: ex
  character ( len = 50 ) :: line
  inquire(file='pos.dat', exist=ex)
  if (do_benchmark) ex=.false.
  if (ex) then
     write(*,*) 'Restarts with the last configuration'
     open(1,file='pos.dat',position='append')
     !***check if the last configuration is complete ***!
     i=0
     do
        backspace(1)
        read(1, '(a)') line
        if (len_trim(line)==0) exit
        i=i+1
        backspace(1)
     enddo
     if (i==N_b+1) then !*** config complete
        read(1,*) T,t_0
        do j=1,N_b
           read(1,*) X(:,j)
        enddo
        close(1)
     else !*** config not complete
        write(*,*) 'configuration is not complete : restart with the previous one'
        backspace(1)
        backspace(1)
        read(1, '(a)') line
        backspace(1)
        write(1, '(2a)') line,'restart'
        close(1)
        open(1,file='pos.dat',position='append')
        backspace(1)
        do j=0,N_b-1
           read(1,*) X(:,N_b-j)
           backspace(1)
           backspace(1)
        enddo
        read(1,*) T,t_0
        close(1)
     endif
  else
     write(*,*) 'Creates an initial configuration'
     open(1,file='pos.dat')
     write(1,*) '# Positions'
     write(1,*) '# SEC_XY'
     write(1,*) '# Box=',Box
     write(1,*) '# N_b=',N_b
     write(1,*) '# r=',r
     write(1,*) '# eta=', eta
     write(1,*) '# Lambda_0=',lambda_0
     write(1,*) '# L=',n_dpl
     write(1,*) '# N_dpl_extr=',N_dpl_extr
     write(1,*) '#'
     write(1,'(f8.5,f9.5,i8,e13.6)') Box, N_b, r
     close(1)
     call init_pos !intial condition
  endif
 !*** create the neighbor table ***!
  do j=1,N_cell(2)
     do i=1,N_cell(1)
        do c_2=-1,1
           do c_1=-1,1
              n_1=i+(j-1)*N_cell(1)
              i_2=modulo(i+c_1-1,N_cell(1))+1
              j_2=modulo(j+c_2-1,N_cell(2))+1
              n_2=i_2+(j_2-1)*N_cell(1)
              Cell_neighb(c_1+2+(c_2+1)*3,n_1)=n_2
           enddo
        enddo
     enddo
  enddo
  !*** end neighbor table ***!

  !*** put the initial condition in the Cells ***!
  do i=1,N_b
     do j=1,2
        i_cell(j)=floor((X(j,i)+Box(j)/2)/l_cell(j))+1
        if (i_cell(j)==N_cell(j)+1) then
           i_cell(j)=1
           X(j,i)=X(j,i)-Box(j)
        endif
     enddo
     n_1=i_cell(1)+(i_cell(2)-1)*N_cell(1)
     Cell_ocp(n_1)=Cell_ocp(n_1)+1
     do j=1,2
        Cell(j,Cell_ocp(n_1),n_1)=X(j,i)+Box(j)/2-(real(i_cell(j))-0.5D0)*l_cell(j)
     enddo
  enddo
  !*** end fill cells ***!
end subroutine init_cell
