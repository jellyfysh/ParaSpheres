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
