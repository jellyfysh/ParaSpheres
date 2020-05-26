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
module alg_parameters
  use phys_parameters
  implicit none
  real*8, parameter :: lambda_0=0.0768D0/sqrt(real(N_b))     !expected mean freepath, used to compute L and N_extr
  real*8, parameter :: n_dpl=3.125D0*sqrt(real(N_b))      !L=n_dpl*lambda_0 (n_dpl ~ <# displ.>, 10% error for big L)
  real*8 :: L=lambda_0*n_dpl                          !total displacement ell
  integer*8 :: N_dpl_tot=1D04*N_b                  !total # displ. for the run (-1 for infinity)
  integer*8,parameter :: N_dpl_extr=1D03*N_b          !extraction rate of configurations in number of displacement, writes in 'pos.dat'
  integer*8,parameter :: N_extr=(N_dpl_extr-1)/(n_dpl+1)+1        !convert N_dpl_extr in ~ # of move
end module alg_parameters
