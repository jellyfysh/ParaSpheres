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
!
! This code was written by E. P. Bernard. If you use it in published work, please
! cite the following reference
! (see [Bernard2011] in References.bib):
! Etienne P. Bernard, Werner Krauth
! Two-Step Melting in Two Dimensions: First-Order Liquid-Hexatic Transition,
! Physical Review Letters 107, 155704, (2011)
!*************************************************************************!
!
!     Event_chain_SEC-XY  Version 28/04/2010
!
!   - Monodisperse hard disks in the (N,V,T) ensemble (faster than (N,P,T))
!   - Disks go straight in (ox) then (oy) then (ox)... until the total displacement
!     L is reached.
!   - ~20 times faster than the local Metropolis 
!   - Periodic boundary conditions
!   - Extract the configurations, centered in (0,0), each N_extr move (in 'pos.dat')
!   - Irreversible (faster) 
!   - Disks do not have any number (cannot compute diffusion constant)
!   - Disks positions are relative to the center of their cell (stored in cell)
!   - Very stable, overlapping is not possible even with low precision
!   - Restart from last valid conf of 'pos.dat'.
!   - On a unique Core2 Q9550 @ 2.83Ghz N=64², ~ 28E9 displacements/hour (gfortran -03)
!                                              ~ 32E9 displacements/hour (ifort -O3)
!     (less for very big system.. memory..) N=1024² ~ 1E10 dpl/hour
!   - Extract position and backup each N_extr
!   - Extract Correl_o each N_extr_correl
!
!*************************************************************************!

module phys_parameters
  implicit none
  real*4, parameter :: pi=acos(-1.D0)
  integer*8, parameter :: N_b=1024**2                        !disks number
  real*4, parameter :: eta=0.708D0                         !volume fraction
  real*4, parameter,dimension(1:2) :: Box=(/1.D0,1.D0/)    !Box's size
  real*4, parameter :: r=sqrt(eta*Box(1)*Box(2)/(N_b*pi))  !radius
end module phys_parameters

module alg_parameters
  use phys_parameters
  implicit none
  logical, parameter :: comp_pos=.false.
  logical            :: comp_psi=.false.
  logical, parameter :: comp_pressure=.false.
  logical, parameter :: comp_orientation=.false.
  logical, parameter :: comp_translation=.false.
  logical, parameter :: comp_mu=.false.
  logical, parameter :: comp_track=.true.
  integer*8 :: N_pair=10 !(N_correl=N_pair*N_b)
  real*4, parameter :: lambda_0=0.0768D0/sqrt(real(N_b))  !expected mean freepath, used to compute L and N_extr
  real*4, parameter :: n_dpl=3.125D0*sqrt(real(N_b))      !L=n_dpl*lambda_0 (n_dpl ~ <# displ.>, 10% error for big L)
  real*4 :: L=lambda_0*n_dpl                              !total displacement ell
  integer*8 :: N_dpl_tot=1D08                         !total # displ. for the run (-1 for infinity)
  integer*8,parameter :: N_dpl_extr=1D08              !extraction rate of configurations in number of displacement, writes in 'pos.dat'
  integer*8,parameter :: N_dpl_extr_correl=1D03*N_b      
  integer*8,parameter :: N_dpl_extr_p=5D01*N_b
  integer*8,parameter :: N_dpl_extr_mu=5D00*N_b
  integer*8,parameter :: N_dpl_extr_track=1D03*N_b
  integer*8,parameter :: N_extr=(N_dpl_extr-1)/(n_dpl+1)+1
  integer*8,parameter :: N_extr_correl=(N_dpl_extr_correl-1)/(n_dpl+1)+1
  integer*8,parameter :: N_extr_p=(N_dpl_extr_p-1)/(n_dpl+1)+1
  integer*8,parameter :: N_extr_mu=(N_dpl_extr_mu-1)/(n_dpl+1)+1
  integer*8,parameter :: N_extr_track=(N_dpl_extr_track-1)/(n_dpl+1)+1
  integer*8, dimension(1:N_extr/N_extr_correl+1) :: T_block
end module alg_parameters

module pot_chim
  implicit none
  real*4 :: P_inser=0
  integer*8 :: N_attempt=0
  integer*8 :: N_mu=1
end module pot_chim

module pressure
  implicit none
  integer*8,parameter :: N_smpl_g=100
  integer*8, dimension(1:N_smpl_g) :: N_g=0
  integer*8 :: N_conf_g=0
  real*4 :: r_max_g=0.1D0
end module pressure

module cell_parameters
  use phys_parameters
  integer*8, parameter,dimension(1:2) :: N_cell=(/1,1/)*nint(sqrt(real(N_b))*7.D0/8.D0)   
  real*4, parameter,dimension(1:2) :: l_cell=Box/N_cell                                    
  integer*8, parameter :: n_cell_max=floor((l_cell(1)/r+2.D0)*(l_cell(2)/r+2.D0)/(2.D0*sqrt(3.D0)))
  real*4, dimension(1:2,1:n_cell_max,1:N_cell(1)*N_cell(2)) :: Cell
  integer*2, dimension(1:N_cell(1)*N_cell(2)) :: Cell_ocp=0
  integer*4,dimension(9,N_cell(1)*N_cell(2)) :: Cell_neighb
  integer*8,dimension(1:2) :: i_cell
end module cell_parameters

module track_parameters
  use cell_parameters
  implicit none 
  integer*8, dimension(1:n_cell_max,1:N_cell(1)*N_cell(2)) :: CelltoN=0
  integer*8, dimension(1:2,1:N_b) :: NtoCell=0
  integer*8, parameter :: N_track=1D3
  integer*8, parameter :: T_track=1D3
  real*8, dimension(1:N_track,1:T_track) :: dx_track=0
  integer*8, dimension(1:N_track) :: cur_track=1 !sur le premier
  real*8, dimension(1:N_track) :: x_track=0
  real*8, dimension(1:T_track) :: dx_correl=0
  integer*8, dimension(1:T_track) :: N_dx_correl=0
  real*8, dimension(1:2) :: x_g=0
  real*8 :: dx_av
end module track_parameters

module order_param
  use phys_parameters
  use  alg_parameters
  implicit none
  real*4, dimension(1:2,1:N_extr/N_extr_correl+1) :: PSI_block
  real*4, dimension(1:2) :: PSI
  real*4, dimension(1:2,1:N_b) :: psi_k
  real*4 :: PSI_2
  integer*8, parameter, dimension(1:2) :: N_cell_o=(/1,1/)*nint(sqrt(real(N_b))/6)
  real*4,parameter, dimension(1:2) :: l_cell_o=Box/N_cell_o
  integer*8,parameter :: n_cell_max_o=floor((l_cell_o(1)/r+2.D0)*(l_cell_o(2)/r+2.D0)/(2.D0*sqrt(3.D0)))
  real*4, dimension(1:2,1:n_cell_max_o,1:N_cell_o(1)*N_cell_o(2)) :: Cell_o
  integer*8, dimension(1:N_cell_o(1)*N_cell_o(2)) :: Cell_ocp_o
  integer*8, dimension(9,N_cell_o(1)*N_cell_o(2)) :: Cell_neighb_o
  integer*8 ::  neighb_max=20
end module order_param

module or_parameters
  use phys_parameters
  implicit none
  integer*8, parameter,dimension(1:2) :: N_smpl_or=(/1,1/)*nint(sqrt(real(N_b)))
  real*4, dimension(0:N_smpl_or(1),0:N_smpl_or(2)) :: Correl_or=0
  integer*8, dimension(0:N_smpl_or(1),0:N_smpl_or(2)) :: N_or=0
  real*4, dimension(1:2) :: d_r_or
end module or_parameters

module tra_parameters
  implicit none
  real*4 :: rho_max_tra
  real*4, dimension(1:2) :: d_r_tra
  integer*8, dimension(1:2) :: N_smpl_tra
  integer*8, allocatable, dimension(:,:) :: N_tra
  real*4, allocatable, dimension(:,:) :: g_r_tra
  integer*8 :: N_conf_tra
end module tra_parameters

module others
  use phys_parameters
  real*4 :: time,t_0
  integer*8 :: T
  real*4, dimension(1:2,1:N_b) :: X
end module others

program event_chain
  use order_param
  use phys_parameters
  use cell_parameters
  use alg_parameters
  use others
  use pressure
  use track_parameter
  implicit none 
  real*4 :: L_cur,L_max,L_min
  real*4 :: dx,dy,g
  logical :: ok
  integer*8 :: i,j,k,k_cur,k_min,cell_cur,cell_min,disk,conf
  L_max=min(minval(Box)/2,minval(l_cell))-2*r
  if (L_max<0) then 
     write(*,*) "L_max<0, too much cells" 
     goto 9
  endif
  if (comp_orientation .or. comp_translation) comp_psi=.true.
  T=0
  t_0=0
  r_max_g=r_max_g*r
  call init_cell
  if (comp_psi) call init_psi
  if (comp_orientation)  call init_or
  if (comp_translation)  call init_tra
  write(*,'("N_b=",i7,"  Eta=",d10.4,"  L=",d8.2)') N_b,eta,L
  write(*,'("r=",d10.4"  Box=",2d11.4,"  N_cell=",2i4,"  L_max/r=",d10.4)') r,Box,N_cell,L_max/r
  write(*,*)
  i=0
  ok=.true.
  conf=0
  do
     i=i+1 !new move
     call pick_move(cell_cur,k_cur)      !initial disk
     if (T>N_dpl_tot .and. N_dpl_tot/=-1) exit !condition to end the program
     L_cur=L
     ok=(.not. ok)
     j=0
     if (ok) then  !disk k_cur goes in (ox) direction
        do   !loop over collisions
           j=j+1
           L_min=L
           call expl_cell_x(cell_cur,k_cur,L_min,k_min,cell_min)
           L_min=max(L_min,0.D0)
           dx=Cell(1,k_cur,cell_cur)+min(L_min,L_cur,L_max)
           call refresh_cell_x(dx,cell_cur,k_cur,k_min,cell_min)
           if (L_max<min(L_min,L_cur)) then
              L_cur=L_cur-L_max
              j=j-1
              cycle !Lmax<L_min the move continues with the same disk and the same direction
           elseif (L_min<L_cur) then
              L_cur=L_cur-L_min
              k_cur=k_min
              cell_cur=cell_min
              cycle !Lcur>0, the move continues
           else
              exit !L_cur=0, end of move
           endif
        enddo
     else   !disk k_cur goes in (oy) direction
        do
           j=j+1
           L_min=L
           call expl_cell_y(cell_cur,k_cur,L_min,k_min,cell_min)
           L_min=max(L_min,0.D0)
           dy=Cell(2,k_cur,cell_cur)+min(L_min,L_cur,L_max)
           call refresh_cell_y(dy,cell_cur,k_cur,k_min,cell_min)
           if (L_max<min(L_min,L_cur)) then
              L_cur=L_cur-L_max
              j=j-1
              cycle
           elseif (L_min<L_cur) then
              L_cur=L_cur-L_min
              k_cur=k_min
              cell_cur=cell_min
              cycle
           else
              exit
           endif
        enddo
     endif
     T=T+j
     !*** extraction of observables***!
     if (comp_psi) then
        if (mod(i,N_extr_correl)==0) then
           disk=0
           do j=1,N_cell(1)*N_cell(2)
              do k=1,Cell_ocp(j)
                 disk=disk+1
                 X(1,disk)=Cell(1,k,j)-box(1)/2+(modulo(j-1,N_cell(1))+0.5)*l_cell(1)
                 X(2,disk)=Cell(2,k,j)-box(2)/2+(ceiling(real(j)/N_cell(1))-0.5)*l_cell(2)
              enddo
           enddo
           call order_parameter
           conf=conf+1
           PSI_block(:,conf)= PSI
           T_block(conf)=T
           if (comp_orientation) call fill_or
           if (comp_translation) call fill_tra
        endif
     endif

     if (comp_track) then
        if (mod(i,N_extr_track)==0) then 
           dx_av=L*i/T
           call write_track
        endif
     endif

     if (comp_pressure) then
        if (mod(i,N_extr_p)==0) then 
           call fill_g_pres
        endif
     endif

     if (comp_mu) then
        if (mod(i,N_extr_mu)==0) then 
           call fill_mu
        endif
     endif

     if (mod(i,N_extr)==0) then
        call cpu_time(time)
        time=time+t_0
        write(*,'("T=",i14,"    CPU-time=",f12.2)') T,time
        if (comp_pos) call write_pos
        if (comp_orientation) call write_or
        if (comp_translation) call write_tra
        if (comp_psi) then
           open(1,file='psi.dat',position='append')
           do j=1,conf
              write(1,*) PSI_block(:,j),T_block(j)
           enddo
           close(1)
           conf=0
        endif
        if (comp_pressure) call write_g_pres
        if (comp_mu) call write_mu
     endif
     !***end extraction ***!
  enddo
9 continue
  call cpu_time(time)
  time=time+t_0
  write(*,'("T=",i14,"    CPU-time=",f12.2)') T,time
end program event_chain

subroutine write_track
  use track_parameters
  implicit none
  open(1,file='correl_track.dat',position='append')
  do dt=1,T_track
     write(1,*) dx_correl(dt)/N_dx_correl(dt)
  =dx_correl(dt)+dx_track(i,cur_track(i))*dx_track(i,pos)
  =N_dx_correl(dt)+1
  
  
  write(1,*)
  close(1)
end subroutine write_track

subroutine fill_mu
  use pot_chim
  use cell_parameters
  use alg_parameters
  implicit none
  real*4, dimension(1:2)::dx,x_cur
  integer*8 :: cell_cur,cell_act
  integer*8 :: i,j
  real*4 :: g
  attempt:do i=1,N_mu*N_b
     call random_number(g)
     cell_cur=floor(g*N_cell(1)*N_cell(2))+1
     call random_number(g)
     x_cur(1)=(g-.5)*l_cell(1)
     call random_number(g)
     x_cur(2)=(g-.5)*l_cell(2)
     cell_act=Cell_neighb(5,cell_cur)
     do j=1,Cell_ocp(cell_act)
        dx=Cell(:,j,cell_act)-x_cur
        if (dx(1)**2+dx(2)**2<4*r**2) cycle attempt
     enddo
     cell_act=Cell_neighb(2,cell_cur)
     do j=1,Cell_ocp(cell_act)
        dx=Cell(:,j,cell_act)-x_cur+(/0.,-l_cell(2)/)
        if (dx(1)**2+dx(2)**2<4*r**2) cycle attempt
     enddo
     cell_act=Cell_neighb(4,cell_cur)
     do j=1,Cell_ocp(cell_act)
        dx=Cell(:,j,cell_act)-x_cur+(/-l_cell(1),0./)
        if (dx(1)**2+dx(2)**2<4*r**2) cycle attempt
     enddo
     cell_act=Cell_neighb(6,cell_cur)
     do j=1,Cell_ocp(cell_act)
        dx=Cell(:,j,cell_act)-x_cur+(/l_cell(1),0./)
        if (dx(1)**2+dx(2)**2<4*r**2) cycle attempt
     enddo
     cell_act=Cell_neighb(8,cell_cur)
     do j=1,Cell_ocp(cell_act)
        dx=Cell(:,j,cell_act)-x_cur+(/0.,l_cell(2)/)
        if (dx(1)**2+dx(2)**2<4*r**2) cycle attempt
     enddo
     cell_act=Cell_neighb(1,cell_cur)
     do j=1,Cell_ocp(cell_act)
        dx=Cell(:,j,cell_act)-x_cur+(/-l_cell(1),-l_cell(2)/)
        if (dx(1)**2+dx(2)**2<4*r**2) cycle attempt
     enddo
     cell_act=Cell_neighb(3,cell_cur)
     do j=1,Cell_ocp(cell_act)
        dx=Cell(:,j,cell_act)-x_cur+(/l_cell(1),-l_cell(2)/)
        if (dx(1)**2+dx(2)**2<4*r**2) cycle attempt
     enddo
     cell_act=Cell_neighb(7,cell_cur)
     do j=1,Cell_ocp(cell_act)
        dx=Cell(:,j,cell_act)-x_cur+(/-l_cell(1),l_cell(2)/)
        if (dx(1)**2+dx(2)**2<4*r**2) cycle attempt
     enddo
     cell_act=Cell_neighb(9,cell_cur)
     do j=1,Cell_ocp(cell_act)
        dx=Cell(:,j,cell_act)-x_cur+(/l_cell(1),l_cell(2)/)
        if (dx(1)**2+dx(2)**2<4*r**2) cycle attempt
     enddo     
     P_inser=P_inser+1
  enddo attempt
  N_attempt=N_attempt+N_mu*N_b
end subroutine fill_mu

subroutine write_mu
  use pot_chim
  use phys_parameters
  use others
  implicit none
  open(1,file='mu.dat')
  write(1,*) '# Probabiliy of insertion'
  write(1,*) '# lenght unit is r'
  write(1,*) '# eta=', eta
  write(1,*) '# Box=', Box
  write(1,*) '# N_b=', N_b
  write(1,*) '# r=', r
  write(1,*) '# T=',real(T) 
  write(1,*) '# N_attempt=',N_attempt 
  write(1,*) '# N_insert=',P_inser 
  write(1,*) P_inser/N_attempt
  close(1)
end subroutine write_mu

subroutine write_g_pres
  use phys_parameters
  use pressure
  use others
  implicit none
  integer*8 :: i
  real*4 :: r_1,r_2,dS
  real*4 :: d_r_g
  d_r_g=r_max_g/N_smpl_g
  open(1,file='g_r_pres.dat')
  write(1,*) '# Short distance correlation function'
  write(1,*) '# lenght unit is r'
  write(1,*) '# eta=', eta
  write(1,*) '# Box=', Box
  write(1,*) '# N_b=', N_b
  write(1,*) '# r=', r
  write(1,*) '# dr=', d_r_g
  write(1,*) '# N_smpl_g=', N_smpl_g
  write(1,*) '# T=',real(T) 
  do i=1,N_smpl_g
     r_1=real(i-1)*d_r_g+2*r
     r_2=real(i)*d_r_g+2*r
     dS=(r_2**2-r_1**2)/r**2
     write(1,*) (real(i)-.5D0)*d_r_g/r,real(N_g(i))/(N_b*N_conf_g)/(eta*dS)
  enddo
  close(1)
end subroutine write_g_pres

subroutine fill_g_pres
  use phys_parameters
  use cell_parameters
  use pressure
  implicit none
  integer*8 :: i,j,k_main,k_expl,cell_main,cell_expl,pos
  real*4, dimension(1:2) :: x_main,x_expl,dx
  real*4 :: dist
  N_conf_g=N_conf_g+1
  do cell_main=1,N_cell(1)*N_cell(2)
     do k_main=1,Cell_ocp(cell_main)
        x_main=Cell(:,k_main,cell_main)
        do j=-1,1
           do i=-1,1
              cell_expl=Cell_neighb(i+3*j+5,cell_main)
              do k_expl=1,Cell_ocp(cell_expl)
                 x_expl=Cell(:,k_expl,cell_expl)
                 dx=x_expl-x_main
                 dx(1)=dx(1)+i*l_cell(1)
                 dx(2)=dx(2)+j*l_cell(2)
                 dist=sqrt(dx(1)**2+dx(2)**2)-2*r
                 if (dist>=r_max_g .or. dist<-r) cycle
                 pos=floor(dist/r_max_g*N_smpl_g)+1
                 if (pos<1) pos=1
                 N_g(pos)=N_g(pos)+1
              enddo
           enddo
        enddo
     enddo
  enddo
end subroutine fill_g_pres

subroutine fill_tra
  use tra_parameters
  use order_param
  use others
  use alg_parameters
  implicit none
  real*4 :: alpha,angle_param,g,correl,theta,rho
  integer*8 :: i,j,k,rho_int,theta_int
  real*4,dimension(1:2) :: dx
  N_conf_tra=N_conf_tra+1
  alpha=angle_param(PSI(1),PSI(2)) !-pi/2,3*pi/2
  alpha=(alpha-pi)/6  !-pi/4,pi/12
  do i=1,N_b
     k=0
     do
        call random_number(g)
        j=floor(g*N_b)+1
        if (i==j) cycle
        k=k+1
        if (k==N_pair+1) exit
        dx=X(:,j)-X(:,i)
        call inbox(dx)
        rho=sqrt(dx(1)**2+dx(2)**2)
        if (rho>=rho_max_tra) cycle
        theta=atan(dx(2)/dx(1))-alpha
        if (theta<0.D0) theta=theta+2*pi
        rho_int=floor(rho/rho_max_tra*N_smpl_tra(1))+1     !rho app [0,rho_max[
        theta_int=floor(theta/pi*N_smpl_tra(2)*6)+1    !theta app [0,2*pi[
        N_tra(theta_int,rho_int)=N_tra(theta_int,rho_int)+1
     enddo
  enddo
end subroutine fill_tra

subroutine fill_or
  use or_parameters
  use order_param
  use others
  use alg_parameters
  implicit none
  real*4 :: alpha,angle_param,g,correl
  integer*8 :: i,j,k,x_int,y_int
  real*4,dimension(1:2) :: dx
  alpha=angle_param(PSI(1),PSI(2)) !-pi/2,3*pi/2
  alpha=(alpha-pi)/6  !-pi/4,pi/12
  do i=1,N_b
     k=0
     do
        call random_number(g)
        j=floor(g*N_b)+1
        k=k+1
        if (k==N_pair+1) exit
        dx=X(:,j)-X(:,i)
        call inbox(dx)
        x_int=floor(abs(dx(1))/d_r_or(1))
        y_int=floor(abs(dx(2))/d_r_or(2))
        correl=dot_product(psi_k(:,i),psi_k(:,j))
        Correl_or(x_int,y_int)=Correl_or(x_int,y_int)+correl
        N_or(x_int,y_int)=N_or(x_int,y_int)+1
     enddo
  enddo
end subroutine fill_or

subroutine write_or
  use or_parameters
  use phys_parameters
  implicit none
  integer*8 :: i,j
  real*4,dimension(1:2) :: dx
  open(1,file='pair_correl_sq.dat')
  write(1,*) '# Density distribution of hard spheres'
  write(1,*) '# Orientation distribution of hard spheres'
  write(1,*) '# Cartesian coordinates'
  write(1,*) '# Not oriented'
  write(1,*) '# lenght unit is r'
  write(1,*) '# eta=', eta
  write(1,*) ' #Box=', Box
  write(1,*) ' #N_b=', N_b
  write(1,*) ' #r=', r
  write(1,*) ' #N_smpl_or=', N_smpl_or
  do i=0,N_smpl_or(1)-1
     dx(1)=(real(i)+.5D0)/N_smpl_or(1)*Box(1)/2
     do j=0,N_smpl_or(2)-1
        dx(2)=(real(j)+.5D0)/N_smpl_or(2)*Box(2)/2
        write(1,'(2f8.1,es13.5,es13.5)') dx/r,Correl_or(i,j),real(N_or(i,j))
     enddo
     write(1,*)
  enddo
  close(1)
end subroutine write_or

subroutine write_tra
  use tra_parameters
  use phys_parameters
  use alg_parameters
  implicit none
  integer*8 :: i,j,k,numb
  real*4,dimension(1:2) :: dx
  real*4 ::  rho,theta
  open(1,file='density.dat')
  write(1,*) '# Density distribution of hard spheres'
  write(1,*) '# Polar coordinates'
  write(1,*) '# lenght unit is r'
  write(1,*) '# eta=', eta
  write(1,*) ' #Box=', Box
  write(1,*) ' #N_b=', N_b
  write(1,*) ' #r=', r
  write(1,*) ' #N_smpl=', N_smpl_tra
  g_r_tra=0
  do i=1,N_smpl_tra(1)
     rho=(real(i)-0.5D0)*d_r_tra(1)
     do j=1,N_smpl_tra(2)
        theta=(real(j)-0.5D0)*d_r_tra(2)
        numb=N_tra(j,i)
        do k=1,5
           numb=numb+N_tra(j+k*2*N_smpl_tra(2),i)
        enddo
        do k=1,6
           numb=numb+N_tra(2*k*N_smpl_tra(2)+1-j,i)
        enddo
        g_r_tra(j,i)=real(numb*(N_b-1))/(rho*d_r_tra(1)*d_r_tra(2)*N_conf_tra*N_pair*N_b*N_b)/12
!        write(1,'(f8.2,f11.4,f8.4)') rho/r*cos(theta), rho/r*sin(theta), g_r_tra(j,i)
        write(1,'(f6.2)') g_r_tra(j,i)
     enddo
     write(1,*)
  enddo
  close(1)
end subroutine write_tra

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
  use track_parameters
  use others
  implicit none
  integer*8 :: i,i_2,j_2,j,c_1,c_2,n_1,n_2
  logical :: ex
  character ( len = 50 ) :: line
  inquire(file='pos.dat', exist=ex)
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
     if (comp_pos) then
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
     endif
     if (comp_psi) then
        open(1,file='psi.dat')
        write(1,*) '# Orientational order parameter'
        write(1,*) '# Vornoi construction'
        write(1,*) '# Box=',Box
        write(1,*) '# N_b=',N_b
        write(1,*) '# r=',r
        write(1,*) '# eta=',r**2*N_b*acos(-1.)/Box(1)/Box(2)
        close(1)
     endif
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
     NtoCell(:,i)=(/n_1,Cell_ocp(n_1)/)
     CelltoN(n_1,Cell_ocp(n_1))=i
  enddo
  !*** end fill cells ***!
end subroutine init_cell

subroutine init_tra
  use phys_parameters
  use tra_parameters
  implicit none
  rho_max_tra=minval(Box)/2
  d_r_tra(1)=r/5
  d_r_tra(2)=d_r_tra(1)/rho_max_tra
  N_smpl_tra(1)=floor(rho_max_tra/d_r_tra(1))
  d_r_tra(1)=rho_max_tra/N_smpl_tra(1)
  N_smpl_tra(2)=floor(pi/6/d_r_tra(2))
  d_r_tra(2)=pi/6/N_smpl_tra(2)
  allocate(N_tra(1:N_smpl_tra(2)*12,1:N_smpl_tra(1)))
  allocate(g_r_tra(1:N_smpl_tra(2),1:N_smpl_tra(1)))
  N_tra=0
  N_conf_tra=0
end subroutine init_tra

subroutine init_psi
  use order_param
  implicit none
  integer*8, dimension(1:2) :: i_cell_o 
  integer*8 :: i,j,c_1,c_2,n_1,n_2
 !*** create the neighbor table ***!
  do j=1,N_cell_o(2)
     do i=1,N_cell_o(1)
        do c_2=-1,1
           do c_1=-1,1
              n_1=i+(j-1)*N_cell_o(1)
              i_cell_o(1)=modulo(i+c_1-1,N_cell_o(1))+1
              i_cell_o(2)=modulo(j+c_2-1,N_cell_o(2))+1
              n_2=i_cell_o(1)+(i_cell_o(2)-1)*N_cell_o(1)
              Cell_neighb_o(c_1+2+(c_2+1)*3,n_1)=n_2
           enddo
        enddo
     enddo
  enddo
  !*** end neighbor table ***!
end subroutine init_psi

subroutine init_or
  use or_parameters
  use order_param
  implicit none
  d_r_or(1)=Box(1)/2.D0/N_smpl_or(1)
  d_r_or(2)=Box(2)/2.D0/N_smpl_or(2)
 end subroutine init_or

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
        write(1,'(f8.5,f9.5)') Cell(1,k,j)-box(1)/2+(modulo(j-1,N_cell(1))+0.5)*l_cell(1),&
             & Cell(2,k,j)-box(2)/2+((j-1)/N_cell(1)+0.5)*l_cell(2)
     enddo
  enddo
  close(1)
end subroutine write_pos

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
  real*4, dimension(1:2) ::dx=(/1.D0,0.D0/),dy=(/0.5D0,1.D0/)
  logical :: ex
  inquire(file='init.dat', exist=ex)
  if (ex) then
     open(1,file='init.dat')
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

subroutine pick_move(cell_cur,k_cur)
!******************************************************************!
!
!Picks one disk in the cells. equi-probability for each disk
!
!******************************************************************!
  use alg_parameters
  use phys_parameters
  use cell_parameters
  implicit none
  integer*8 :: k_cur,cell_cur
  real*4 ::g
  do
     call random_number(g)
     cell_cur=floor(g*N_cell(1)*N_cell(2))+1
     call random_number(g)
     k_cur=floor(g*n_cell_max)+1
     if (k_cur>cell_ocp(cell_cur)) cycle
     exit
  enddo
end subroutine pick_move

subroutine expl_cell_x(cell_cur,k_cur,L_min,k_min,cell_min)
!******************************************************************!
! 
! Explore the disks neighbor of k_cur, 6 cell-neighb to explore (2,3,5,6,8,9)
! Faster without the loop
!
!******************************************************************!
  use phys_parameters
  use cell_parameters
  implicit none
  real*4, dimension(1:2)::dx,x_cur
  integer*8 :: cell_cur,cell_act,cell_min
  real*4 :: L_min,L_coll
  integer*8 :: i,j,k_cur,k_min
  x_cur=Cell(:,k_cur,cell_cur)
  cell_act=Cell_neighb(6,cell_cur)
  do j=1,Cell_ocp(cell_act)
     dx(2)=Cell(2,j,cell_act)-x_cur(2)
     if (abs(dx(2))<2*r) then
        dx(1)=Cell(1,j,cell_act)-x_cur(1)+l_cell(1)
        L_coll=dx(1)-sqrt(4*r**2-dx(2)**2)
        if (L_coll<L_min) then
           L_min=L_coll
           k_min=j
           cell_min=cell_act
        endif
     endif
  enddo
  cell_act=Cell_neighb(2,cell_cur)
  do j=1,Cell_ocp(cell_act)
     dx(2)=Cell(2,j,cell_act)-x_cur(2)-l_cell(2)
     if (abs(dx(2))<2*r) then
        dx(1)=Cell(1,j,cell_act)-x_cur(1)
        if (dx(1)>0.D0) then
           L_coll=dx(1)-sqrt(4*r**2-dx(2)**2)
           if (L_coll<L_min) then
              L_min=L_coll
              k_min=j
              cell_min=cell_act
           endif
        endif
     endif
  enddo
  cell_act=Cell_neighb(8,cell_cur)
  do j=1,Cell_ocp(cell_act)
     dx(2)=Cell(2,j,cell_act)-x_cur(2)+l_cell(2)
     if (abs(dx(2))<2*r) then
        dx(1)=Cell(1,j,cell_act)-x_cur(1)
        if (dx(1)>0.D0) then
           L_coll=dx(1)-sqrt(4*r**2-dx(2)**2)
           if (L_coll<L_min) then
              L_min=L_coll
              k_min=j
              cell_min=cell_act
           endif
        endif
     endif
  enddo
  cell_act=Cell_neighb(5,cell_cur)
  do j=1,Cell_ocp(cell_act)
     dx(2)=Cell(2,j,cell_act)-x_cur(2)
     if (abs(dx(2))<2*r) then
        dx(1)=Cell(1,j,cell_act)-x_cur(1)
        if (dx(1)>0.D0) then
           L_coll=dx(1)-sqrt(4*r**2-dx(2)**2)
           if (L_coll<L_min) then
              L_min=L_coll
              k_min=j
              cell_min=cell_act
           endif
        endif
     endif
  enddo
  cell_act=Cell_neighb(3,cell_cur)
  do j=1,Cell_ocp(cell_act)
     dx(2)=Cell(2,j,cell_act)-x_cur(2)-l_cell(2)
     if (abs(dx(2))<2*r) then
        dx(1)=Cell(1,j,cell_act)-x_cur(1)+l_cell(1)
        L_coll=dx(1)-sqrt(4*r**2-dx(2)**2)
        if (L_coll<L_min) then
           L_min=L_coll
           k_min=j
           cell_min=cell_act
        endif
     endif
  enddo  
  cell_act=Cell_neighb(9,cell_cur)
  do j=1,Cell_ocp(cell_act)
     dx(2)=Cell(2,j,cell_act)-x_cur(2)+l_cell(2)
     if (abs(dx(2))<2*r) then
        dx(1)=Cell(1,j,cell_act)-x_cur(1)+l_cell(1)
        L_coll=dx(1)-sqrt(4*r**2-dx(2)**2)
        if (L_coll<L_min) then
           L_min=L_coll
           k_min=j
           cell_min=cell_act
        endif
     endif
  enddo
end subroutine expl_cell_x

subroutine expl_cell_y(cell_cur,k_cur,L_min,k_min,cell_min)
!******************************************************************!
! 
! Explore the disks neighbor of k_cur, 6 cell-neighb to explore (4,5,6,7,8,9)
! Faster without the loop
!
!******************************************************************!
  use phys_parameters
  use cell_parameters
  implicit none
  real*4, dimension(1:2)::dx,x_cur
  integer*8 :: cell_cur,cell_act,cell_min
  real*4 :: L_min,L_coll
  integer*8 :: i,j,k_cur,k_min
  x_cur=Cell(:,k_cur,cell_cur)  
  cell_act=Cell_neighb(8,cell_cur)
  do j=1,Cell_ocp(cell_act)
     dx(1)=Cell(1,j,cell_act)-x_cur(1)       
     if (abs(dx(1))<2*r) then
        dx(2)=Cell(2,j,cell_act)-x_cur(2)+l_cell(2)
        L_coll=dx(2)-sqrt(4*r**2-dx(1)**2)
        if (L_coll<L_min) then
           L_min=L_coll
           k_min=j
           cell_min=cell_act
        endif
     endif
  enddo
  cell_act=Cell_neighb(6,cell_cur)
  do j=1,Cell_ocp(cell_act)
     dx(1)=Cell(1,j,cell_act)-x_cur(1)+l_cell(1)       
     if (abs(dx(1))<2*r) then
        dx(2)=Cell(2,j,cell_act)-x_cur(2)
        if (dx(2)>0.D0) then
           L_coll=dx(2)-sqrt(4*r**2-dx(1)**2)
           if (L_coll<L_min) then
              L_min=L_coll
              k_min=j
              cell_min=cell_act
           endif
        endif
     endif
  enddo
  cell_act=Cell_neighb(4,cell_cur)
  do j=1,Cell_ocp(cell_act)
     dx(1)=Cell(1,j,cell_act)-x_cur(1)-l_cell(1)       
     if (abs(dx(1))<2*r) then
        dx(2)=Cell(2,j,cell_act)-x_cur(2)
        if (dx(2)>0.D0) then
           L_coll=dx(2)-sqrt(4*r**2-dx(1)**2)
           if (L_coll<L_min) then
              L_min=L_coll
              k_min=j
              cell_min=cell_act
           endif
        endif
     endif
  enddo
  cell_act=Cell_neighb(5,cell_cur)
  do j=1,Cell_ocp(cell_act)
     dx(1)=Cell(1,j,cell_act)-x_cur(1)       
     if (abs(dx(1))<2*r) then
        dx(2)=Cell(2,j,cell_act)-x_cur(2)
        if (dx(2)>0.D0) then           
           L_coll=dx(2)-sqrt(4*r**2-dx(1)**2)
           if (L_coll<L_min) then
              L_min=L_coll
              k_min=j
              cell_min=cell_act
           endif
        endif
     endif
  enddo
  cell_act=Cell_neighb(9,cell_cur)
  do j=1,Cell_ocp(cell_act)
     dx(1)=Cell(1,j,cell_act)-x_cur(1)+l_cell(1)   
     if (abs(dx(1))<2*r) then
        dx(2)=Cell(2,j,cell_act)-x_cur(2)+l_cell(2)
        L_coll=dx(2)-sqrt(4*r**2-dx(1)**2)
        if (L_coll<L_min) then
           L_min=L_coll
           k_min=j
           cell_min=cell_act
        endif
     endif
  enddo  
  cell_act=Cell_neighb(7,cell_cur)
  do j=1,Cell_ocp(cell_act)
     dx(1)=Cell(1,j,cell_act)-x_cur(1)-l_cell(1)
     if (abs(dx(1))<2*r) then
        dx(2)=Cell(2,j,cell_act)-x_cur(2)+l_cell(2)
        L_coll=dx(2)-sqrt(4*r**2-dx(1)**2)
        if (L_coll<L_min) then
           L_min=L_coll
           k_min=j
           cell_min=cell_act
        endif
     endif
  enddo
end subroutine expl_cell_y

subroutine refresh_cell_x(dx,cell_cur,k_cur)
!******************************************************************!
! 
! Compute the new cell of k_cur
!
!******************************************************************!
  use phys_parameters
  use cell_parameters
  use track_parameters
  implicit none
  real*4 :: dx
  real*8 :: correl
  integer*8 :: cell_new,cell_cur,k_cur,i,j,dt,pos

  i=CelltoN(k_cur,cell_cur)
  if (i < = N_track) then
     x_track(:,i)=x_track(:,i)+(/dx,0/)
     dx_track(i,cur_track(i))=dx
     do dt=1,T_track-1
        pos=cur_track(i)-dt
        if (pos<1) pos=pos+T_track
        dx_correl(dt)=dx_correl(dt)+dx_track(i,cur_track(i))*dx_track(i,pos)
        N_dx_correl(dt)=N_dx_correl(dt)+1
        cur_track(i)=cur_track(i)+1
        if (cur_track(i)>T_track) cur_track(i)=1
     enddo
  endif
  if (dx>l_cell(1)/2) then                                      !disk k_cur is now in cell_new=Cell_neighb(6,cell_cur)

     dx=dx-l_cell(1)                                            ! position in the new cell
     cell_new=Cell_neighb(6,cell_cur)                           !new cell
     Cell_ocp(cell_new)=Cell_ocp(cell_new)+1                    ! refresh # of disk in the new cell
     Cell(1,Cell_ocp(cell_new),cell_new)=dx                     ! x position in the new cell
     Cell(2,Cell_ocp(cell_new),cell_new)=Cell(2,k_cur,cell_cur) ! y position in the new cell

     !track i
     NtoCell(:,i)=(/Cell_ocp(cell_new),cell_new/) 
     CelltoN(Cell_ocp(cell_new),cell_new)=i       
     !end track i

     !swap
     Cell(:,k_cur,cell_cur)=Cell(:,Cell_ocp(cell_cur),cell_cur) !swap the last disk in the old cell to take k_cur's vacant place 

     !track j
     j=CelltoN(Cell_ocp(cell_cur),cell_cur) !track
     NtoCell(:,j)=(/k_cur,cell_cur/) !track disk
     CelltoN(k_cur,cell_cur)=j       !track disk
     !end track j

     Cell_ocp(cell_cur)=Cell_ocp(cell_cur)-1                    ! refresh # of disks in the old cell
     !end swap

     cell_cur=cell_new                                          ! The cell of the "current" disk is now cell_new....
     k_cur=Cell_ocp(cell_new)                                   ! and its place is the last one


  else                                                          !disk k_cur stays in cell_cur
     Cell(1,k_cur,cell_cur)=dx
  endif
end subroutine refresh_cell_x

subroutine refresh_cell_y(dy,cell_cur,k_cur)
!******************************************************************!
! 
! Compute the new cell of k_cur
!
!******************************************************************!
  use phys_parameters
  use cell_parameters
  implicit none
  real*4 :: dy
  integer*8:: cell_new,cell_cur,k_cur

  i=CelltoN(k_cur,cell_cur)
  if (i < = N_track) then
     x_track(:,i)=x_track(:,i)+(/0,dx/)
  endif

  if (dy>l_cell(2)/2) then
     dy=dy-l_cell(2)
     cell_new=Cell_neighb(8,cell_cur)
     Cell_ocp(cell_new)=Cell_ocp(cell_new)+1 
     Cell(1,Cell_ocp(cell_new),cell_new)=Cell(1,k_cur,cell_cur)
     Cell(2,Cell_ocp(cell_new),cell_new)=dy

     !track i
     NtoCell(:,i)=(/Cell_ocp(cell_new),cell_new/) 
     CelltoN(Cell_ocp(cell_new),cell_new)=i       
     !end track i

     !swap
     Cell(:,k_cur,cell_cur)=Cell(:,Cell_ocp(cell_cur),cell_cur) !swap the last disk in the old cell to take k_cur's vacant place 

     !track j
     j=CelltoN(Cell_ocp(cell_cur),cell_cur) !track
     NtoCell(:,j)=(/k_cur,cell_cur/) !track disk
     CelltoN(k_cur,cell_cur)=j       !track disk
     !end track j

     Cell_ocp(cell_cur)=Cell_ocp(cell_cur)-1
     cell_cur=cell_new
     k_cur=Cell_ocp(cell_new)
  else
     Cell(2,k_cur,cell_cur)=dy
  endif
end subroutine refresh_cell_y

subroutine order_parameter
!***********************************************!
!
! Computes orientational order parameter
! 
! Uses Voronoi contruction
!
! Uses Cells => O(N)
!
!***********************************************!
  use or_parameters
  use order_param
  use others
  implicit none
  integer   ( kind = 4 )  node_num
  real      ( kind = 8 ), allocatable, dimension ( :, : ) :: node_xy
  integer   ( kind = 4 ), allocatable, dimension ( :, : ) :: triangle_neighbor
  integer   ( kind = 4 ), allocatable, dimension ( :, : ) :: triangle_node
  integer   ( kind = 4 )  triangle_num
  integer   ( kind = 4 )  triangle_order
  integer*8, allocatable, dimension(:) :: N_neighbor
  integer*8, allocatable, dimension(:,:) :: neighbor
  real*4, dimension(2) :: dx
  integer*8, dimension(1:2) :: i_cell_c 
  integer*8 :: s_1,s_2,s_3,c,i,j,k,disk
  complex*8 :: c_dummy, c2_dummy, xi=(0.D0,1.D0)
  real*4 :: angle
  integer*8 :: c_1,c_2,n_1,n_2,cell_cur
  
  !*** Fill Cells***!
  Cell_ocp_o=0
  Cell_o=0
  do i=1,N_b
     do j=1,2
        i_cell_c(j)=floor((X(j,i)+Box(j)/2)/l_cell_o(j))+1
        if (i_cell_c(j)==N_cell_o(j)+1) then
           i_cell_c(j)=1
           X(j,i)=X(j,i)-Box(j)
        endif
     enddo
     n_1=i_cell_c(1)+(i_cell_c(2)-1)*N_cell_o(1)
     Cell_ocp_o(n_1)=Cell_ocp_o(n_1)+1
     do j=1,2   
        Cell_o(j,Cell_ocp_o(n_1),n_1)=X(j,i)+Box(j)/2-(real(i_cell_c(j))-0.5D0)*l_cell_o(j)
     enddo
  enddo
  !*** End fill Cells***!
  
  PSI=0.D0
  PSI_2=0.D0
  disk=0
  do cell_cur=1,N_cell_o(1)*N_cell_o(2)
     !*** fill node_xy ***!
     node_num=0
     do i=1,9
        node_num=node_num+Cell_ocp_o(Cell_neighb_o(i,cell_cur))
     enddo
     node_num=node_num+4
     allocate(node_xy(1:2,1:node_num))
     allocate(N_neighbor(1:node_num))
     allocate(neighbor(1:node_num,1:neighb_max))
     N_neighbor=0
     neighbor=0
     node_xy=0
     k=0
     do j=1,Cell_ocp_o(cell_cur)
        k=k+1
        node_xy(:,k)=Cell_o(:,j,cell_cur)
     enddo
     do i=1,4
        do j=1,Cell_ocp_o(Cell_neighb_o(i,cell_cur))
           k=k+1
           node_xy(:,k)=Cell_o(:,j,Cell_neighb_o(i,cell_cur))+&
                & (/modulo(i-1,3)-1,(i-1)/3-1/)*l_cell_o
        enddo
     enddo
     do i=6,9
        do j=1,Cell_ocp_o(Cell_neighb_o(i,cell_cur))
           k=k+1
           node_xy(:,k)=Cell_o(:,j,Cell_neighb_o(i,cell_cur))+&
                &(/modulo(i-1,3)-1,(i-1)/3-1/)*l_cell_o
        enddo
     enddo
     node_xy(:,k+1)=2*(/-l_cell_o(1),-l_cell_o(2)/)
     node_xy(:,k+2)=2*(/-l_cell_o(1),l_cell_o(2)/)
     node_xy(:,k+3)=2*(/l_cell_o(1),-l_cell_o(2)/)
     node_xy(:,k+4)=2*(/l_cell_o(1),l_cell_o(2)/)
     !*** end fill node_xy ***!
     !***  Determine the Delaunay triangulation ***!
     triangle_order = 3
     allocate ( triangle_node(triangle_order,3*node_num) )
     allocate ( triangle_neighbor(triangle_order,3*node_num) )
     triangle_node=0
     triangle_neighbor=0
     triangle_num=0
     call dtris2 ( node_num, node_xy, triangle_num, triangle_node, & 
          triangle_neighbor ) !!pas de tri de node_xy apparement !!!
     do i=1,triangle_num
        s_1=triangle_node(1,i)
        s_2=triangle_node(2,i)
        s_3=triangle_node(3,i)
        N_neighbor(s_1)=min(N_neighbor(s_1)+1,neighb_max)
        N_neighbor(s_2)=min(N_neighbor(s_2)+1,neighb_max)
        N_neighbor(s_3)=min(N_neighbor(s_3)+1,neighb_max)
        neighbor(s_1,N_neighbor(s_1))=s_2
        neighbor(s_2,N_neighbor(s_2))=s_3
        neighbor(s_3,N_neighbor(s_3))=s_1
     enddo
     !***  End determine the Delaunay triangulation ***!
     !*** compute orientationnal order parameter for cell_cur***!
     do i=1,Cell_ocp_o(cell_cur)
        disk=disk+1
        c_dummy=(0.D0,0.D0)
        if (N_neighbor(i) /= 0) then
           do j=1,N_neighbor(i)
              dx=node_xy(:,neighbor(i,j))-node_xy(:,i)
              c2_dummy=6.D0*xi*angle(dx(1),dx(2))
              c_dummy=c_dummy + exp(c2_dummy)
           enddo
           c_dummy=c_dummy/N_neighbor(i)
        else
           c_dummy=(0.D0,0.D0)
        endif
        X(1,disk)= node_xy(1,i)-box(1)/2+(modulo(cell_cur-1,N_cell_o(1))+0.5)*l_cell_o(1)
        X(2,disk)= node_xy(2,i)-box(2)/2+(ceiling(real(cell_cur)/N_cell_o(1))-0.5)*l_cell_o(2)!attention
        psi_k(1,disk)=real(c_dummy)
        psi_k(2,disk)=imag(c_dummy)
        PSI=PSI+psi_k(:,disk)
        PSI_2=PSI_2+psi_k(1,disk)**2+psi_k(2,disk)**2
     enddo
     !*** end compute orientationnal order parameter for cell_cur***!
     deallocate(node_xy)
     deallocate(N_neighbor)
     deallocate(neighbor)
     deallocate(triangle_node)
     deallocate(triangle_neighbor)
  enddo
  PSI=PSI/N_b
  Correl_or(0,0)=Correl_or(0,0)+PSI_2
  N_or(0,0)=N_or(0,0)+N_b
  PSI_2=PSI_2/N_b
end subroutine order_parameter

subroutine inbox(dx)
  use phys_parameters
  implicit none
  real*4, dimension(1:2) :: dx
  integer :: i
  do i=1,2
     dx(i)=mod(dx(i),Box(i))
     if (dx(i)<-Box(i)/2.) dx(i)=dx(i)+Box(i)
     if (dx(i)>Box(i)/2.)  dx(i)=dx(i)-Box(i)
  enddo
end subroutine inbox

real*4 function angle(dx,dy) !-pi/2,pi/2
  use phys_parameters
  implicit none
  real*4 :: dx,dy
  if (dx/=0) then
     angle=atan(dy/dx)
  else
     if (dy>0) then
        angle=pi/2
     else
        angle=-pi/2
     endif
  endif
end function angle

real*4 function angle_param(dx,dy) !
  use phys_parameters
  implicit none
  real*4 :: dx,dy
  if (dx/=0) then
     angle_param=atan(dy/dx)
  else
     if (dy>0) then
        angle_param=pi/2
     else
        angle_param=-pi/2
     endif
  endif
  if (dx<0) angle_param=angle_param+pi
end function angle_param

subroutine ch_cap ( c )

!*****************************************************************************80
!
!! CH_CAP capitalizes a single character.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character C, the character to capitalize.
!
  implicit none

  character              c
  integer   ( kind = 4 ) itemp

  itemp = ichar ( c )

  if ( 97 <= itemp .and. itemp <= 122 ) then
    c = char ( itemp - 32 )
  end if

  return
end
function ch_eqi ( c1, c2 )

!*****************************************************************************80
!
!! CH_EQI is a case insensitive comparison of two characters for equality.
!
!  Example:
!
!    CH_EQI ( 'A', 'a' ) is .TRUE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C1, C2, the characters to compare.
!
!    Output, logical CH_EQI, the result of the comparison.
!
  implicit none

  logical   ch_eqi
  character c1
  character c1_cap
  character c2
  character c2_cap

  c1_cap = c1
  c2_cap = c2

  call ch_cap ( c1_cap )
  call ch_cap ( c2_cap )

  if ( c1_cap == c2_cap ) then
    ch_eqi = .true.
  else
    ch_eqi = .false.
  end if

  return
end
subroutine ch_to_digit ( c, digit )

!*****************************************************************************80
!
!! CH_TO_DIGIT returns the integer value of a base 10 digit.
!
!  Example:
!
!     C   DIGIT
!    ---  -----
!    '0'    0
!    '1'    1
!    ...  ...
!    '9'    9
!    ' '    0
!    'X'   -1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C, the decimal digit, '0' through '9' or blank
!    are legal.
!
!    Output, integer ( kind = 4 ) DIGIT, the corresponding integer value.  
!    If C was 'illegal', then DIGIT is -1.
!
  implicit none

  character              c
  integer   ( kind = 4 ) digit

  if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then

    digit = ichar ( c ) - 48

  else if ( c == ' ' ) then

    digit = 0

  else

    digit = -1

  end if

  return
end
function diaedg ( x0, y0, x1, y1, x2, y2, x3, y3 )

!*****************************************************************************80
!
!! DIAEDG chooses a diagonal edge.
!
!  Discussion:
!
!    The routine determines whether 0--2 or 1--3 is the diagonal edge
!    that should be chosen, based on the circumcircle criterion, where
!    (X0,Y0), (X1,Y1), (X2,Y2), (X3,Y3) are the vertices of a simple
!    quadrilateral in counterclockwise order.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 February 2001
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms, 
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X0, Y0, X1, Y1, X2, Y2, X3, Y3, the
!    coordinates of the vertices of a quadrilateral, given in
!    counter clockwise order.
!
!    Output, integer ( kind = 4 ) DIAEDG, chooses a diagonal:
!    +1, if diagonal edge 02 is chosen;
!    -1, if diagonal edge 13 is chosen;
!     0, if the four vertices are cocircular.
!
  implicit none

  real    ( kind = 8 ) ca
  real    ( kind = 8 ) cb
  integer ( kind = 4 ) diaedg
  real    ( kind = 8 ) dx10
  real    ( kind = 8 ) dx12
  real    ( kind = 8 ) dx30
  real    ( kind = 8 ) dx32
  real    ( kind = 8 ) dy10
  real    ( kind = 8 ) dy12
  real    ( kind = 8 ) dy30
  real    ( kind = 8 ) dy32
  real    ( kind = 8 ) s
  real    ( kind = 8 ) tol
  real    ( kind = 8 ) tola
  real    ( kind = 8 ) tolb
  real    ( kind = 8 ) x0
  real    ( kind = 8 ) x1
  real    ( kind = 8 ) x2
  real    ( kind = 8 ) x3
  real    ( kind = 8 ) y0
  real    ( kind = 8 ) y1
  real    ( kind = 8 ) y2
  real    ( kind = 8 ) y3

  tol = 100.0D+00 * epsilon ( tol )

  dx10 = x1 - x0
  dy10 = y1 - y0
  dx12 = x1 - x2
  dy12 = y1 - y2
  dx30 = x3 - x0
  dy30 = y3 - y0
  dx32 = x3 - x2
  dy32 = y3 - y2

  tola = tol * max ( abs ( dx10 ), abs ( dy10 ), abs ( dx30 ), abs ( dy30 ) )
  tolb = tol * max ( abs ( dx12 ), abs ( dy12 ), abs ( dx32 ), abs ( dy32 ) )

  ca = dx10 * dx30 + dy10 * dy30
  cb = dx12 * dx32 + dy12 * dy32

  if ( tola < ca .and. tolb < cb ) then

    diaedg = -1

  else if ( ca < -tola .and. cb < -tolb ) then

    diaedg = 1

  else

    tola = max ( tola, tolb )
    s = ( dx10 * dy30 - dx30 * dy10 ) * cb + ( dx32 * dy12 - dx12 * dy32 ) * ca

    if ( tola < s ) then
      diaedg = -1
    else if ( s < -tola ) then
      diaedg = 1
    else
      diaedg = 0
    end if

  end if

  return
end
subroutine dtable_data_read ( input_file_name, m, n, table )

!*****************************************************************************80
!
!! DTABLE_DATA_READ reads data from a DTABLE file.
!
!  Discussion:
!
!    The file may contain more than N points, but this routine will
!    return after reading N of them.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILE_NAME, the name of the input file.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Output, real ( kind = 8 ) TABLE(M,N), the table data.
!
  implicit none

  integer   ( kind = 4 )   m
  integer   ( kind = 4 )   n

  integer   ( kind = 4 )   ierror
  character ( len = * )    input_file_name
  integer   ( kind = 4 )   input_status
  integer   ( kind = 4 )   input_unit
  integer   ( kind = 4 )   j
  character ( len = 255 )  line
  real      ( kind = 8 )   table(m,n)
  real      ( kind = 8 )   x(m)

  ierror = 0

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_file_name, status = 'old', &
    iostat = input_status )

  if ( input_status /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DTABLE_DATA_READ - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the input file "' // &
      trim ( input_file_name ) // '" on unit ', input_unit
    stop
  end if

  j = 0

  do while ( j < n )

    read ( input_unit, '(a)', iostat = input_status ) line

    if ( input_status /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DTABLE_DATA_READ - Fatal error!'
      write ( *, '(a)' ) '  Error while reading lines of data.'
      write ( *, '(a,i8)' ) '  Number of values expected per line M = ', m
      write ( *, '(a,i8)' ) '  Number of data lines read, J =         ', j
      write ( *, '(a,i8)' ) '  Number of data lines needed, N =       ', n
      stop
    end if

    if ( line(1:1) == '#' .or. len_trim ( line ) == 0 ) then
      cycle
    end if

    call s_to_r8vec ( line, m, x, ierror )

    if ( ierror /= 0 ) then
      cycle
    end if

    j = j + 1

    table(1:m,j) = x(1:m)

  end do

  close ( unit = input_unit )

  return
end
subroutine dtable_header_read ( input_file_name, m, n )

!*****************************************************************************80
!
!! DTABLE_HEADER_READ reads the header from a DTABLE file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILE_NAME, the name of the input file.
!
!    Output, integer ( kind = 4 ) M, spatial dimension.
!
!    Output, integer ( kind = 4 ) N, the number of points. 
!
  implicit none

  character ( len = * )  input_file_name
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) n

  call file_column_count ( input_file_name, m )

  if ( m <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DTABLE_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  There was some kind of I/O problem while trying'
    write ( *, '(a)' ) '  to count the number of data columns in'
    write ( *, '(a)' ) '  the file "' // trim ( input_file_name ) // '".'
    stop
  end if

  call file_row_count ( input_file_name, n )

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DTABLE_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  There was some kind of I/O problem while trying'
    write ( *, '(a)' ) '  to count the number of data rows in'
    write ( *, '(a)' ) '  the file "' // trim ( input_file_name ) // '".'
    stop
  end if

  return
end
subroutine dtris2 ( point_num, point_xy, tri_num, tri_vert, tri_nabe )

!*****************************************************************************80
!
!! DTRIS2 constructs a Delaunay triangulation of 2D vertices.
!
!  Discussion:
!
!    The routine constructs the Delaunay triangulation of a set of 2D vertices
!    using an incremental approach and diagonal edge swaps.  Vertices are
!    first sorted in lexicographically increasing (X,Y) order, and
!    then are inserted one at a time from outside the convex hull.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 August 2001
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms, 
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of vertices.
!
!    Input/output, real ( kind = 8 ) POINT_XY(2,POINT_NUM), the coordinates 
!    of the vertices.  On output, the vertices have been sorted into 
!    dictionary order.
!
!    Output, integer ( kind = 4 ) TRI_NUM, the number of triangles in the 
!    triangulation; TRI_NUM is equal to 2*POINT_NUM - NB - 2, where NB is the 
!    number of boundary vertices.
!
!    Output, integer ( kind = 4 ) TRI_VERT(3,TRI_NUM), the nodes that make up
!    each triangle.  The elements are indices of POINT_XY.  The vertices of the 
!    triangles are in counter clockwise order.
!
!    Output, integer ( kind = 4 ) TRI_NABE(3,TRI_NUM), the triangle neighbor 
!    list.  Positive elements are indices of TIL; negative elements are used 
!    for links of a counter clockwise linked list of boundary edges; 
!    LINK = -(3*I + J-1) where I, J = triangle, edge index; TRI_NABE(J,I) refers
!    to the neighbor along edge from vertex J to J+1 (mod 3).
!
  implicit none

  integer ( kind = 4 ) point_num

  real    ( kind = 8 ) cmax
  integer ( kind = 4 ) e
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) indx(point_num)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) ledg
  integer ( kind = 4 ) lr
  integer ( kind = 4 ) lrline
  integer ( kind = 4 ) ltri
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) n
  real    ( kind = 8 ) point_xy(2,point_num)
  integer ( kind = 4 ) redg
  integer ( kind = 4 ) rtri
  integer ( kind = 4 ) stack(point_num)
  integer ( kind = 4 ) t
  real    ( kind = 8 ) tol
  integer ( kind = 4 ) top
  integer ( kind = 4 ) tri_nabe(3,point_num*2)
  integer ( kind = 4 ) tri_num
  integer ( kind = 4 ) tri_vert(3,point_num*2)

  tol = 100.0D+00 * epsilon ( tol )

  ierr = 0
!
!  Sort the vertices by increasing (x,y).
!
  call r82vec_sort_heap_index_a ( point_num, point_xy, indx )

  call r82vec_permute ( point_num, indx, point_xy )
!
!  Make sure that the data points are "reasonably" distinct.
!
  m1 = 1

  do i = 2, point_num

    m = m1
    m1 = i

    k = 0

    do j = 1, 2

      cmax = max ( abs ( point_xy(j,m) ), abs ( point_xy(j,m1) ) )

      if ( tol * ( cmax + 1.0D+00 ) &
           < abs ( point_xy(j,m) - point_xy(j,m1) ) ) then
        k = j
        exit
      end if

    end do

    if ( k == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DTRIS2 - Fatal error!'
      write ( *, '(a,i8)' ) '  Fails for point number I = ', i
      write ( *, '(a,i8)' ) '  M = ', m
      write ( *, '(a,i8)' ) '  M1 = ', m1
      write ( *, '(a,2g14.6)' ) '  X,Y(M)  = ', point_xy(1,m), point_xy(2,m)
      write ( *, '(a,2g14.6)' ) '  X,Y(M1) = ', point_xy(1,m1), point_xy(2,m1)
      ierr = 224
      return
    end if

  end do
!
!  Starting from points M1 and M2, search for a third point M that
!  makes a "healthy" triangle (M1,M2,M)
!
  m1 = 1
  m2 = 2
  j = 3

  do

    if ( point_num < j ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DTRIS2 - Fatal error!'
      ierr = 225
      return
    end if

    m = j

    lr = lrline ( point_xy(1,m), point_xy(2,m), point_xy(1,m1), &
      point_xy(2,m1), point_xy(1,m2), point_xy(2,m2), 0.0D+00 )

    if ( lr /= 0 ) then
      exit
    end if

    j = j + 1

  end do
!
!  Set up the triangle information for (M1,M2,M), and for any other
!  triangles you created because points were collinear with M1, M2.
!
  tri_num = j - 2

  if ( lr == -1 ) then

    tri_vert(1,1) = m1
    tri_vert(2,1) = m2
    tri_vert(3,1) = m
    tri_nabe(3,1) = -3

    do i = 2, tri_num

      m1 = m2
      m2 = i+1
      tri_vert(1,i) = m1
      tri_vert(2,i) = m2
      tri_vert(3,i) = m
      tri_nabe(1,i-1) = -3 * i
      tri_nabe(2,i-1) = i
      tri_nabe(3,i) = i - 1

    end do

    tri_nabe(1,tri_num) = -3 * tri_num - 1
    tri_nabe(2,tri_num) = -5
    ledg = 2
    ltri = tri_num

  else

    tri_vert(1,1) = m2
    tri_vert(2,1) = m1
    tri_vert(3,1) = m
    tri_nabe(1,1) = -4

    do i = 2, tri_num
      m1 = m2
      m2 = i+1
      tri_vert(1,i) = m2
      tri_vert(2,i) = m1
      tri_vert(3,i) = m
      tri_nabe(3,i-1) = i
      tri_nabe(1,i) = -3 * i - 3
      tri_nabe(2,i) = i - 1
    end do

    tri_nabe(3,tri_num) = -3 * tri_num
    tri_nabe(2,1) = -3 * tri_num - 2
    ledg = 2
    ltri = 1

  end if
!
!  Insert the vertices one at a time from outside the convex hull,
!  determine visible boundary edges, and apply diagonal edge swaps until
!  Delaunay triangulation of vertices (so far) is obtained.
!
  top = 0

  do i = j+1, point_num

    m = i
    m1 = tri_vert(ledg,ltri)

    if ( ledg <= 2 ) then
      m2 = tri_vert(ledg+1,ltri)
    else
      m2 = tri_vert(1,ltri)
    end if

    lr = lrline ( point_xy(1,m), point_xy(2,m), point_xy(1,m1), &
      point_xy(2,m1), point_xy(1,m2), point_xy(2,m2), 0.0D+00 )

    if ( 0 < lr ) then
      rtri = ltri
      redg = ledg
      ltri = 0
    else
      l = -tri_nabe(ledg,ltri)
      rtri = l / 3
      redg = mod(l,3) + 1
    end if

    call vbedg ( point_xy(1,m), point_xy(2,m), point_num, point_xy, tri_num, &
      tri_vert, tri_nabe, ltri, ledg, rtri, redg )

    n = tri_num + 1
    l = -tri_nabe(ledg,ltri)

    do

      t = l / 3
      e = mod ( l, 3 ) + 1
      l = -tri_nabe(e,t)
      m2 = tri_vert(e,t)

      if ( e <= 2 ) then
        m1 = tri_vert(e+1,t)
      else
        m1 = tri_vert(1,t)
      end if

      tri_num = tri_num + 1
      tri_nabe(e,t) = tri_num
      tri_vert(1,tri_num) = m1
      tri_vert(2,tri_num) = m2
      tri_vert(3,tri_num) = m
      tri_nabe(1,tri_num) = t
      tri_nabe(2,tri_num) = tri_num - 1
      tri_nabe(3,tri_num) = tri_num + 1
      top = top + 1

      if ( point_num < top ) then
        ierr = 8
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DTRIS2 - Fatal error!'
        write ( *, '(a)' ) '  Stack overflow.'
        return
      end if

      stack(top) = tri_num

      if ( t == rtri .and. e == redg ) then
        exit
      end if

    end do

    tri_nabe(ledg,ltri) = -3 * n - 1
    tri_nabe(2,n) = -3 * tri_num - 2
    tri_nabe(3,tri_num) = -l
    ltri = n
    ledg = 2

    call swapec ( m, top, ltri, ledg, point_num, point_xy, tri_num, &
      tri_vert, tri_nabe, stack, ierr )

    if ( ierr /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DTRIS2 - Fatal error!'
      write ( *, '(a)' ) '  Error return from SWAPEC.'
      return
    end if

  end do
!
!  Now account for the sorting that we did.
!
  do i = 1, 3
    do j = 1, tri_num
      tri_vert(i,j) = indx ( tri_vert(i,j) )
    end do
  end do

  call perm_inverse ( point_num, indx )

  call r82vec_permute ( point_num, indx, point_xy )

  return
end
subroutine file_column_count ( input_file_name, column_num )

!*****************************************************************************80
!
!! FILE_COLUMN_COUNT counts the number of columns in the first line of a file.
!
!  Discussion:
!
!    The file is assumed to be a simple text file.
!
!    Most lines of the file is presumed to consist of COLUMN_NUM words,
!    separated by spaces.  There may also be some blank lines, and some 
!    comment lines,
!    which have a "#" in column 1.
!
!    The routine tries to find the first non-comment non-blank line and
!    counts the number of words in that line.
!
!    If all lines are blanks or comments, it goes back and tries to analyze
!    a comment line.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 June 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILE_NAME, the name of the file.
!
!    Output, integer ( kind = 4 ) COLUMN_NUM, the number of columns in 
!    the file.
!
  implicit none

  integer   ( kind = 4 )   column_num
  logical                  got_one
  character ( len = * )    input_file_name
  integer   ( kind = 4 )   input_status
  integer   ( kind = 4 )   input_unit
  character ( len = 255 )  line
!
!  Open the file.
!
  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_file_name, status = 'old', &
    form = 'formatted', access = 'sequential', iostat = input_status )

  if ( input_status /= 0 ) then
    column_num = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COLUMN_COUNT - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the input file "' &
      // trim ( input_file_name ) // '" on unit ', input_unit
    return
  end if
!
!  Read one line, but skip blank lines and comment lines.
!
  got_one = .false.

  do

    read ( input_unit, '(a)', iostat = input_status ) line

    if ( input_status /= 0 ) then
      exit
    end if

    if ( len_trim ( line ) == 0 ) then
      cycle
    end if

    if ( line(1:1) == '#' ) then
      cycle
    end if

    got_one = .true.
    exit

  end do

  if ( .not. got_one ) then

    rewind ( input_unit )

    do

      read ( input_unit, '(a)', iostat = input_status ) line

      if ( input_status /= 0 ) then
        exit
      end if

      if ( len_trim ( line ) == 0 ) then
        cycle
      end if

      got_one = .true.
      exit

    end do

  end if

  close ( unit = input_unit )

  if ( .not. got_one ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COLUMN_COUNT - Warning!'
    write ( *, '(a)' ) '  The file does not seem to contain any data.'
    column_num = -1
    return
  end if

  call s_word_count ( line, column_num )

  return
end
subroutine file_row_count ( input_file_name, row_num )

!*****************************************************************************80
!
!! FILE_ROW_COUNT counts the number of row records in a file.
!
!  Discussion:
!
!    It does not count lines that are blank, or that begin with a
!    comment symbol '#'.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILE_NAME, the name of the input file.
!
!    Output, integer ( kind = 4 ) ROW_NUM, the number of rows found.
!
  implicit none

  integer   ( kind = 4 )   bad_num
  integer   ( kind = 4 )   comment_num
  integer   ( kind = 4 )   ierror
  character ( len = * )    input_file_name
  integer   ( kind = 4 )   input_status
  integer   ( kind = 4 )   input_unit
  character ( len = 255 )  line
  integer   ( kind = 4 )   record_num
  integer   ( kind = 4 )   row_num

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_file_name, status = 'old', &
    iostat = input_status )

  if ( input_status /= 0 ) then
    row_num = -1;
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_ROW_COUNT - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the input file "' // &
      trim ( input_file_name ) // '" on unit ', input_unit
    stop
  end if

  comment_num = 0
  row_num = 0
  record_num = 0
  bad_num = 0

  do

    read ( input_unit, '(a)', iostat = input_status ) line

    if ( input_status /= 0 ) then
      ierror = record_num
      exit
    end if

    record_num = record_num + 1

    if ( line(1:1) == '#' ) then
      comment_num = comment_num + 1
      cycle
    end if

    if ( len_trim ( line ) == 0 ) then
      comment_num = comment_num + 1
      cycle
    end if

    row_num = row_num + 1

  end do

  close ( unit = input_unit )

  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is an integer between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IUNIT, the free unit number.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  logical              lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

      inquire ( unit = i, opened = lopen, iostat = ios )

      if ( ios == 0 ) then
        if ( .not. lopen ) then
          iunit = i
          return
        end if
      end if

    end if

  end do

  return
end
function i4_modp ( i, j )

!*****************************************************************************80
!
!! I4_MODP returns the nonnegative remainder of I4 division.
!
!  Discussion:
!
!    If
!      NREM = I4_MODP ( I, J )
!      NMULT = ( I - NREM ) / J
!    then
!      I = J * NMULT + NREM
!    where NREM is always nonnegative.
!
!    The MOD function computes a result with the same sign as the
!    quantity being divided.  Thus, suppose you had an angle A,
!    and you wanted to ensure that it was between 0 and 360.
!    Then mod(A,360) would do, if A was positive, but if A
!    was negative, your result would be between -360 and 0.
!
!    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Example:
!
!        I     J     MOD I4_MODP    Factorization
!
!      107    50       7       7    107 =  2 *  50 + 7
!      107   -50       7       7    107 = -2 * -50 + 7
!     -107    50      -7      43   -107 = -3 *  50 + 43
!     -107   -50      -7      43   -107 =  3 * -50 + 43
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the number to be divided.
!
!    Input, integer ( kind = 4 ) J, the number that divides I.
!
!    Output, integer ( kind = 4 ) I4_MODP, the nonnegative remainder when I is
!    divided by J.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) j
  integer ( kind = 4 ) value

  if ( j == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_MODP - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal divisor J = ', j
    stop
  end if

  value = mod ( i, j )

  if ( value < 0 ) then
    value = value + abs ( j )
  end if

  i4_modp = value

  return
end
function i4_sign ( x )

!*****************************************************************************80
!
!! I4_SIGN evaluates the sign of an I4.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X, the number whose sign is desired.
!
!    Output, integer ( kind = 4 ) I4_SIGN, the sign of X:
!
  implicit none

  integer ( kind = 4 ) i4_sign
  integer ( kind = 4 ) x

  if ( x < 0 ) then
    i4_sign = -1
  else
    i4_sign = +1
  end if

  return
end
function i4_wrap ( ival, ilo, ihi )

!*****************************************************************************80
!
!! I4_WRAP forces an I4 to lie between given limits by wrapping.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Example:
!
!    ILO = 4, IHI = 8
!
!    I  Value
!
!    -2     8
!    -1     4
!     0     5
!     1     6
!     2     7
!     3     8
!     4     4
!     5     5
!     6     6
!     7     7
!     8     8
!     9     4
!    10     5
!    11     6
!    12     7
!    13     8
!    14     4
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 August 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IVAL, a value.
!
!    Input, integer ( kind = 4 ) ILO, IHI, the desired bounds.
!
!    Output, integer ( kind = 4 ) I4_WRAP, a "wrapped" version of the value.
!
  implicit none

  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) value
  integer ( kind = 4 ) wide

  jlo = min ( ilo, ihi )
  jhi = max ( ilo, ihi )

  wide = jhi - jlo + 1

  if ( wide == 1 ) then
    value = jlo
  else
    value = jlo + i4_modp ( ival - jlo, wide )
  end if

  i4_wrap = value

  return
end
subroutine i4mat_transpose_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! I4MAT_TRANSPOSE_PRINT_SOME prints some of the transpose of an I4MAT.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of I4 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ILO, JLO, the first row and column to print.
!
!    Input, integer IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer   ( kind = 4 ), parameter :: incx = 10
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) n

  integer   ( kind = 4 ) a(m,n)
  character ( len = 8 )  ctemp(incx)
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) i2
  integer   ( kind = 4 ) i2hi
  integer   ( kind = 4 ) i2lo
  integer   ( kind = 4 ) ihi
  integer   ( kind = 4 ) ilo
  integer   ( kind = 4 ) inc
  integer   ( kind = 4 ) j
  integer   ( kind = 4 ) j2hi
  integer   ( kind = 4 ) j2lo
  integer   ( kind = 4 ) jhi
  integer   ( kind = 4 ) jlo
  character ( len = * )  title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  do i2lo = max ( ilo, 1 ), min ( ihi, m ), incx

    i2hi = i2lo + incx - 1
    i2hi = min ( i2hi, m )
    i2hi = min ( i2hi, ihi )

    inc = i2hi + 1 - i2lo

    write ( *, '(a)' ) ' '

    do i = i2lo, i2hi
      i2 = i + 1 - i2lo
      write ( ctemp(i2), '(i8)' ) i
    end do

    write ( *, '(''  Row '',10a8)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Col'
    write ( *, '(a)' ) ' '

    j2lo = max ( jlo, 1 )
    j2hi = min ( jhi, n )

    do j = j2lo, j2hi

      do i2 = 1, inc

        i = i2lo - 1 + i2

        write ( ctemp(i2), '(i8)' ) a(i,j)

      end do

      write ( *, '(i5,1x,10a8)' ) j, ( ctemp(i), i = 1, inc )

    end do

  end do

  write ( *, '(a)' ) ' '

  return
end
subroutine i4vec_indicator ( n, a )

!*****************************************************************************80
!
!! I4VEC_INDICATOR sets an I4VEC to the indicator vector.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Output, integer ( kind = 4 ) A(N), the array to be initialized.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i

  do i = 1, n
    a(i) = i
  end do

  return
end
subroutine i4vec_sort_heap_index_a ( n, a, indx )

!*****************************************************************************80
!
!! I4VEC_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    The sorting is not actually carried out.  Rather an index array is
!    created which defines the sorting.  This array may be used to sort
!    or index the array, or to sort or index related arrays keyed on the
!    original array.
!
!    Once the index array is computed, the sorting can be carried out
!    "implicitly:
!
!      A(INDX(1:N)) is sorted,
!
!    or explicitly, by the call
!
!      call i4vec_permute ( n, indx, a )
!
!    after which A(1:N) is sorted.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, integer ( kind = 4 ) A(N), an array to be index-sorted.
!
!    Output, integer ( kind = 4 ) INDX(N), the sort index.  The
!    I-th element of the sorted array is A(INDX(I)).
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) aval
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) indxt
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l

  if ( n < 1 ) then
    return
  end if

  do i = 1, n
    indx(i) = i
  end do

  if ( n == 1 ) then
    return
  end if

  l = n / 2 + 1
  ir = n

  do

    if ( 1 < l ) then

      l = l - 1
      indxt = indx(l)
      aval = a(indxt)

    else

      indxt = indx(ir)
      aval = a(indxt)
      indx(ir) = indx(1)
      ir = ir - 1

      if ( ir == 1 ) then
        indx(1) = indxt
        exit
      end if

    end if

    i = l
    j = l + l

    do while ( j <= ir )

      if ( j < ir ) then
        if ( a(indx(j)) < a(indx(j+1)) ) then
          j = j + 1
        end if
      end if

      if ( aval < a(indx(j)) ) then
        indx(i) = indx(j)
        i = j
        j = j + j
      else
        j = ir + 1
      end if

    end do

    indx(i) = indxt

  end do

  return
end
subroutine itable_write0 ( output_file_name, m, n, table )

!*****************************************************************************80
!
!! ITABLE_WRITE0 writes an ITABLE file with no headers.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OUTPUT_FILE_NAME, the output file name.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, integer ( kind = 4 ) TABLE(M,N), the table data.
!
  implicit none

  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) n

  integer   ( kind = 4 ) j
  character ( len = * )  output_file_name
  integer   ( kind = 4 ) output_status
  integer   ( kind = 4 ) output_unit
  character ( len = 30 ) string
  integer   ( kind = 4 ) table(m,n)
!
!  Open the file.
!
  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_file_name, &
    status = 'replace', iostat = output_status )

  if ( output_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ITABLE_WRITE0 - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the output file "' // &
      trim ( output_file_name ) // '" on unit ', output_unit
    output_unit = -1
    stop
  end if
!
!  Create a format string.
!
  write ( string, '(a1,i8,a4)' ) '(', m, 'i10)'
!
!  Write the data.
!
  do j = 1, n
    write ( output_unit, string ) table(1:m,j)
  end do
!
!  Close the file.
!
  close ( unit = output_unit )

  return
end
function lrline ( xu, yu, xv1, yv1, xv2, yv2, dv )

!*****************************************************************************80
!
!! LRLINE determines if a point is left of, right or, or on a directed line.
!
!  Discussion:
!
!    The directed line is parallel to, and at a signed distance DV from
!    a directed base line from (XV1,YV1) to (XV2,YV2).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 July 2001
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XU, YU, the coordinates of the point whose
!    position relative to the directed line is to be determined.
!
!    Input, real ( kind = 8 ) XV1, YV1, XV2, YV2, the coordinates of two points
!    that determine the directed base line.
!
!    Input, real ( kind = 8 ) DV, the signed distance of the directed line
!    from the directed base line through the points (XV1,YV1) and (XV2,YV2).
!    DV is positive for a line to the left of the base line.
!
!    Output, integer ( kind = 4 ) LRLINE, the result:
!    +1, the point is to the right of the directed line;
!     0, the point is on the directed line;
!    -1, the point is to the left of the directed line.
!
  implicit none

  real    ( kind = 8 ) dv
  real    ( kind = 8 ) dx
  real    ( kind = 8 ) dxu
  real    ( kind = 8 ) dy
  real    ( kind = 8 ) dyu
  integer ( kind = 4 ) lrline
  real    ( kind = 8 ) t
  real    ( kind = 8 ) tol
  real    ( kind = 8 ) tolabs
  real    ( kind = 8 ) xu
  real    ( kind = 8 ) xv1
  real    ( kind = 8 ) xv2
  real    ( kind = 8 ) yu
  real    ( kind = 8 ) yv1
  real    ( kind = 8 ) yv2

  tol = 100.0D+00 * epsilon ( tol )

  dx = xv2 - xv1
  dy = yv2 - yv1
  dxu = xu - xv1
  dyu = yu - yv1

  tolabs = tol * max ( abs ( dx ), abs ( dy ), abs ( dxu ), &
    abs ( dyu ), abs ( dv ) )

  t = dy * dxu - dx * dyu + dv * sqrt ( dx * dx + dy * dy )

  if ( tolabs < t ) then
    lrline = 1
  else if ( -tolabs <= t ) then
    lrline = 0
  else
    lrline = -1
  end if

  return
end
subroutine perm_check ( n, p, base, ierror )

!*****************************************************************************80
!
!! PERM_CHECK checks that a vector represents a permutation.
!
!  Discussion:
!
!    The routine verifies that each of the integers from BASE to
!    to BASE+N-1 occurs among the N entries of the permutation.
!
!    Set the input quantity BASE to 0, if P is a 0-based permutation,
!    or to 1 if P is a 1-based permutation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries.
!
!    Input, integer ( kind = 4 ) P(N), the array to check.
!
!    Input, integer ( kind = 4 ) BASE, the index base.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, the array represents a permutation.
!    nonzero, the array does not represent a permutation.  The smallest
!    missing value is equal to IERROR.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) base
  integer ( kind = 4 ) find
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) p(n)
  integer ( kind = 4 ) seek

  ierror = 0

  do seek = base, base + n - 1

    ierror = 1

    do find = 1, n
      if ( p(find) == seek ) then
        ierror = 0
        exit
      end if
    end do

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PERM_CHECK - Fatal error!'
      write ( *, '(a)' ) '  The input array does not represent'
      write ( *, '(a)' ) '  a proper permutation.'
      stop
    end if

  end do

  return
end
subroutine perm_inverse ( n, p )

!*****************************************************************************80
!
!! PERM_INVERSE inverts a permutation "in place".
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 January 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects being permuted.
!
!    Input/output, integer ( kind = 4 ) P(N), the permutation, in standard 
!    index form.  On output, P describes the inverse permutation
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ), parameter :: base = 1
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i0
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i4_sign
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) is
  integer ( kind = 4 ) p(n)

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PERM_INVERSE - Fatal error!'
    write ( *, '(a,i8)' ) '  Input value of N = ', n
    stop
  end if

  call perm_check ( n, p, base, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PERM_INVERSE - Fatal error!'
    write ( *, '(a)' ) '  PERM_CHECK rejects this permutation.'
    stop
  end if


  is = 1

  do i = 1, n

    i1 = p(i)

    do while ( i < i1 )
      i2 = p(i1)
      p(i1) = -i2
      i1 = i2
    end do

    is = - i4_sign ( p(i) )
    p(i) = is * abs ( p(i) )

  end do

  do i = 1, n

    i1 = - p(i)

    if ( 0 <= i1 ) then

      i0 = i

      do

        i2 = p(i1)
        p(i1) = i0

        if ( i2 < 0 ) then
          exit
        end if

        i0 = i1
        i1 = i2

      end do

    end if

  end do

  return
end
subroutine r82vec_permute ( n, p, a )

!*****************************************************************************80
!
!! R82VEC_PERMUTE permutes an R82VEC in place.
!
!  Discussion:
!
!    An R82VEC is an array of pairs of R8 values.
!
!    The same logic can be used to permute an array of objects of any 
!    arithmetic type, or an array of objects of any complexity.  The only
!    temporary storage required is enough to store a single object.  The number
!    of data movements made is N + the number of cycles of order 2 or more,
!    which is never more than N + N/2.
!
!  Example:
!
!    Input:
!
!      N = 5
!      P = (   2,    4,    5,    1,    3 )
!      A = ( 1.0,  2.0,  3.0,  4.0,  5.0 )
!          (11.0, 22.0, 33.0, 44.0, 55.0 )
!
!    Output:
!
!      A    = (  2.0,  4.0,  5.0,  1.0,  3.0 )
!             ( 22.0, 44.0, 55.0, 11.0, 33.0 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects.
!
!    Input, integer ( kind = 4 ) P(N), the permutation.  P(I) = J means
!    that the I-th element of the output array should be the J-th
!    element of the input array.  
!
!    Input/output, real ( kind = 8 ) A(2,N), the array to be permuted.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: dim_num = 2

  real    ( kind = 8 ) a(dim_num,n)
  real    ( kind = 8 ) a_temp(dim_num)
  integer ( kind = 4 ), parameter :: base = 1
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iget
  integer ( kind = 4 ) iput
  integer ( kind = 4 ) istart
  integer ( kind = 4 ) p(n)

  call perm_check ( n, p, base, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R82VEC_PERMUTE - Fatal error!'
    write ( *, '(a)' ) '  PERM_CHECK rejects this permutation.'
    stop
  end if
!
!  Search for the next element of the permutation that has not been used.
!
  do istart = 1, n

    if ( p(istart) < 0 ) then

      cycle

    else if ( p(istart) == istart ) then

      p(istart) = - p(istart)
      cycle

    else

      a_temp(1:dim_num) = a(1:dim_num,istart)
      iget = istart
!
!  Copy the new value into the vacated entry.
!
      do

        iput = iget
        iget = p(iget)

        p(iput) = - p(iput)

        if ( iget < 1 .or. n < iget ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R82VEC_PERMUTE - Fatal error!'
          write ( *, '(a)' ) '  A permutation index is out of range.'
          write ( *, '(a,i8,a,i8)' ) '  P(', iput, ') = ', iget
          stop
        end if

        if ( iget == istart ) then
          a(1:dim_num,iput) = a_temp(1:dim_num)
          exit
        end if

        a(1:dim_num,iput) = a(1:dim_num,iget)

      end do

    end if

  end do
!
!  Restore the signs of the entries.
!
  p(1:n) = - p(1:n)

  return
end
subroutine r82vec_sort_heap_index_a ( n, a, indx )

!*****************************************************************************80
!
!! R82VEC_SORT_HEAP_INDEX_A ascending index heaps an R82VEC.
!
!  Discussion:
!
!    An R82VEC is an array of R82's.
!
!    The sorting is not actually carried out.  Rather an index array is
!    created which defines the sorting.  This array may be used to sort
!    or index the array, or to sort or index related arrays keyed on the
!    original array.
!
!    Once the index array is computed, the sorting can be carried out
!    "implicitly:
!
!      A(1:2,INDX(1:N)) is sorted,
!
!    or explicitly, by the call
!
!      call r82vec_permute ( n, indx, a )
!
!    after which A(1:2,I), I = 1 to N is sorted.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = 8 ) A(2,N), an array to be index-sorted.
!
!    Output, integer ( kind = 4 ) INDX(N), the sort index.  The
!    I-th element of the sorted array is A(1:2,INDX(I)).
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: dim_num = 2

  real    ( kind = 8 ) a(dim_num,n)
  real    ( kind = 8 ) aval(dim_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) indxt
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l

  if ( n < 1 ) then
    return
  end if

  do i = 1, n
    indx(i) = i
  end do

  if ( n == 1 ) then
    return
  end if

  l = n / 2 + 1
  ir = n

  do

    if ( 1 < l ) then

      l = l - 1
      indxt = indx(l)
      aval(1:dim_num) = a(1:dim_num,indxt)

    else

      indxt = indx(ir)
      aval(1:dim_num) = a(1:dim_num,indxt)
      indx(ir) = indx(1)
      ir = ir - 1

      if ( ir == 1 ) then
        indx(1) = indxt
        exit
      end if

    end if

    i = l
    j = l + l

    do while ( j <= ir )

      if ( j < ir ) then
        if (   a(1,indx(j)) <  a(1,indx(j+1)) .or. &
             ( a(1,indx(j)) == a(1,indx(j+1)) .and. &
               a(2,indx(j)) <  a(2,indx(j+1)) ) ) then
          j = j + 1
        end if
      end if

      if (   aval(1) <  a(1,indx(j)) .or. &
           ( aval(1) == a(1,indx(j)) .and. &
             aval(2) <  a(2,indx(j)) ) ) then
        indx(i) = indx(j)
        i = j
        j = j + j
      else
        j = ir + 1
      end if

    end do

    indx(i) = indxt

  end do

  return
end
subroutine r8mat_transpose_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
!
!    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer   ( kind = 4 ), parameter :: incx = 5
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) n

  real      ( kind = 8 ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) i2
  integer   ( kind = 4 ) i2hi
  integer   ( kind = 4 ) i2lo
  integer   ( kind = 4 ) ihi
  integer   ( kind = 4 ) ilo
  integer   ( kind = 4 ) inc
  integer   ( kind = 4 ) j
  integer   ( kind = 4 ) j2hi
  integer   ( kind = 4 ) j2lo
  integer   ( kind = 4 ) jhi
  integer   ( kind = 4 ) jlo
  character ( len = * )  title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  do i2lo = max ( ilo, 1 ), min ( ihi, m ), incx

    i2hi = i2lo + incx - 1
    i2hi = min ( i2hi, m )
    i2hi = min ( i2hi, ihi )

    inc = i2hi + 1 - i2lo

    write ( *, '(a)' ) ' '

    do i = i2lo, i2hi
      i2 = i + 1 - i2lo
      write ( ctemp(i2), '(i8,6x)' ) i
    end do

    write ( *, '(''  Row   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Col'
    write ( *, '(a)' ) ' '

    j2lo = max ( jlo, 1 )
    j2hi = min ( jhi, n )

    do j = j2lo, j2hi

      do i2 = 1, inc
        i = i2lo - 1 + i2
        write ( ctemp(i2), '(g14.6)' ) a(i,j)
      end do

      write ( *, '(i5,1x,5a14)' ) j, ( ctemp(i), i = 1, inc )

    end do

  end do

  write ( *, '(a)' ) ' '

  return
end
subroutine s_blank_delete ( s )

!*****************************************************************************80
!
!! S_BLANK_DELETE removes blanks from a string, left justifying the remainder.
!
!  Discussion:
!
!    All TAB characters are also removed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string to be transformed.
!
  implicit none

  character              c
  integer   ( kind = 4 ) get
  integer   ( kind = 4 ) put
  integer   ( kind = 4 ) nchar
  character ( len = * )  s
  character, parameter :: TAB = char ( 9 )

  put = 0
  nchar = len_trim ( s )

  do get = 1, nchar

    c = s(get:get)

    if ( c /= ' ' .and. c /= TAB ) then
      put = put + 1
      s(put:put) = c
    end if

  end do

  s(put+1:nchar) = ' '

  return
end
subroutine s_to_r8 ( s, dval, ierror, length )

!*****************************************************************************80
!
!! S_TO_R8 reads an R8 from a string.
!
!  Discussion:
!
!    The routine will read as many characters as possible until it reaches
!    the end of the string, or encounters a character which cannot be
!    part of the number.
!
!    Legal input is:
!
!       1 blanks,
!       2 '+' or '-' sign,
!       2.5 blanks
!       3 integer part,
!       4 decimal point,
!       5 fraction part,
!       6 'E' or 'e' or 'D' or 'd', exponent marker,
!       7 exponent sign,
!       8 exponent integer part,
!       9 exponent decimal point,
!      10 exponent fraction part,
!      11 blanks,
!      12 final comma or semicolon,
!
!    with most quantities optional.
!
!  Example:
!
!    S                 DVAL
!
!    '1'               1.0
!    '     1   '       1.0
!    '1A'              1.0
!    '12,34,56'        12.0
!    '  34 7'          34.0
!    '-1E2ABCD'        -100.0
!    '-1X2ABCD'        -1.0
!    ' 2E-1'           0.2
!    '23.45'           23.45
!    '-4.2E+2'         -420.0
!    '17d2'            1700.0
!    '-14e-2'         -0.14
!    'e2'              100.0
!    '-12.73e-9.23'   -12.73 * 10.0**(-9.23)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string containing the
!    data to be read.  Reading will begin at position 1 and
!    terminate at the end of the string, or when no more
!    characters can be read to form a legal real.  Blanks,
!    commas, or other nonnumeric data will, in particular,
!    cause the conversion to halt.
!
!    Output, real ( kind = 8 ) DVAL, the value read from the string.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no errors occurred.
!    1, 2, 6 or 7, the input number was garbled.  The
!    value of IERROR is the last type of input successfully
!    read.  For instance, 1 means initial blanks, 2 means
!    a plus or minus sign, and so on.
!
!    Output, integer ( kind = 4 ) LENGTH, the number of characters read
!    to form the number, including any terminating
!    characters such as a trailing comma or blanks.
!
  implicit none

  character              c
  logical                ch_eqi
  real      ( kind = 8 ) dval
  integer   ( kind = 4 ) ierror
  integer   ( kind = 4 ) ihave
  integer   ( kind = 4 ) isgn
  integer   ( kind = 4 ) iterm
  integer   ( kind = 4 ) jbot
  integer   ( kind = 4 ) jsgn
  integer   ( kind = 4 ) jtop
  integer   ( kind = 4 ) length
  integer   ( kind = 4 ) nchar
  integer   ( kind = 4 ) ndig
  real      ( kind = 8 ) rbot
  real      ( kind = 8 ) rexp
  real      ( kind = 8 ) rtop
  character ( len = * )  s

  nchar = len_trim ( s )

  ierror = 0
  dval = 0.0D+00
  length = -1
  isgn = 1
  rtop = 0
  rbot = 1
  jsgn = 1
  jtop = 0
  jbot = 1
  ihave = 1
  iterm = 0

  do

    length = length + 1

    if ( nchar < length+1 ) then
      exit
    end if

    c = s(length+1:length+1)
!
!  Blank character.
!
    if ( c == ' ' ) then

      if ( ihave == 2 ) then

      else if ( ihave == 6 .or. ihave == 7 ) then
        iterm = 1
      else if ( 1 < ihave ) then
        ihave = 11
      end if
!
!  Comma.
!
    else if ( c == ',' .or. c == ';' ) then

      if ( ihave /= 1 ) then
        iterm = 1
        ihave = 12
        length = length + 1
      end if
!
!  Minus sign.
!
    else if ( c == '-' ) then

      if ( ihave == 1 ) then
        ihave = 2
        isgn = -1
      else if ( ihave == 6 ) then
        ihave = 7
        jsgn = -1
      else
        iterm = 1
      end if
!
!  Plus sign.
!
    else if ( c == '+' ) then

      if ( ihave == 1 ) then
        ihave = 2
      else if ( ihave == 6 ) then
        ihave = 7
      else
        iterm = 1
      end if
!
!  Decimal point.
!
    else if ( c == '.' ) then

      if ( ihave < 4 ) then
        ihave = 4
      else if ( 6 <= ihave .and. ihave <= 8 ) then
        ihave = 9
      else
        iterm = 1
      end if
!
!  Scientific notation exponent marker.
!
    else if ( ch_eqi ( c, 'E' ) .or. ch_eqi ( c, 'D' ) ) then

      if ( ihave < 6 ) then
        ihave = 6
      else
        iterm = 1
      end if
!
!  Digit.
!
    else if (  ihave < 11 .and. lle ( '0', c ) .and. lle ( c, '9' ) ) then

      if ( ihave <= 2 ) then
        ihave = 3
      else if ( ihave == 4 ) then
        ihave = 5
      else if ( ihave == 6 .or. ihave == 7 ) then
        ihave = 8
      else if ( ihave == 9 ) then
        ihave = 10
      end if

      call ch_to_digit ( c, ndig )

      if ( ihave == 3 ) then
        rtop = 10.0D+00 * rtop + real ( ndig, kind = 8 )
      else if ( ihave == 5 ) then
        rtop = 10.0D+00 * rtop + real ( ndig, kind = 8 )
        rbot = 10.0D+00 * rbot
      else if ( ihave == 8 ) then
        jtop = 10 * jtop + ndig
      else if ( ihave == 10 ) then
        jtop = 10 * jtop + ndig
        jbot = 10 * jbot
      end if
!
!  Anything else is regarded as a terminator.
!
    else
      iterm = 1
    end if
!
!  If we haven't seen a terminator, and we haven't examined the
!  entire string, go get the next character.
!
    if ( iterm == 1 ) then
      exit
    end if

  end do
!
!  If we haven't seen a terminator, and we have examined the
!  entire string, then we're done, and LENGTH is equal to NCHAR.
!
  if ( iterm /= 1 .and. length+1 == nchar ) then
    length = nchar
  end if
!
!  Number seems to have terminated.  Have we got a legal number?
!  Not if we terminated in states 1, 2, 6 or 7!
!
  if ( ihave == 1 .or. ihave == 2 .or. ihave == 6 .or. ihave == 7 ) then
    ierror = ihave
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'S_TO_R8 - Serious error!'
    write ( *, '(a)' ) '  Illegal or nonnumeric input:'
    write ( *, '(a)' ) '    ' // trim ( s )
    return
  end if
!
!  Number seems OK.  Form it.
!
  if ( jtop == 0 ) then
    rexp = 1.0D+00
  else
    if ( jbot == 1 ) then
      rexp = 10.0D+00 ** ( jsgn * jtop )
    else
      rexp = 10.0D+00 ** ( real ( jsgn * jtop, kind = 8 ) &
        / real ( jbot, kind = 8 ) )
    end if
  end if

  dval = real ( isgn, kind = 8 ) * rexp * rtop / rbot

  return
end
subroutine s_to_r8vec ( s, n, rvec, ierror )

!*****************************************************************************80
!
!! S_TO_R8VEC reads an R8VEC from a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string to be read.
!
!    Input, integer ( kind = 4 ) N, the number of values expected.
!
!    Output, real ( kind = 8 ) RVEC(N), the values read from the string.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no errors occurred.
!    -K, could not read data for entries -K through N.
!
  implicit none

  integer   ( kind = 4 ) n

  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) ierror
  integer   ( kind = 4 ) ilo
  integer   ( kind = 4 ) lchar
  real      ( kind = 8 ) rvec(n)
  character ( len = * )  s

  i = 0
  ierror = 0
  ilo = 1

  do while ( i < n )

    i = i + 1

    call s_to_r8 ( s(ilo:), rvec(i), ierror, lchar )

    if ( ierror /= 0 ) then
      ierror = -i
      exit
    end if

    ilo = ilo + lchar

  end do

  return
end
subroutine s_word_count ( s, nword )

!*****************************************************************************80
!
!! S_WORD_COUNT counts the number of "words" in a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string to be examined.
!
!    Output, integer ( kind = 4 ) NWORD, the number of "words" in the string.
!    Words are presumed to be separated by one or more blanks.
!
  implicit none

  logical                blank
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) lens
  integer   ( kind = 4 ) nword
  character ( len = * )  s

  nword = 0
  lens = len ( s )

  if ( lens <= 0 ) then
    return
  end if

  blank = .true.

  do i = 1, lens

    if ( s(i:i) == ' ' ) then
      blank = .true.
    else if ( blank ) then
      nword = nword + 1
      blank = .false.
    end if

  end do

  return
end
subroutine swapec ( i, top, btri, bedg, point_num, point_xy, tri_num, &
  tri_vert, tri_nabe, stack, ierr )

!*****************************************************************************80
!
!! SWAPEC swaps diagonal edges until all triangles are Delaunay.
!
!  Discussion:
!
!    The routine swaps diagonal edges in a 2D triangulation, based on
!    the empty circumcircle criterion, until all triangles are Delaunay,
!    given that I is the index of the new vertex added to the triangulation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 July 2001
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms, 
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the index of the new vertex.
!
!    Input/output, integer ( kind = 4 ) TOP, the index of the top of the stack.
!    On output, TOP is zero.
!
!    Input/output, integer ( kind = 4 ) BTRI, BEDG; on input, if positive, are 
!    the triangle and edge indices of a boundary edge whose updated indices
!    must be recorded.  On output, these may be updated because of swaps.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, real ( kind = 8 ) POINT_XY(2,POINT_NUM), the coordinates
!    of the points.
!
!    Input, integer ( kind = 4 ) TRI_NUM, the number of triangles.
!
!    Input/output, integer ( kind = 4 ) TRI_VERT(3,TRI_NUM), the triangle 
!    incidence list.  May be updated on output because of swaps.
!
!    Input/output, integer ( kind = 4 ) TRI_NABE(3,TRI_NUM), the triangle 
!    neighbor list; negative values are used for links of the counter-clockwise 
!    linked list of boundary edges;  May be updated on output because of swaps.
!      LINK = -(3*I + J-1) where I, J = triangle, edge index.
!
!    Workspace, integer ( kind = 4 ) STACK(MAXST); on input, entries 1 through
!    TOP contain the indices of initial triangles (involving vertex I)
!    put in stack; the edges opposite I should be in interior;  entries
!    TOP+1 through MAXST are used as a stack.
!
!    Output, integer ( kind = 4 ) IERR is set to 8 for abnormal return.
!
  implicit none

  integer ( kind = 4 ) point_num
  integer ( kind = 4 ) tri_num

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) bedg
  integer ( kind = 4 ) btri
  integer ( kind = 4 ) c
  integer ( kind = 4 ) diaedg
  integer ( kind = 4 ) e
  integer ( kind = 4 ) ee
  integer ( kind = 4 ) em1
  integer ( kind = 4 ) ep1
  integer ( kind = 4 ) f
  integer ( kind = 4 ) fm1
  integer ( kind = 4 ) fp1
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) l
  integer ( kind = 4 ) r
  integer ( kind = 4 ) s
  integer ( kind = 4 ) stack(point_num)
  integer ( kind = 4 ) swap
  integer ( kind = 4 ) t
  integer ( kind = 4 ) top
  integer ( kind = 4 ) tri_nabe(3,tri_num)
  integer ( kind = 4 ) tri_vert(3,tri_num)
  integer ( kind = 4 ) tt
  integer ( kind = 4 ) u
  real    ( kind = 8 ) point_xy(2,point_num)
  real    ( kind = 8 ) x
  real    ( kind = 8 ) y
!
!  Determine whether triangles in stack are Delaunay, and swap
!  diagonal edge of convex quadrilateral if not.
!
  x = point_xy(1,i)
  y = point_xy(2,i)

  do

    if ( top <= 0 ) then
      exit
    end if

    t = stack(top)
    top = top - 1

    if ( tri_vert(1,t) == i ) then
      e = 2
      b = tri_vert(3,t)
    else if ( tri_vert(2,t) == i ) then
      e = 3
      b = tri_vert(1,t)
    else
      e = 1
      b = tri_vert(2,t)
    end if

    a = tri_vert(e,t)
    u = tri_nabe(e,t)

    if ( tri_nabe(1,u) == t ) then
      f = 1
      c = tri_vert(3,u)
    else if ( tri_nabe(2,u) == t ) then
      f = 2
      c = tri_vert(1,u)
    else
      f = 3
      c = tri_vert(2,u)
    end if

    swap = diaedg ( x, y, point_xy(1,a), point_xy(2,a), point_xy(1,c), &
      point_xy(2,c), point_xy(1,b), point_xy(2,b) )

    if ( swap == 1 ) then

      em1 = i4_wrap ( e - 1, 1, 3 )
      ep1 = i4_wrap ( e + 1, 1, 3 )
      fm1 = i4_wrap ( f - 1, 1, 3 )
      fp1 = i4_wrap ( f + 1, 1, 3 )

      tri_vert(ep1,t) = c
      tri_vert(fp1,u) = i
      r = tri_nabe(ep1,t)
      s = tri_nabe(fp1,u)
      tri_nabe(ep1,t) = u
      tri_nabe(fp1,u) = t
      tri_nabe(e,t) = s
      tri_nabe(f,u) = r

      if ( 0 < tri_nabe(fm1,u) ) then
        top = top + 1
        stack(top) = u
      end if

      if ( 0 < s ) then

        if ( tri_nabe(1,s) == u ) then
          tri_nabe(1,s) = t
        else if ( tri_nabe(2,s) == u ) then
          tri_nabe(2,s) = t
        else
          tri_nabe(3,s) = t
        end if

        top = top + 1

        if ( point_num < top ) then
          ierr = 8
          return
        end if

        stack(top) = t

      else

        if ( u == btri .and. fp1 == bedg ) then
          btri = t
          bedg = e
        end if

        l = - ( 3 * t + e - 1 )
        tt = t
        ee = em1

        do while ( 0 < tri_nabe(ee,tt) )

          tt = tri_nabe(ee,tt)

          if ( tri_vert(1,tt) == a ) then
            ee = 3
          else if ( tri_vert(2,tt) == a ) then
            ee = 1
          else
            ee = 2
          end if

        end do

        tri_nabe(ee,tt) = l

      end if

      if ( 0 < r ) then

        if ( tri_nabe(1,r) == t ) then
          tri_nabe(1,r) = u
        else if ( tri_nabe(2,r) == t ) then
          tri_nabe(2,r) = u
        else
          tri_nabe(3,r) = u
        end if

      else

        if ( t == btri .and. ep1 == bedg ) then
          btri = u
          bedg = f
        end if

        l = - ( 3 * u + f - 1 )
        tt = u
        ee = fm1

        do while ( 0 < tri_nabe(ee,tt) )

          tt = tri_nabe(ee,tt)

          if ( tri_vert(1,tt) == b ) then
            ee = 3
          else if ( tri_vert(2,tt) == b ) then
            ee = 1
          else
            ee = 2
          end if

        end do

        tri_nabe(ee,tt) = l

      end if

    end if

  end do

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    May 31 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 )  ampm
  integer   ( kind = 4 ) d
  character ( len = 8 )  date
  integer   ( kind = 4 ) h
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer   ( kind = 4 ) n
  integer   ( kind = 4 ) s
  character ( len = 10 ) time
  integer   ( kind = 4 ) values(8)
  integer   ( kind = 4 ) y
  character ( len = 5 )  zone

  call date_and_time ( date, time, zone, values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
subroutine vbedg ( x, y, point_num, point_xy, tri_num, tri_vert, tri_nabe, &
  ltri, ledg, rtri, redg )

!*****************************************************************************80
!
!! VBEDG determines which boundary edges are visible to a point.
!
!  Discussion:
!
!    The point (X,Y) is assumed to be outside the convex hull of the
!    region covered by the 2D triangulation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 August 2001
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms, 
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, Y, the coordinates of a point outside
!    the convex hull of the current triangulation.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, real ( kind = 8 ) POINT_XY(2,POINT_NUM), the coordinates 
!    of the vertices.
!
!    Input, integer ( kind = 4 ) TRI_NUM, the number of triangles.
!
!    Input, integer ( kind = 4 ) TRI_VERT(3,TRI_NUM), the triangle incidence 
!    list.
!
!    Input, integer ( kind = 4 ) TRI_NABE(3,TRI_NUM), the triangle neighbor 
!    list; negative values are used for links of a counter clockwise linked 
!    list of boundary edges;
!      LINK = -(3*I + J-1) where I, J = triangle, edge index.
!
!    Input/output, integer ( kind = 4 ) LTRI, LEDG.  If LTRI /= 0 then these 
!    values are assumed to be already computed and are not changed, else they 
!    are updated.  On output, LTRI is the index of boundary triangle to the 
!    left of the leftmost boundary triangle visible from (X,Y), and LEDG is 
!    the boundary edge of triangle LTRI to the left of the leftmost boundary
!    edge visible from (X,Y).  1 <= LEDG <= 3.
!
!    Input/output, integer ( kind = 4 ) RTRI.  On input, the index of the 
!    boundary triangle to begin the search at.  On output, the index of the 
!    rightmost boundary triangle visible from (X,Y).
!
!    Input/output, integer ( kind = 4 ) REDG, the edge of triangle RTRI that 
!    is visible from (X,Y).  1 <= REDG <= 3.
!
  implicit none

  integer ( kind = 4 ) point_num
  integer ( kind = 4 ) tri_num

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) e
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) l
  logical              ldone
  integer ( kind = 4 ) ledg
  integer ( kind = 4 ) lr
  integer ( kind = 4 ) lrline
  integer ( kind = 4 ) ltri
  real    ( kind = 8 ) point_xy(2,point_num)
  integer ( kind = 4 ) redg
  integer ( kind = 4 ) rtri
  integer ( kind = 4 ) t
  integer ( kind = 4 ) tri_nabe(3,tri_num)
  integer ( kind = 4 ) tri_vert(3,tri_num)
  real    ( kind = 8 ) x
  real    ( kind = 8 ) y
!
!  Find the rightmost visible boundary edge using links, then possibly
!  leftmost visible boundary edge using triangle neighbor information.
!
  if ( ltri == 0 ) then
    ldone = .false.
    ltri = rtri
    ledg = redg
  else
    ldone = .true.
  end if

  do

    l = -tri_nabe(redg,rtri)
    t = l / 3
    e = mod ( l, 3 ) + 1
    a = tri_vert(e,t)

    if ( e <= 2 ) then
      b = tri_vert(e+1,t)
    else
      b = tri_vert(1,t)
    end if

    lr = lrline ( x, y, point_xy(1,a), point_xy(2,a), point_xy(1,b), &
      point_xy(2,b), 0.0D+00 )

    if ( lr <= 0 ) then
      exit
    end if

    rtri = t
    redg = e

  end do

  if ( ldone ) then
    return
  end if

  t = ltri
  e = ledg

  do

    b = tri_vert(e,t)
    e = i4_wrap ( e-1, 1, 3 )

    do while ( 0 < tri_nabe(e,t) )

      t = tri_nabe(e,t)

      if ( tri_vert(1,t) == b ) then
        e = 3
      else if ( tri_vert(2,t) == b ) then
        e = 1
      else
        e = 2
      end if

    end do

    a = tri_vert(e,t)

    lr = lrline ( x, y, point_xy(1,a), point_xy(2,a), point_xy(1,b), &
       point_xy(2,b), 0.0D+00 )

    if ( lr <= 0 ) then
      exit
    end if

  end do

  ltri = t
  ledg = e

  return
end subroutine vbedg

