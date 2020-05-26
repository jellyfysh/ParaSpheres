!JeLLyFysh/ParaSpheres - Multithreaded event-chain Monte Carlo with local times
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
!     Event_chain_SEC-XY  Version 28/04/2010
!   - Monodisperse hard disks in the (N,V,T) ensemble (faster than (N,P,T))
!   - Disks go straight in (ox) then (oy) then (ox)... until the total displacement
!     L is reached.
!   - ~20 times faster than the local Metropolis in # of displaced disk (x100)
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
!
!   - Extract position and backup each N_extr
! 
!*************************************************************************!

module phys_parameters
  implicit none
  real*4, parameter :: pi=acos(-1.D0)                !pi
  integer*8, parameter :: N_b=256**2             !number of disk
  real*4, parameter :: eta=0.708D0              !density (volume fraction)
  real*4, parameter,dimension(1:2) :: Box=(/1.D0,1.D0/)    !Box's size
  real*4, parameter :: r=sqrt(eta*Box(1)*Box(2)/(N_b*pi))    !radius
  real*4, dimension(1:2,1:N_b) :: X
end module phys_parameters

module cell_parameters
  use phys_parameters
  integer*8, parameter,dimension(1:2) :: N_cell=(/1,1/)*nint(sqrt(real(N_b))*7.D0/8.D0)    !Cells number (1024 ~ 28*28 etc)
  real*4, parameter,dimension(1:2) :: l_cell=Box/N_cell        !cells size
  integer*8, parameter :: n_cell_max=floor((l_cell(1)/r+2.D0)*(l_cell(2)/r+2.D0)/(2.D0*sqrt(3.D0)))!>max # of disks/cell
  real*4, dimension(1:2,1:n_cell_max,1:N_cell(1)*N_cell(2)) :: Cell  !positions of the disks/ center of their cell
  integer*2, dimension(1:N_cell(1)*N_cell(2)) :: Cell_ocp=0          !# of disks for each cell
  integer*4,dimension(9,N_cell(1)*N_cell(2)) :: Cell_neighb          !Cell neighbor's table
  integer*8,dimension(1:2) :: i_cell
end module cell_parameters

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

module others
  real*4 :: time,t_0
  integer*8 :: T
end module others

program event_chain
  use phys_parameters
  use cell_parameters
  use alg_parameters
  use others
  implicit none 
  real*4 :: L_cur,L_max,L_min     !total displacement still to complete, maximum displacement, minimum collision lenght
  real*4 :: dx,dy,g          !cpu-time, dispacement in x and y
  logical :: ok
  integer*8 :: i,j,k,k_cur,k_min,cell_cur,cell_min,disk,conf
  L_max=min(minval(Box)/2,minval(l_cell))-2*r
  if (L_max<0) then 
     write(*,*) "L_max<0, too much cells" 
     goto 9
  endif
  T=0
  t_0=0
  call init_cell
  write(*,'("N_b=",i8,"  Eta=",d10.4,"  L=",d8.2)') N_b,eta,L
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
     if (mod(i,N_extr)==0 .and. N_extr/=-1) then
        call cpu_time(time)
        time=time+t_0
        write(*,'("T=",i14,"    CPU-time=",f12.2)') T,time
        call write_pos
     endif
     !***end extraction ***!
  enddo
9 continue
end program event_chain

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
  implicit none
  real*4 :: dx
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
  if (dy>l_cell(2)/2) then
     dy=dy-l_cell(2)
     cell_new=Cell_neighb(8,cell_cur)
     Cell_ocp(cell_new)=Cell_ocp(cell_new)+1 
     Cell(1,Cell_ocp(cell_new),cell_new)=Cell(1,k_cur,cell_cur)
     Cell(2,Cell_ocp(cell_new),cell_new)=dy
     Cell(:,k_cur,cell_cur)=Cell(:,Cell_ocp(cell_cur),cell_cur)
     Cell_ocp(cell_cur)=Cell_ocp(cell_cur)-1
     cell_cur=cell_new
     k_cur=Cell_ocp(cell_new)
  else
     Cell(2,k_cur,cell_cur)=dy
  endif
end subroutine refresh_cell_y

