program event_chain
  use phys_parameters
  use cell_parameters
  use alg_parameters
  use others
  use benchmark
  implicit none
  real*8 :: L_cur,L_max,L_min     !total displacement still to complete, maximum displacement, minimum collision lenght
  real*8 :: dx,dy,g          !cpu-time, dispacement in x and y
  logical :: ok
  integer*8 :: i,j,k,k_cur,k_min,cell_cur,cell_min,disk,conf
  call resetseed() ! comment this line for different random number for each run
  call init_benchmark
  L_max=min(minval(Box)/2,minval(l_cell))-2*r
  if (L_max<0) then
     write(*,*) "L_max<0, too much cells"
     goto 9
  endif
  T=0
  t_0=0
  call init_cell
  write(*,'("Max disks/cell",i8)') n_cell_max
  write(*,'("N_b=",i8,"  Eta=",d10.4,"  L=",d8.2)') N_b,eta,L
  write(*,'("r=",d10.4"  Box=",2d11.4,"  N_cell=",2i4,"  L_max/r=",d10.4)') r,Box,N_cell,L_max/r
  write(*,*)
  i=0
  ok=.true.
  conf=0
  k_cur=-1
  cell_cur=-1
  do
     i=i+1 !new move
     call pick_move(cell_cur,k_cur)      !initial disk
     if (T>N_dpl_tot .and. N_dpl_tot/=-1) exit !condition to end the program
     L_cur=L
     ok=(.not. ok)
     if (do_benchmark) ok=.true.
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
              if (do_benchmark) call precision_check(L_cur,L_max)
              L_cur=L_cur-L_max
              j=j-1
              cycle !Lmax<L_min the move continues with the same disk and the same direction
           elseif (L_min<L_cur) then
              if (do_benchmark) call precision_check(L_cur,L_min)
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
              if (do_benchmark) call precision_check(L_cur,L_max)
              L_cur=L_cur-L_max
              j=j-1
              cycle
           elseif (L_min<L_cur) then
              if (do_benchmark) call precision_check(L_cur,L_min)
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
  write(*,'("Number of chains=",i14)') i-1
  if (create_init) then
    call write_init
    write(*,*) 'Writing new initial configuration'
  endif
  if (do_benchmark) then
    write(*,'("Final cell=",i7,"  final active particle=(",d10.4,", ",d10.4,")")') cell_cur, &
          Cell(1,k_cur,cell_cur), Cell(2,k_cur,cell_cur)
    write(*,'("Global position = (",d12.6,", ",d12.6,")")') Cell(1,k_cur,cell_cur)-box(1)/2+ &
          (modulo(cell_cur-1,N_cell(1))+0.5)*l_cell(1),&
          Cell(2,k_cur,cell_cur)-box(2)/2+((cell_cur-1)/N_cell(1)+0.5)*l_cell(2)
    write(*,*) 'Writing number of particles in each cell.'
    call write_hist_cells
  endif
9 continue
end program event_chain
