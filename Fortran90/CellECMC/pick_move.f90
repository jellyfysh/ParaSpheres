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
