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
