subroutine refresh_cell_y(dy,cell_cur,k_cur)
!******************************************************************!
!
! Compute the new cell of k_cur
!
!******************************************************************!
  use phys_parameters
  use cell_parameters
  implicit none
  real*8 :: dy
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
