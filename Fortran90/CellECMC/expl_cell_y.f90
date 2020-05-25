subroutine expl_cell_y(cell_cur,k_cur,L_min,k_min,cell_min)
!******************************************************************!
!
! Explore the disks neighbor of k_cur, 6 cell-neighb to explore (4,5,6,7,8,9)
! Faster without the loop
!
!******************************************************************!
  use phys_parameters
  use cell_parameters
  use benchmark
  implicit none
  real*8, dimension(1:2)::dx,x_cur
  integer*8 :: cell_cur,cell_act,cell_min
  real*8 :: L_min,L_coll
  integer*8 :: i,j,k_cur,k_min
  x_cur=Cell(:,k_cur,cell_cur)
  cell_act=Cell_neighb(8,cell_cur)
  do j=1,Cell_ocp(cell_act)
     dx(1)=Cell(1,j,cell_act)-x_cur(1)
     if (do_benchmark) call precision_check(2*r,abs(dx(1)))
     if (abs(dx(1))<2*r) then
        dx(2)=Cell(2,j,cell_act)-x_cur(2)+l_cell(2)
        L_coll=dx(2)-sqrt(4*r**2-dx(1)**2)
        if (L_coll<L_min) then
           if (do_benchmark) call precision_check(L_min,L_coll)
           L_min=L_coll
           k_min=j
           cell_min=cell_act
        endif
     endif
  enddo
  cell_act=Cell_neighb(6,cell_cur)
  do j=1,Cell_ocp(cell_act)
     dx(1)=Cell(1,j,cell_act)-x_cur(1)+l_cell(1)
     if (do_benchmark) call precision_check(2*r,abs(dx(1)))
     if (abs(dx(1))<2*r) then
        dx(2)=Cell(2,j,cell_act)-x_cur(2)
        if (dx(2)>0.D0) then
           L_coll=dx(2)-sqrt(4*r**2-dx(1)**2)
           if (L_coll<L_min) then
              if (do_benchmark) call precision_check(L_min,L_coll)
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
     if (do_benchmark) call precision_check(2*r,abs(dx(1)))
     if (abs(dx(1))<2*r) then
        dx(2)=Cell(2,j,cell_act)-x_cur(2)
        if (dx(2)>0.D0) then
           L_coll=dx(2)-sqrt(4*r**2-dx(1)**2)
           if (L_coll<L_min) then
              if (do_benchmark) call precision_check(L_min,L_coll)
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
     if (do_benchmark) call precision_check(2*r,abs(dx(1)))
     if (abs(dx(1))<2*r) then
        dx(2)=Cell(2,j,cell_act)-x_cur(2)
        if (dx(2)>0.D0) then
           L_coll=dx(2)-sqrt(4*r**2-dx(1)**2)
           if (L_coll<L_min) then
              if (do_benchmark) call precision_check(L_min,L_coll)
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
     if (do_benchmark) call precision_check(2*r,abs(dx(1)))
     if (abs(dx(1))<2*r) then
        dx(2)=Cell(2,j,cell_act)-x_cur(2)+l_cell(2)
        L_coll=dx(2)-sqrt(4*r**2-dx(1)**2)
        if (L_coll<L_min) then
           if (do_benchmark) call precision_check(L_min,L_coll)
           L_min=L_coll
           k_min=j
           cell_min=cell_act
        endif
     endif
  enddo
  cell_act=Cell_neighb(7,cell_cur)
  do j=1,Cell_ocp(cell_act)
     dx(1)=Cell(1,j,cell_act)-x_cur(1)-l_cell(1)
     if (do_benchmark) call precision_check(2*r,abs(dx(1)))
     if (abs(dx(1))<2*r) then
        dx(2)=Cell(2,j,cell_act)-x_cur(2)+l_cell(2)
        L_coll=dx(2)-sqrt(4*r**2-dx(1)**2)
        if (L_coll<L_min) then
           if (do_benchmark) call precision_check(L_min,L_coll)
           L_min=L_coll
           k_min=j
           cell_min=cell_act
        endif
     endif
  enddo
end subroutine expl_cell_y
