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
