module test_MCMC
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use test_functions
  use random_uniform
  use MHMCMC
  use OMP_LIB
  implicit none
  private

  public:: collect_MCMCtests

contains

!> Collect all exported unit tests
subroutine collect_MCMCtests(testsuite)
  implicit none
  !> Collection of tests
  type(unittest_type), allocatable, intent(out):: testsuite(:)

  testsuite = [ &
    new_unittest("step_pars", test_step_pars), &
    new_unittest("step_pars_real", test_step_pars_real), &
    new_unittest("mcmc_output_type", test_mcmc_output_type), &
    new_unittest("mcmc_len0", test_run_mcmc_len0), &
    new_unittest("mcmc_len1000", test_run_mcmc_len1000), &
    new_unittest("mcmc_nchains1_len0", test_run_parallel_mcmc_nchains1_len0), &
    new_unittest("mcmc_nchains4_len1000", test_run_parallel_mcmc_nchains4_len0), &
    new_unittest("mcmc_nchains4omp_len1000", test_run_parallel_mcmc_nchains4_enforceomp_len0) &
    ]

end subroutine collect_MCMCtests


subroutine test_step_pars(error)
  implicit none
  type(error_type), allocatable, intent(out):: error
  ! generate new proposal PARS given current point PARS0 and other specs
  ! example point x, y = 0, 0 with independent variances
  integer, parameter:: npars = 2
  double precision, dimension(npars):: pars0, pars
  logical:: multivariate = .true.
  double precision, dimension(npars, npars):: covariance
  double precision:: beta = 0.05d0  ! recommended constant, only relevant in multivariate phase
  double precision:: opt_scaling = 5.67/dble(npars)  ! recommended constant 
  double precision:: par_minstepsize = 0.001d0 ! ?? what does this do ?
  type(UNIF_VECTOR):: random_uniform
  call random_uniform%initialize()
  call init_PI()
  pars0 = (/0d0, 0d0/)
  covariance = reshape(source = [1d0, 0d0, 0d0, 1d0], shape = [2, 2])  ! uncorellated with large variance
  call step_pars(pars0, pars, npars, multivariate, covariance, beta, opt_scaling, par_minstepsize, random_uniform)
  ! expect pars is close to, but not equal pars0
  call check(error, (pars(1)/=0d0) .and. (pars(2)/=0d0))
  ! Note step_pars does not respect bouns 0..1 of normed space !!  as per originnal cardamom, checked and rejected later
  ! in bounds
  ! call check(error, (pars(1)>=0d0) .and. (pars(1)<=1d0))
  ! call check(error, (pars(2)>=0d0) .and. (pars(2)<=1d0))
end subroutine test_step_pars

subroutine test_step_pars_real(error)
  implicit none
  type(error_type), allocatable, intent(out):: error
  ! generate new proposal PARS given current point PARS0 and other specs
  ! example point x, y = 0, 0 with independent variances
  double precision, dimension(PI_xy%n_pars):: pars0, pars
  logical:: multivariate = .true.
  double precision, dimension(PI_xy%n_pars, PI_xy%n_pars):: covariance
  double precision:: beta = 0.05d0  ! recommended constant, only relevant in multivariate phase
  double precision:: opt_scaling  ! recommended constant 
  double precision:: par_minstepsize = 0.001d0 ! ?? what does this do ?
  type(UNIF_VECTOR):: random_uniform
  opt_scaling = 5.67/dble(PI_xy%n_pars)
  call random_uniform%initialize()
  call init_PI()
  pars0 = (/0d0, 3.0d0/)  ! within bounds of PI_xy
  covariance = reshape(source = [1d0, 0d0, 0d0, 1d0], shape = [2, 2])  ! uncorellated with large variance
  call step_pars_real(pars0, pars, PI_xy, multivariate, covariance, beta, opt_scaling, par_minstepsize, random_uniform)
  ! expect pars is close to, but not equal pars0
  call check(error, (pars(1)/=pars0(1)) .and. (pars(2)/=pars0(2)))
  ! Note step_pars does not respect bouns 0..1 of normed space !!  as per originnal cardamom, checked and rejected later
  ! in bounds
  !call check(error, (pars(1)>=PI_xy%parmin(1)) .and. (pars(1)<=PI_xy%parmax(1)))
  !call check(error, (pars(2)>=PI_xy%parmin(2)) .and. (pars(2)<=PI_xy%parmax(2)))
end subroutine test_step_pars_real

subroutine test_adapt_step_size(error)
  implicit none
  type(error_type), allocatable, intent(out):: error
  ! ...
  ! expect either a covariance matrix exists or pi%use_multivariate is set to false
  ! TODO should be called update_covariance
end subroutine test_adapt_step_size

subroutine test_increment_covariance_matrix(error)
  implicit none
  type(error_type), allocatable, intent(out):: error
  ! call increment_covariance_matrix(history, mean, N_local, N_history, cov)
  ! expect changes to mean, cov
  ! and n_local +=1 
end subroutine

subroutine test_cholesky_factor(error)
  implicit none
  type(error_type), allocatable, intent(out):: error
  !call cholesky_factor(N, cov, posdef)
end subroutine

subroutine test_covariance_matrix(error)
  implicit none
  type(error_type), allocatable, intent(out):: error
  !call covariance_matrix(history, mean, npars, N_history, cov)
  ! expect output at mean, cov
end subroutine

subroutine test_mcmc_output_type(error)
  use MHMCMC, only: MCMC_OUTPUT
  implicit none
  type(error_type), allocatable, intent(out):: error
  ! test declare an object of type MCMC_OUTPUT
  type(MCMC_OUTPUT):: MCOUT
end subroutine

! TODO code to run before tests, such as init_pi

subroutine test_run_mcmc_len0(error)
  implicit none
  type(error_type), allocatable, intent(out):: error
  type(mcmc_output):: mcout
  type(mcmc_options):: mcopt  ! filled with defaults only
  ! zero length run : takes expected input arguments, setup works, 
  ! outputs/writes unchanged state
  ! all on defaults, without optional arguments
  call init_pi()
  mcopt%maxiter = 0
  call run_mcmc(ll_normal, pi_xy, mcopt, mcout)
  ! expect values in mcout : a random initial state (within given parameter bounds) 
  ! and its loglikelihood
  ! call check(error, mcout%ll > 0d0 , .true. )
  write(*,*) mcout%ll, mcout%pars
  call check(error, mcout%pars(1) >= pi_xy%parmin(1) .and. mcout%pars(1) <= pi_xy%parmax(1))
  call check(error, mcout%pars(2) >= pi_xy%parmin(2) .and. mcout%pars(2) <= pi_xy%parmax(2) )
  ! given starting values which happen to be the correct (min loglik) solution, 
  ! expect after zero run steps mcout has these initial values and 
  ! lolglik is low.
  !call check(error, mcout%pars(1), x_ideal )
  !call check(error, mcout%pars(2), y_ideal)
  !call check(error, MCOUT% ll, 0.0 )
end subroutine 

subroutine test_run_mcmc_len1000(error)
  implicit none
  type(error_type), allocatable, intent(out):: error
  type(mcmc_output):: mcout
  type(mcmc_options):: mcopt  ! filled with defaults only
  ! zero length run : takes expected input arguments, setup works, 
  ! outputs/writes unchanged state
  ! all on defaults, without optional arguments
  call init_pi()
  mcopt%maxiter = 1000
  call run_mcmc(ll_normal, pi_xy, mcopt, mcout)
  ! expect values in mcout : loglikelihood and parameters in bounds
  ! call check(error, mcout%ll > 0d0 )
  call check(error, mcout%pars(1) >= pi_xy%parmin(1) .and. mcout%pars(1) <= pi_xy%parmax(1))
  call check(error, mcout%pars(2) >= pi_xy%parmin(2) .and. mcout%pars(2) <= pi_xy%parmax(2) )
  ! this is an easy problem (quadratic potential), expect convergence on solution in a short run
  !write(*,*) mcout%ll, mcout%pars
  !call check(error, mcout%pars(1), x_ideal )
  !call check(error, mcout%pars(2), y_ideal)
  !call check(error, MCOUT%ll, 50.0 )
end subroutine 

subroutine test_run_parallel_mcmc_nchains1_len0(error)
  implicit none
  integer, parameter :: nchains =1
  type(error_type), allocatable, intent(out):: error
  type(mcmc_output), dimension(:), allocatable :: mcout
  type(mcmc_output) :: mcout1
  type(mcmc_options):: mcopt  ! filled with defaults only
  ! zero length run : takes expected input arguments, setup works, 
  ! outputs/writes unchanged state
  ! all on defaults, without optional arguments
  call init_pi()
  mcopt%maxiter = 0
  call run_parallel_mcmc(ll_normal, pi_xy, mcopt, mcout)
  ! expect values in mcout : a random initial state (within given parameter bounds) 
  ! and its loglikelihood
  mcout1 = mcout(1)
  call check(error, mcout1%pars(1) >= pi_xy%parmin(1) .and. mcout1%pars(1) <= pi_xy%parmax(1))
  call check(error, mcout1%pars(2) >= pi_xy%parmin(2) .and. mcout1%pars(2) <= pi_xy%parmax(2) )
end subroutine 

subroutine test_run_parallel_mcmc_nchains4_len0(error)
  implicit none
  integer, parameter :: nchains = 4
  type(error_type), allocatable, intent(out):: error
  type(mcmc_output), dimension(:), allocatable :: mcout
  type(mcmc_output) :: mcout1
  type(mcmc_options):: mcopt  ! filled with defaults only
  integer :: i
  ! zero length run : takes expected input arguments, setup works, 
  ! outputs/writes unchanged state
  ! all on defaults, without optional arguments
  call init_pi()
  mcopt%maxiter = 0
  call run_parallel_mcmc(ll_normal, pi_xy, mcopt, mcout, nchains=nchains)
  ! expect values in mcout : a random initial state (within given parameter bounds) 
  ! and its loglikelihood
  do i=1, nchains
    mcout1 = mcout(i)
    call check(error, mcout1%pars(1) >= pi_xy%parmin(1) .and. mcout1%pars(1) <= pi_xy%parmax(1))
    call check(error, mcout1%pars(2) >= pi_xy%parmin(2) .and. mcout1%pars(2) <= pi_xy%parmax(2) )
  end do
end subroutine 


subroutine test_run_parallel_mcmc_nchains4_enforceomp_len0(error)
  implicit none
  integer, parameter :: nchains = 4
  type(error_type), allocatable, intent(out):: error
  type(mcmc_output), dimension(:), allocatable :: mcout
  type(mcmc_options):: mcopt  ! filled with defaults only
  type(mcmc_output) :: mcout1
  integer :: i
  ! zero length run : takes expected input arguments, setup works, 
  ! outputs/writes unchanged state
  ! all on defaults, without optional arguments
  call omp_set_num_threads(4)
  call init_pi()
  mcopt%maxiter = 0
  call run_parallel_mcmc(ll_normal, pi_xy, mcopt, mcout, nchains = nchains)
  ! expect values in mcout : a random initial state (within given parameter bounds) 
  ! and its loglikelihood
  do i=1, nchains
    mcout1 = mcout(i)
    call check(error, mcout1%pars(1) >= pi_xy%parmin(1) .and. mcout1%pars(1) <= pi_xy%parmax(1))
    call check(error, mcout1%pars(2) >= pi_xy%parmin(2) .and. mcout1%pars(2) <= pi_xy%parmax(2) )
  end do 
end subroutine 

subroutine test_run_parallel_mcmc_len1000(error)
  implicit none
  type(error_type), allocatable, intent(out):: error
  type(MCMC_OUTPUT):: MCOUT
  type(MCMC_OPTIONS):: MCOPT 
  ! all on defaults, without optional arguments
  MCOPT%MAXITER = 1000
  call run_mcmc(ll_normal, PI_xy, MCOPT, MCOUT)
  ! Expect x and y close to true values were found
  ! and best loglik is close to 0
  ! TODO "is close"  helper
  call check(error, MCOUT%pars(1), x_ideal )
  call check(error, MCOUT%pars(2), y_ideal)
  call check(error, MCOUT%ll, 0.0 )
  ! This is an easy problem, expect convergence was reached long before MAXITER.
end subroutine 


end module test_MCMC
