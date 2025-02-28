module test_MCMC
  use testdrive, only : new_unittest, unittest_type, error_type, check
  implicit none
  private

  public:: collect_MCMCtests

contains

!> Collect all exported unit tests
subroutine collect_MCMCtests(testsuite)
  !> Collection of tests
  type(unittest_type), allocatable, intent(out):: testsuite(:)

  testsuite = [ &
    new_unittest("valid", test_valid) &
    ]

end subroutine collect_MCMCtests

subroutine test_valid(error)
  type(error_type), allocatable, intent(out):: error
  ! ...
end subroutine test_valid


!TODO move to common/math test
subroutine test_random_multivariate(error)
  type(error_type), allocatable, intent(out):: error
  ! make covariance matrix
  ! make npars
  ! make vector mu (mean) to center the step on 
  ! get vector rn
  call random_multivariate(npars, 1, covariance_matrix, mu, rn)
end subroutine test_random_multivariate 

!TODO move to common
subroutine test_random_normal(error)
  type(error_type), allocatable, intent(out):: error
  ! make covariance matrix
  ! make npars
  ! make vector mu (mean) to center the step on 
  ! get output value result
  call random_normal(result)
end subroutine test_random_normal

subroutine test_step(error)
  type(error_type), allocatable, intent(out):: error
  ! step() takes vector of pars, modifies it

end subroutine test_step

subroutine test_adapt_step_size(error)
  type(error_type), allocatable, intent(out):: error
  ! ...
  ! expect either a covariance matrix exists or pi%use_multivariate is set to false
  ! TODO should be called update_covariance
end subroutine test_adapt_step_size

subroutine test_increment_covariance_matrix(error)
  type(error_type), allocatable, intent(out):: error
  call increment_covariance_matrix(history, mean, N_local, N_history, cov)
  ! expect changes to mean, cov
  ! and n_local +=1 
end subroutine

subroutine test_cholesky_factor(error)
  type(error_type), allocatable, intent(out):: error
  call cholesky_factor(N, cov, posdef)
end subroutine

subroutine test_covariance_matrix(error)
  type(error_type), allocatable, intent(out):: error
  call covariance_matrix(history, mean, npars, N_history, cov)
  ! expect output at mean, cov
end subroutine

subroutine test_run_mcmc_len0(error)
  type(error_type), allocatable, intent(out):: error
  ! zero length run : takes expected input arguments, setup works, 
  ! outputs/writes unchanged state
  
  ! all on defaults, without optional arguments
  call run_mcmc(quadratic_loglik, PI, MCOPT, MCOUT)
  ! expect values in MCOUT : a random initial state (within given parameter bounds) 
  ! and its loglikelihood
  call check(error, MCOUT%ll > 0 )
  call check(error, MCOUT%p(1) >= x_lower .and. MCOUT%p(1) <=x_upper )
  call check(error, MCOUT%p(2) >= y_lower .and. MCOUT%p(2) <=y_upper )
  ! Given starting values which happen to be the correct (min loglik) solution, 
  ! Expect after zero run steps MCOUT has these initial values and 
  ! lolglik is low.
  call check(error, MCOUT%p(1), x_ideal )
  call check(error, MCOUT%p(2), y_ideal)
  call check(error, MCOUT% ll, 0.0 )
end subroutine 

subroutine test_run_mcmc_len1000(error)
  type(error_type), allocatable, intent(out):: error
  ! all on defaults, without optional arguments
  call run_mcmc(quadratic_loglik, PI, MCOPT, MCOUT)
  ! Expect x and y close to true values were found
  ! and best loglik is low.
  ! TODO "is close"  helper
  call check(error, MCOUT%p(1), x_ideal )
  call check(error, MCOUT%p(2), y_ideal)
  call check(error, MCOUT% ll, 0.0 )
  ! This is an easy problem, expect convergence was reached long before MAXITER.
end subroutine 

subroutine test_run_multichain_len0(error)
  type(error_type), allocatable, intent(out):: error
  ! zero length run : takes expected input arguments, setup works, 
  ! outputs/writes unchanged state
  !TODO set nchain and call main function 
  ! all on defaults, without optional arguments
  call run_mcmc(quadratic_loglik, PI, MCOPT, MCOUT)
  ! expect values in MCOUT : a random initial state (within given parameter bounds) 
  ! and its loglikelihood
  call check(error, MCOUT%ll > 0 )
  call check(error, MCOUT%p(1) >= x_lower .and. MCOUT%p(1) <=x_upper )
  call check(error, MCOUT%p(2) >= y_lower .and. MCOUT%p(2) <=y_upper )
  ! Given starting values which happen to be the correct (min loglik) solution, 
  ! Expect after zero run steps MCOUT has these initial values and 
  ! lolglik is low.
  call check(error, MCOUT%p(1), x_ideal )
  call check(error, MCOUT%p(2), y_ideal)
  call check(error, MCOUT% ll, 0.0 )
end subroutine 

subroutine test_run_multichain_len1000(error)
  type(error_type), allocatable, intent(out):: error
  ! all on defaults, without optional arguments
  call run_mcmc(quadratic_loglik, PI, MCOPT, MCOUT)
  ! Expect x and y close to true values were found
  ! and best loglik is low.
  ! TODO "is close"  helper
  call check(error, MCOUT%p(1), x_ideal )
  call check(error, MCOUT%p(2), y_ideal)
  call check(error, MCOUT% ll, 0.0 )
  ! This is an easy problem, expect convergence was reached long before MAXITER.
end subroutine 


end module test_MCMC
