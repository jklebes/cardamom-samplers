module test_DEMCz
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use DEMCz_module
  implicit none
  private

  public:: collect_DEMCztests

contains

! A test function E = A*(x-x_0)^2 + B*(y-y_0)^2 
pure real function ll_normal(pars) result(res)
real, dimension(:), intent(in):: pars
real :: x, y ! the pars to fit
real :: x_0, y_0 ! The correct, energy/loglikelihood-minimizing answer will be x=x0, y=y0
real :: A, B 
x = pars(1)
y = pars(2)
x_0 = 5.1
y_0 = 5.0
A = 1.0
B = 1.6 ! covariance matrix expected to have inversely proportional entries on diagonal
! and zeros on off-diagonal for x-y correlation 
res = A*(x-x_0)**2 + B*(y-y_0)**2
end function

!> Collect all exported unit tests
subroutine collect_DEMCztests(testsuite)
  !> Collection of tests
  type(unittest_type), allocatable, intent(out):: testsuite(:)

  testsuite = [ &
    new_unittest("addition", test_test), &
    new_unittest("Options_type", test_MCO), &
    new_unittest("random_int", test_random_int), &
    new_unittest("metropolis_choice", test_metropolis_choice), &
    new_unittest("metropolis_stochastic", test_metropolis_stochastic) &
    ]

end subroutine collect_DEMCztests

subroutine test_test(error)
  implicit none
  type(error_type), allocatable, intent(out):: error

  call check(error, 1+2 == 3)
  if (allocated(error)) return

  ! equivalent to the above
  call check(error, 1+2, 3)
  if (allocated(error)) return
end subroutine test_test

subroutine test_MCO(error)
  type(error_type), allocatable, intent(out):: error
  !Having imported DEMCz_module, we should have a type(DEMCzOpt) 
  ! and access a module-level object holding default options
  type(DEMCzOPT) :: options
end subroutine

! TODO move to common
subroutine test_random_int(error)
  use DEMCz_module, only: random_int
  implicit none
  type(error_type), allocatable, intent(out):: error
  integer:: r
  r = random_int(2)
  !! should choose from (1, 2) (inclusive)
  call check(error, r == 1 .or. r==2)
end subroutine test_random_int

subroutine test_metropolis_choice(error)
!> Metropolis chioce function takes two log(!) likelihoods
!> and returns logical
!> The first argument is the new/proposed loglikelihood
  use DEMCz_module, only: metropolis_choice
  implicit none
  type(error_type), allocatable, intent(out):: error
  ! certain acceptance of state with probability 1 vs 0
  call check(error, metropolis_choice(log(1.0), log(1e-15)),.true. )
  ! certain rejection of state with probability 0 vs 1
  call check(error, metropolis_choice(log(1e-15), log(1.0)),.false. )
end subroutine test_metropolis_choice

subroutine test_metropolis_stochastic(error)
  !> Metropolis chioce function takes two log(!) likelihoods
  !> and returns logical
  !> The first argument is the new/proposed loglikelihood
    use DEMCz_module, only: metropolis_choice
    implicit none
    type(error_type), allocatable, intent(out):: error
    real :: new_loglikelihood, old_loglikelihood, accept_ratio
    integer :: i, N, accept_count
    ! Something with a likelihood l1 = 1/2 l2 should be acceped
    ! 50% of the time.  
    new_loglikelihood = log(.3) 
    old_loglikelihood = log(.6)
    ! get acceptance N times  - expect about 50% true 
    N= 500
    accept_count = 0
    do i=1, N
      if (metropolis_choice(new_loglikelihood, old_loglikelihood)) then
        accept_count = accept_count +1 
      end if 
    end do
    accept_ratio  = accept_count / real(N)
    write (*,*) accept_count
    write (*,*) accept_ratio

    call check(error, accept_ratio > .4  .and. accept_ratio < .6)
  end subroutine test_metropolis_stochastic

  subroutine test_DEMCz_runs(error)
    use DEMCz_module, only: DEMCz, PARINFO  
    implicit none
    type(error_type), allocatable, intent(out):: error
    ! test the main DEMCz function just runs when given a function
    ! it returns / modifies no info ; it writes to file

    type(PARINFO) :: PI
    type(DEMCzOPT) :: options
    type(MCMC_OUTPUT) :: DEMCzOUT 

    ! A struct with parameter limits
    PI%n_pars = 2
    PI%parmin(1) = -3.0
    PI%parmax(1) = 110
    PI%parmin(2) = 2.5
    PI%parmin(2) = 7.5
   
    ! single chain
    options%n_steps = 100
    options%MAXITER =1000
    options%N_chains = 1
    options%differential_weight = 0.8

    call DEMCz(ll_normal, PI, Options, DEMCzOUT)

    ! multi-chain
    options%N_chains = 4
    call DEMCz(ll_normal, PI, Options, DEMCzOUT)
    
  end subroutine test_DEMCz_runs

end module test_DEMCz
