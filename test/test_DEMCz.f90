module test_DEMCz
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use DEMCz_module
  implicit none
  private

  public:: collect_DEMCztests

contains

!> Collect all exported unit tests
subroutine collect_DEMCztests(testsuite)
  !> Collection of tests
  type(unittest_type), allocatable, intent(out):: testsuite(:)

  testsuite = [ &
    new_unittest("addition", test_test), &
    new_unittest("random_int", test_random_int), &
    new_unittest("metropolis_choice", test_metropolis_choice) &
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

subroutine test_random_int(error)
  use DEMCz_module, only: random_int
  implicit none
  type(error_type), allocatable, intent(out):: error
  integer:: r
  r = random_int(2)
  !! should choose from (1, 2) (inclusive)
  call check(error, r > 0)
  call check(error, r < 3)
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

end module test_DEMCz
