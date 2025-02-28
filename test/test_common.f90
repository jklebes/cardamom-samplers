module test_common
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use samplers_shared
  implicit none
  private

  public:: collect_commontests

  integer, parameter :: dp = kind(0.0d0)

contains

!> Collect all exported unit tests
subroutine collect_commontests(testsuite)
  !> Collection of tests
  type(unittest_type), allocatable, intent(out):: testsuite(:)

  testsuite = [ &
    new_unittest("infini", test_is_infinity)  &
    ]

end subroutine collect_commontests

subroutine test_is_infinity(error)
  type(error_type), allocatable, intent(out):: error
  double precision :: P, log_P
  P = 0.0_dp
  log_P = log(P) 
  call check(error, is_infinity(log_P), .true.) 
  P = 0.01_dp
  log_P = log(P) 
  call check(error, is_infinity(log_P), .false.) 
end subroutine test_is_infinity

end module test_common
