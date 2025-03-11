module test_math
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use math_functions
  implicit none
  private

  public:: collect_mathtests

  integer, parameter :: dp = kind(0.0d0)

contains

!> Collect all exported unit tests
subroutine collect_mathtests(testsuite)
  !> Collection of tests
  type(unittest_type), allocatable, intent(out):: testsuite(:)

  !testsuite = [ &
    !new_unittest("infini", test_is_infinity)  &
  !  ]

end subroutine collect_mathtests


end module test_math
