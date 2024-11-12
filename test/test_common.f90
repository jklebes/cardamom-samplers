module test_common
  use testdrive, only : new_unittest, unittest_type, error_type, check
  implicit none
  private

  public:: collect_commontests

contains

!> Collect all exported unit tests
subroutine collect_commontests(testsuite)
  !> Collection of tests
  type(unittest_type), allocatable, intent(out):: testsuite(:)

  testsuite = [ &
    new_unittest("valid", test_valid)  &
    ]

end subroutine collect_commontests

subroutine test_valid(error)
  type(error_type), allocatable, intent(out):: error
  ! ...
end subroutine test_valid

end module test_common
