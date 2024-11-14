module test_common
  use testdrive, only : new_unittest, unittest_type, error_type, check
  !use common
  implicit none
  private

  public:: collect_commontests

contains

!> Collect all exported unit tests
subroutine collect_commontests(testsuite)
  !> Collection of tests
  type(unittest_type), allocatable, intent(out):: testsuite(:)

  testsuite = [ &
    new_unittest("infini", test_infini)  &
    ]

end subroutine collect_commontests

subroutine test_infini(error)
  type(error_type), allocatable, intent(out):: error
  ! we have a constant infini = log(0d0), 
  ! And I expect it was set and 
  ! holds a big number (not error, NaN, underflow, etc)
  ! TODO compare against size of e?? 
  !check(error, infini > 1e10)
end subroutine test_infini

end module test_common
