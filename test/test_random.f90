module test_random
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use random_uniform
  implicit none
  private

  public:: collect_randomtests

  integer, parameter :: dp = kind(0.0d0)

contains

!> Collect all exported unit tests
subroutine collect_randomtests(testsuite)
  !> Collection of tests
  type(unittest_type), allocatable, intent(out):: testsuite(:)

  testsuite = [ &
    new_unittest("fill random_uniform", test_fill_random_uniform),  &
    new_unittest("initialize random uniform", test_initialize),  &
    new_unittest("get_random_uniform", test_get_random_uniform),  &
    new_unittest("next_random_uniform", test_next_random_uniform) &
    ]

end subroutine collect_randomtests

subroutine test_fill_random_uniform(error)
  implicit none
  type(error_type), allocatable, intent(out):: error
  type(UNIF_VECTOR) :: random_uniform
  integer :: n 
  double precision, dimension(:), allocatable :: arr
  call random_uniform%initialize
  n = random_uniform%length
  allocate(arr(n))
  arr = 0d0
  call fill_random_uniform(arr, n)
  !write(*,*) arr
  call check(error, arr(1) > 0d0 .and. arr(1) <= 1.0 )
  call check(error, arr(5) > 0d0 .and. arr(5) <= 1.0 )
  call check(error, arr(n) > 0d0 .and. arr(n) <= 1.0 )
  call check(error, arr(2) /= arr(1) )
end subroutine 

subroutine test_initialize(error)
  implicit none
  type(error_type), allocatable, intent(out):: error
  type(UNIF_VECTOR) :: random_uniform
  call random_uniform%initialize
  call check(error, random_uniform%index , 1 )
  call check(error, allocated(random_uniform%u))
  call check(error, random_uniform%u(2) > 0d0 .and. random_uniform%u(2) <= 1.0 )
  call check(error, random_uniform%u(2) /= random_uniform%u(1) )
  call check(error, random_uniform%u(random_uniform%length) > 0d0 .and. random_uniform%u(random_uniform%length) <= 1.0 )
end subroutine

subroutine test_get_random_uniform(error)
  implicit none
  type(error_type), allocatable, intent(out):: error
  type(UNIF_VECTOR) :: random_uniform
  integer:: n
  double precision , dimension(:), allocatable :: x 
  call random_uniform%initialize
  n = 1
  x = random_uniform%get_random_uniform(n)
  call check(error, x(1) > 0d0 .and. x(1) <= 1.0 )
  n = 12
  x = random_uniform%get_random_uniform(n)
  call check(error, x(1) > 0d0 .and. x(1) <= 1.0 )
  call check(error, x(n) > 0d0 .and. x(n) <= 1.0 )
  !n = random_uniform%length +1 
  !x = random_uniform%get_random_uniform(n) ! TODO not handled,  error / crash expected
end subroutine

subroutine test_next_random_uniform(error)
  implicit none
  type(error_type), allocatable, intent(out):: error
  double precision  :: x 
  integer ::  index
  type(UNIF_VECTOR) :: random_uniform
  call random_uniform%initialize
  index = random_uniform%index
  x = random_uniform%next_random_uniform()
  call check(error, x > 0d0 .and. x <= 1.0 )
  call check(error, random_uniform%index , index + 1)
end subroutine

subroutine test_random_multivariate(error)
  implicit none
  type(error_type), allocatable, intent(out):: error
  ! make covariance matrix
  ! make npars
  ! make vector mu (mean) to center the step on 
  ! get vector rn
  !call random_multivariate(npars, 1, covariance_matrix, mu, rn)
end subroutine test_random_multivariate 


subroutine test_random_normal(error)
  implicit none
  type(error_type), allocatable, intent(out):: error
  ! make covariance matrix
  ! make npars
  ! make vector mu (mean) to center the step on 
  ! get output value result
  !call random_normal(result)
end subroutine test_random_normal


end module test_random
