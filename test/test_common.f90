module test_common
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use random_uniform
  use test_functions
  use samplers_shared
  implicit none
  private

  public:: collect_commontests

  integer, parameter:: dp = kind(0.0d0)

contains

!> Collect all exported unit tests
subroutine collect_commontests(testsuite)
  !> Collection of tests
  type(unittest_type), allocatable, intent(out):: testsuite(:)

  testsuite = [ &
    new_unittest("infini", test_is_infinity), &
    new_unittest("test init pars random", test_init_pars_random ) &
    ]

end subroutine collect_commontests

subroutine test_is_infinity(error)
  type(error_type), allocatable, intent(out):: error
  double precision:: P, log_P
  P = 0.0_dp
  log_P = log(P) 
  call check(error, is_infinity(log_P), .true.) 
  P = 0.01_dp
  log_P = log(P) 
  call check(error, is_infinity(log_P), .false.) 
end subroutine test_is_infinity

subroutine test_init_pars_random(error)
  type(error_type), allocatable, intent(out):: error
  logical, dimension(:), allocatable:: fix_pars_flag 
  double precision, dimension(:), allocatable:: pars0
  type(UNIF_VECTOR):: random_uniform
  call random_uniform%initialize

  if (.not. allocated(pars0)) allocate(pars0(PI_xy%npars))
  ! with no fix_pars_flag
  call init_pars_random(PI_xy, pars0, uniform_random_vector = random_uniform) 
  if (.not. allocated(fix_pars_flag) ) allocate(fix_pars_flag(PI_xy%npars))
  ! with fix_pars_flag all false
  fix_pars_flag = .false.
  call init_pars_random(PI_xy, pars0, fix_pars_flag, random_uniform) 
  ! With fix_pars_flag all true
  fix_pars_flag = .true.
  call init_pars_random(PI_xy, pars0, fix_pars_flag, random_uniform) 
end subroutine

end module test_common
