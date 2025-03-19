module test_math
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use math_functions
  implicit none
  private

  public:: collect_mathtests

  integer, parameter :: dp = kind(0.0d0)

  ! some example arguments for conversion
  integer, parameter :: npars = 5
  double precision, dimension(npars) :: parmax = [1d0, 0d0, 1.20d0, -13.5d0, 0.00000043d0]
  double precision , dimension(npars) :: parmin = [0.00d0, -0.3d0, -1.20d0,  -190.5d0, 0.00000042d0]
  double precision , dimension(npars) :: paradj = [1.0d0, 1.3d0, 2.20d0,  191.5d0, 0.0d0] ! abs(parmin1 if nonpositive 
  double precision , dimension(npars) :: pars1 = [1.0d0, -.01d0, 0.0d0, -100d0, 0.000000421d0]
  double precision , dimension(npars) :: pars_norm1 = [0.0d0, 0.99d0, 0.32d0, 1.0d0, 0.000001d0]
contains

!> Collect all exported unit tests
subroutine collect_mathtests(testsuite)
  !> Collection of tests
  type(unittest_type), allocatable, intent(out):: testsuite(:)

  testsuite = [ &
    new_unittest("approx", test_approx),  &
    new_unittest("par2nor", test_par2nor),  &
    new_unittest("nor2par", test_nor2par),  &
    new_unittest("log_par2nor", test_log_par2nor),  &
    new_unittest("log_nor2par", test_log_nor2par)  &
    ]

end subroutine collect_mathtests

! testing util approx equal
function approx(a, b, rel_tol) result(eq)
  implicit none
  double precision , intent(in) :: a, b
  double precision, intent(in), optional :: rel_tol
  double precision :: rel_tol_default = 1.0e-8
  double precision :: rel_tol_
  logical :: eq
  if (.not. present(rel_tol)) then
    rel_tol_ = rel_tol_default
  else
    rel_tol_ = rel_tol
  endif 
  ! check close to zero
  if (a==0d0) then 
    eq = (abs(b) <= epsilon(0d0))
  else if (b==0d0) then
    eq = (abs(a) <= epsilon(0d0))
  else 
  ! else check relative difference
    eq = (abs(a-b) <= rel_tol_*abs(a))  .or.   (abs(a-b) <= rel_tol*abs(b)) 
  endif
end function


subroutine test_approx(error)
  !! test my approx equals testing utility
  implicit none
  type(error_type), allocatable, intent(out):: error
  call check(error, approx(1d0, 0d0), .false.)
  call check(error, approx(0d0, -1d0), .false.)
  call check(error, approx(0d0, 0d0), .true.)
  call check(error, approx(27d0/9d0, 3d0), .true.)
end subroutine 

subroutine test_par2nor(error)
  implicit none
  type(error_type), allocatable, intent(out):: error
  double precision, dimension(npars) :: pars_norm, pars2
  integer :: i
  pars_norm = par2nor(npars, pars1, parmin, parmax)
  ! expect values between 0 and 1
  call check(error, all(pars_norm>=0d0) .and. all(pars_norm<=1d0))
  ! try converting back 
  pars2 = nor2par(npars, pars_norm, parmin, parmax)
  ! expect approximately equal to the original
  call check(error, all ([(approx(pars1(i), pars2(i)), i=1, npars)] ))
end subroutine 


subroutine test_nor2par(error)
  implicit none
  type(error_type), allocatable, intent(out):: error
  double precision, dimension(npars) :: pars_real
  integer :: i
  pars_real = nor2par(npars, pars_norm1, parmin, parmax)
  ! expect values in bounds parmin, parmax
  call check(error, all(pars_real>=parmin) .and. all(pars_real<=parmax))
end subroutine 


subroutine test_log_par2nor(error)
  implicit none
  type(error_type), allocatable, intent(out):: error
  double precision, dimension(npars) :: pars_norm, pars2
  integer :: i
  pars_norm = log_par2nor(npars, pars1, parmin, parmax, paradj)
  ! expect values between 0 and 1
  call check(error, all(pars_norm>=0d0) .and. all(pars_norm<=1d0))
  ! try converting back 
  pars2 = log_nor2par(npars, pars_norm, parmin, parmax, paradj)
  ! expect approximately equal to the original
  call check(error, all ([(approx(pars1(i), pars2(i)), i=1, npars)] ))
end subroutine 


subroutine test_log_nor2par(error)
  implicit none
  type(error_type), allocatable, intent(out):: error
  double precision, dimension(npars) :: pars_real
  pars_real = log_nor2par(npars, pars_norm1, parmin, parmax, paradj)
  ! expect values in bounds parmin, parmax
  call check(error, all(pars_real>=parmin) .and. all(pars_real<=parmax))

end subroutine 

end module test_math
