module samplers_shared

!> A collection of info about the model's parameters
!> n_pars and min, max bounds as two arrays
!> Closely related to model fct; model fct must take this number 
!> and type of paramters
!> Unlike in previous versions, intended to be intent(in) only
!> TODO we could bundle this in a type with the function:
!> the fortran quasi-object
type PARINFO
integer :: n_pars
double precision, allocatable, dimension(:) :: parmin, parmax, paradj
logical, allocatable, dimension(:):: fix_pars
end type PARINFO



contains 

!! utilities

pure logical function is_infinity(ll) 
!! Check whether a loglikelihood is infinity / likelihood is zero 
!! usually signalling hard reject of the state due to boundary conditions and physical constraints
double precision, intent(in) :: ll
double precision :: l
! undo the log 
l = exp(ll)
is_infinity = ( abs(l) <= tiny ) ! check approx zero 
end function

!! core sampler math

!> Accept or reject 
!> First argument is new/proposed log(!)likelihood, 
!> second is old log likelihood
!> return logical 
logical function metropolis_choice(new_loglikelihood, old_loglikelihood)  ! Really should be in or shared with MCMC
double precision, intent(in):: new_loglikelihood, old_loglikelihood   
double precision :: r  ! draw random number 0 to 1
call random_number(r)
! TODO add optional pregen random
! l1/l2 > r  <=> logl1-logl2 > log(r)
metropolis_choice = ((new_loglikelihood-old_loglikelihood) > log(r) )
end function

!!!! routines for random initialization 

subroutine init_pars_random_pregen(PI, pars0, fix_pars_flag, uniform_random_vector)
    use random_uniform, only : UNIF_VECTOR, next_random_uniform
    implicit none
    type(PARINFO), intent(in) :: PI ! give number, bounds of params
    double precision , dimension(PI%n_pars), intent(inout) :: pars0 ! return random initial values - nonnormalized
    logical, dimension(PI%n_pars), optional :: fix_pars_flag ! flags .true. to keep inidividual pars
    type(UNIF_VECTOR) :: uniform_random_vector ! object supplying pre-generated randoms 0 to 1
    integer :: i

    ! TODO check for need to generate random numbers vector

    if (.not. (present(fix_pars_flag) )) then
        fix_pars_flag(:) = .false.
    endif

    do i = 1, PI%n_pars
       ! parfix = 1 stay at prior value, parfix = 0 randomly search
        ! TODO condense this horrible logic to restart and fix_pars_flag somewhere outside of this fct
       !if (MCO%fixedpars .and. PI%parini(i) /= -9999d0) PI%parfix(i) = 1d0
       ! only assign random parameters if (a) randparini == .true. or (b) PI$parini(n) == -9999)
       !if (MCO%randparini .and. PI%parfix(i) == 0d0 .and. .not.restart_flag) then
!         
        !TODO vectorize
        if (.not. (fix_pars_flag(i) )) then
            ! make sure to give it a random number vector unique to the chain
            ! scale each to the parameter's range
            call log_nor2par(1, uniform_random_vector%next_random_uniform(), PI%parmin(i), PI%parmax(i), pars0(i)) 
       end if

    end do  ! for PI%npar loop
    !TODO test
    !TODO merge , optional uniform_random_vector arg
  end subroutine 

  subroutine init_pars_random(PI, pars0, fix_pars_flag)
    implicit none
    type(PARINFO), intent(in) :: PI !give number, bounds, and potentially current value of params
    double precision , dimension(PI%n_pars), intent(inout) :: pars0 ! return random initial values - nonnormalized
    logical, dimension(PI%n_pars), optional :: fix_pars_flag ! flags .true. to keep inidividual pars
    double precision :: r
    integer :: i
    if (.not. (present(fix_pars_flag) )) then
        fix_pars_flag(:) = .false.
    endif

    do i = 1, PI%n_pars
        if (.not. (fix_pars_flag(i) )) then
            call random_number(r) ! using the fortran intrinsic function
            call log_nor2par(1, r, PI%parmin(i), PI%parmax(i), pars0(i))  
       end if 
    end do  ! for PI%npar loop

  end subroutine 

  subroutine init_latin_square(PI, pars0, n_chains) !TODO 
    use math_functions, only: nor2par_scalar
    implicit none
    type(PARINFO), intent(in) :: PI !give number, bounds, and potentially current value of params
    integer, intent(in) :: n_chains
    double precision , dimension(PI%n_pars, n_chains), intent(out) :: pars0 ! return initial values - nonnormalized
    double precision, dimension(PI%n_pars,n_chains) :: points 
    integer :: i,j
    ! generate N initial points on the (0..1)^N space in latin hypercube distribution... 
    ! Must be run for all chains at once outside of parallel regions


    ! convert to real parameter values space
    do j=1, n_chains
        do i=1, PI%n_pars
            pars0(i,j) = nor2par_scalar( points(i,j), PI%parmin(i), PI%parmax(i))
        end do
    end do

  end subroutine 

  pure logical function bounds_check(PI, PARS)
    type(PARINFO), intent(in) :: PI 
    double precision, dimension(:), intent(in) :: PARS

  ! given real-space params ... convert to lognorm space, check if all in (0,1)
  ! or check directly against real boundary values
  bounds_check = all( (PARS > PI%parmin) .and. (PARS<PI%parmax) )

  end function

end module
