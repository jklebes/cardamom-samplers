module samplers_shared

!> A collection of info about the model's parameters
!> n_pars and min, max bounds as two arrays
!> Closely related to model fct; model fct must take this number 
!> and type of paramters
!> Unlike in previous versions, intended to be intent(in) only
!> TODO we could bundle this in a type with the function:
!> the fortran quasi-object
type PARINFO
integer:: n_pars
real, allocatable, dimension(:):: parmin, parmax
end type PARINFO

!!!! routines for random initialization 
contains 

subroutine init_pars_random_pregen(PI, pars0, uniform_random_vector, uniform, fix_pars_flag)
    implicit none
    type(PARINFO), intent(in) :: PI !give number, bounds, and potentially current value of params
    double precision , dimension(PI%n_pars), intent(inout) :: pars0 ! return random initial values - nonnormalized
    logical, dimension(PI%n_pars), optional :: fix_pars_flag ! flags .true. to keep inidividual pars
    double precision , dimension(:), intent(in) :: uniform_random_vector ! pre-generated (can be faster) randoms 0 to 1
    integer, intent(inout) :: uniform ! our index in uniform_random_vector
    if (.not. (present(fix_pars_flag) )) then
        fix_pars_flag(:) = .false.
    endif

    do i = 1, PI%npars
       ! parfix = 1 stay at prior value, parfix = 0 randomly search
        ! TODO condense this horrible logic to restart and fix_pars_flag somewhere outside of this fct
       !if (MCO%fixedpars .and. PI%parini(i) /= -9999d0) PI%parfix(i) = 1d0
       ! only assign random parameters if (a) randparini == .true. or (b) PI$parini(n) == -9999)
       !if (MCO%randparini .and. PI%parfix(i) == 0d0 .and. .not.restart_flag) then
!         
        if (.not. (fix_pars_flag )) then
            !TODO bad for determinism - access module level shared random number vector.  make 
            ! it local to run_mcmc
            call nor2par(1, uniform_random_vector(uniform), PI%parmin(i), PI%parmax(i), pars0(i))  
            !call log_nor2par(1, uniform_random_vector(uniform), PI%parmin(i), PI%parmax(i), PI%paradj(i), PI%parini(i))
            uniform = uniform+1
       end if

    end do  ! for PI%npar loop

  end subroutine 

  subroutine init_pars_random(PI, pars0, fix_pars_flag)
    implicit none
    type(PARINFO), intent(in) :: PI !give number, bounds, and potentially current value of params
    double precision , dimension(PI%n_pars), intent(inout) :: pars0 ! return random initial values - nonnormalized
    logical, dimension(PI%n_pars), optional :: fix_pars_flag ! flags .true. to keep inidividual pars
    
    if (.not. (present(fix_pars_flag) )) then
        fix_pars_flag(:) = .false.
    endif

    do i = 1, PI%npars
        if (.not. (fix_pars_flag )) then
            !TODO bad for determinism - access module level shared random number vector.  make 
            ! it local to run_mcmc
            call nor2par(1, random_number(), PI%parmin(i), PI%parmax(i), pars0(i))  
       end if ( condition ) then
    end do  ! for PI%npar loop

  end subroutine 

  subroutine init_latin_square(PI, pars0, n_chains) !TODO 
    implicit none
    type(PARINFO), intent(in) :: PI !give number, bounds, and potentially current value of params
    interger, intent(in) :: n_chains
    double precision , dimension(PI%n_pars, n_chains), intent(out) :: pars0 ! return initial values - nonnormalized
   
    ! generate N initial points on the (0..1)^N space in latin hypercube distribution... 
    ! Must be run for all chains at once outside of parallel regions


    ! convert to real parameter values space
    pars0 = nor2par(PI%n_pars, points, PI%parmin(i), PI%parmax(i), pars0)

  end subroutine 


end module
