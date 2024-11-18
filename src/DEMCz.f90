module DEMCz_module

   !!!!!!!!!!!
   ! MC differential evolution z algorithm
   ! jklebes 2024
   ! Implementing ter Braak & Vrugt 2008
   !!!!

   implicit none
   
   ! TODO  expose for testing?
   public 

   public :: DEMCz

   !> A collection of input options to the DEMCz sampler run
   !> contains default values 
   type DEMCzOPT
      integer :: MAXITER ! overall steps, if convergence not reached
      integer :: n_steps ! steps per independent sampling period
      integer :: N_chains ! consider setting OMP env to something compatible
      real :: differential_weight = 0.8 ! differential weight f, [0,2]
      real :: crossover_probability = 0.9 ! crossover probability gamma , [0,1]
      real :: P_target ! termination criteria
   end type DEMCzOPT

   !> A collection of info about the model's parameters
   !> n_pars and min, max bounds as two arrays
   !> Closely related to model fct; model must take this number 
   !> and type of paramters
   !> Unlike in previous versions, intended to be intent(in) only
   !> TODO we could bundle this in a type with the function
   type PARINFO
      integer :: n_pars
      real, allocatable, dimension(:) :: parmin, parmax
   end type PARINFO

   !> Collection of info for output of the sampling rnu
   !> Note output is mainly via file writing
   type MCMC_OUTPUT
      real :: bestll
      real, allocatable, dimension(:) :: bestpars
   end type MCMC_OUTPUT

contains

   !> Main DEMCz sampler subroutine
   !> Write all OMP parallelization at this level only
   !> IN: @todo model_likelihood_write alternative loglikelihood for writing to file
   !> IN: model_loglikelihood a function parameters (array) -> loglikelihood (real)
   !> IN: PI type(ParInfo) collection of parameter bounds
   !> IN: OPT type(DEMCz) collection of sampling options
   !> OUT: MCOUT type(DEMCzOUT) collection of results
   !> Also writes history to file/output stream and progress to console.  
   subroutine DEMCz(model_likelihood_write, model_likelihood, PI, MCO, MCOUT)
    implicit none
      
       interface
      function model_likelihood_write(param_vector) result(ML) 
           implicit none
           real, dimension(:), intent(in):: param_vector
           real :: ML
      end function model_likelihood_write
        end interface
       interface
      function model_likelihood(param_vector) result(ML)
           implicit none
           real, dimension(:), intent(in):: param_vector
           real :: ML
      end function model_likelihood
        end interface
      

      !! inputs
      type(PARINFO), intent(in) :: PI
      type(DEMCzOPT), intent(in) :: MCO
      type(MCMC_OUTPUT), intent(out) :: MCOUT

      !> Matrix X , npars x nchains, holding current state of the n chains
      real, allocatable, dimension(:,:) :: PARS_current 
      ! and their current loglikelihood values
      real, allocatable, dimension(:) :: l0
      ! and their best likelihood values 
      real, allocatable, dimension(:) :: l_best
      real, allocatable, dimension(:,:) :: PARS_best 
 
      !> history Matrix Z , npars x maxiter 
      real, allocatable, dimension(:,:) :: PARS_history
      
      integer :: npars 
      integer :: nchains, nsteps, len_history, MAXITER
      integer :: i,j

      npars = PI%n_pars
      nchains = MCO%N_chains
      nsteps = MCO%N_chains
      MAXITER = MCO%MAXITER
      !TODO check the function takes npars arguments 

      allocate(PARS_current(npars, Nchains))
      ! but we need these to persist between parallel regions
      allocate(l0(nchains))
      allocate(l_best(nchains))
      allocate(PARS_best(npars,nchains))
      allocate(PARS_history(npars, MCO%MAXITER))

      ! number of entries in history matrix so far
      len_history = 0

!$    OMP PARALLEL DO
      do j=1, nchains
         ! choose initial values
         call init_random(PI, PARS_current(:,j), l0(j))

         ! potential burnin steps
         ! write first values to history matrix
         PARS_history(:,j) = PARS_current(:,j)
      end do
!$    OMP END PARALLEL DO
      len_history = len_history+nchains

      do i = 2, MAXITER
!$       OMP PARALLEL DO
         do j=1, nchains
            call step_chain(PARS_current(i, :), l0(i), PARS_history, model_likelihood, & 
            npars, len_history, nsteps)
         end do
         ! write to Z
         PARS_history(:,len_history+j) = PARS_current(:,i)
!$       OMP END PARALLEL DO !!Barrier implicit ?
      ! increment length M (filled so far) of Z
      len_history = len_history+nchains
         ! check convergence

         ! Reorder for best chains ?  Then write to Z later.
      end do

   end subroutine

   !> Initialize the chain's state with random values from 
   !> parameter ranges.
   !> This is not the best initialization; all later sampling is 
   !> bounded by min/max of the chains' random initial
   !> values, plus noise.
   !> Use a latin hypercube init ater.
   subroutine init_random(PI, norpars, ll)
      type(PARINFO), intent(in) :: PI
      real, dimension(:), intent(out) :: norpars
      real, intent(out) :: ll
      do i = 1, PI%npars
         call random_number(norpars(i)) 
      end do
   end subroutine

   !> Evolve the state of one chain by k steps
   subroutine step_chain(X_i, l0, PARS_history, model_likelihood, npars, len_history, n_steps)
    real, dimension(:), intent(inout) :: X_i
    real, dimension(:), allocatable :: vector , prop_vector !internal: save previous state, proposed new state
    real, dimension(:,:), intent(in) :: PARS_history !the matrix Z so far, to read 2 rows from 
    real, intent(inout) :: l0 ! likelihood of previous accepted params
    real             :: l !likelihood of proposed values
    integer, intent(in) :: npars, len_history 
    integer :: R1, R2 !indices of 2 random rows (may be same)
    real :: rand
    integer :: i, n_steps

    interface
    function model_likelihood(param_vector) result(ML)
         implicit none
         real, dimension(:), intent(in):: param_vector !TODO kind double precision ?
         real :: ML
    end function model_likelihood
      end interface

      do i=1, n_steps
      R1 = random_int(len_history)
      R2 = random_int(len_history) ! without replacement - possibly equal R1
      vector = X_i
      call step(prop_vector, vector, PARS_history(R1,:), PARS_history(R2,:))
      l= model_likelihood(prop_vector)
      if (metropolis_Choice(l, l0)) then !We should have this in MCMC common
         X_i = prop_vector
         l0=l
      end if 
      end do
   end subroutine

   subroutine step(vout, v1, v2, v3, gamma, f) 
    real, dimension(:), intent(out) :: vout
    real, dimension(:), intent(in) :: v1, v2, v3
    real, intent(in) :: gamma, f
    ! TODO 
    ! R1, R2 non-identical random indices
    ! mask arrays
      vout = v1 + gamma*(v2 -v3) ! + noise e
   end subroutine

   !> Accept or reject 
   !> First argument is new / proposed log(!) likelihood, 
   !> second is old log likelihood
   !> return logical 
   logical function metropolis_choice(new_loglikelihood, old_loglikelihood) ! Really should be in or shared with MCMC
      real, intent(in) :: new_loglikelihood, old_loglikelihood   
      real r ! draw random number 0 to 1
      call random_number(r)
      ! l1/l2 > r  <=> logl1 - logl2 > log(r)
      metropolis_choice = ((new_loglikelihood-old_loglikelihood) > log(r) )
   end function

   integer function random_int(N)
      integer, intent(in) :: N
      real :: r
      call random_number(r)
      random_int = floor(N*r)+1
   end function

end module DEMCz_module
