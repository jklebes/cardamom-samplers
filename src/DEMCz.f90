module DEMCz_module

   !!!!!!!!!!!
   ! MC differential evolution z algorithm
   ! jklebes 2024
   ! Implementing ter Braak & Vrugt 2008
   !
   ! On normalized parameters: 
   ! All calculation, statistics, and writing is done on "raw" values, 
   ! NOT parameters normalized to (0,1).  
   !
   ! How to use:
   !  Create an object of type PARINFO : number and bounds of parameters 
   !  and DEMCzOPT - options, containing as many or few of the fields as needed, the rest default to the 
   !                 default values in type definiton here
   !  Create an object of type MCMC_OUTPUT to write reults to.
   !  Create a real function loglikelihood taking a vector of n_pars (same as in PARINFO) parameters and 
   !                 returning rel::loglikelihood.
   !  (not implemented yet) Optionally set OMP_NUM_THREADS
   !  Call subroutine DEMCz(fct, parinfo, demczopt, mcmcout) 
   !!!!
   use samplers_shared
   implicit none
   
   ! TODO  expose for testing?
   public 

   public:: DEMCz

   ! TODO add file writing

   !> Settings for this run, as module data
   real :: differential_weight

   !> A collection of input options to the DEMCz sampler run
   !> contains default values 
   type DEMCzOPT
      integer:: MAXITER  ! overall steps, if convergence not reached
      integer:: n_steps  ! steps per independent sampling period
      integer:: N_chains  ! consider setting OMP env to something compatible
      real:: differential_weight = 0.8  ! differential weight gamma, [0, 2]
      real:: crossover_probability = 0.9  ! crossover probability CR, [0, 1]
      real:: P_target  ! termination criteria
   end type DEMCzOPT


   !> Collection of info for output of the sampling run
   !> Note output is mainly via file writing
   type MCMC_OUTPUT
      real:: bestll
      real, allocatable, dimension(:):: bestpars
   end type MCMC_OUTPUT

contains

   !> Main DEMCz sampler subroutine
   !> Write all OMP parallelization at this top level only
   !> IN: model_loglikelihood a function parameters (array) -> loglikelihood (real)
   !> IN: model_loglikelihood_write: an alternative function of same type for writing to file
   !>                                (optional, defaults to using model_loglikelihood )
   !> IN: PI type(ParInfo) collection of parameter bounds
   !> IN: OPT type(DEMCz) collection of sampling options
   !> OUT: MCOUT type(DEMCzOUT) collection of results
   !> Also writes history to file/output stream and progress to console.  
   subroutine DEMCz(model_likelihood, PI, MCO, MCOUT, model_likelihood_write_in)
    implicit none

      !! input and output structs
      type(PARINFO), intent(in):: PI
      type(DEMCzOPT), intent(in):: MCO
      type(MCMC_OUTPUT), intent(out):: MCOUT

      ! the function to minimize ! TODO change name to loglikelihood everywehere
      interface
      function model_likelihood(param_vector) result(ML)
           implicit none
           real, dimension(:), intent(in):: param_vector
           real:: ML
      end function model_likelihood
      end interface
      
      ! optionally give a second function with same shape as model_likelihood , 
      ! for writing to file;
      procedure(model_likelihood), optional :: model_likelihood_write_in
      ! Going forwards this pointer is the alternateive likelihood function for writing to file
      procedure(model_likelihood), pointer :: model_likelihood_write

      !> Matrix X, (npars x nchains), holding current state of the n chains
      real, allocatable, dimension(:,:):: PARS_current 
      ! and their current loglikelihood values
      real, allocatable, dimension(:):: l0 
      ! and their best likelihood values and best pars so far
      real, allocatable, dimension(:):: l_best
      real, allocatable, dimension(:,:):: PARS_best 
 
      !> history Matrix Z, (npars x (nchains * maxiter))
      real, allocatable, dimension(:,:):: PARS_history

      
      integer:: npars, nchains, MAXITER , Ksteps
      integer:: i, j, k, len_history ! counters

      ! Argument processing  !!!!!!!!!!!!!!!
      ! Set functions ...

      if (present(model_likelihood_write_in)) then
         ! if second function given use it for writing to file
         model_likelihood_write => model_likelihood_write_in
      else
         ! default , if no argument given: use same function for calculation and printing
         model_likelihood_write => model_likelihood
      end if
 
      ! TODO check the function(s) take npars arguments !

      ! wrap the model loglikelihood function , which takes raw parameter values, in a new function
      ! which takes normalized values


      ! Extract from types ...
      differential_weight=MCO%differential_weight

      npars = PI%n_pars
      nchains = MCO%N_chains
      MAXITER = MCO%MAXITER
      Ksteps = MCO%n_steps

      ! Allocate arrays

      allocate(PARS_current(npars, Nchains))
      ! but we need these to persist between parallel regions
      allocate(l0(nchains))
      allocate(l_best(nchains))
      allocate(PARS_best(npars, nchains))
      allocate(PARS_history(npars, MCO%MAXITER*nchains))

      ! where we are in filling in the history matrix so far
      len_history = 0

      !!! Initial state

!$    OMP PARALLEL DO
      do j = 1, nchains
         ! choose initial values
         ! TODO better function for initial state : latin square
         call init_random(npars, PARS_current(:,j))
         ! also set the loglikelihoof of the state generated
         l0(j) = model_likelihood(PARS_current(:,j))

         ! potential burnin steps
         ! ... TODO

         ! write first values to history matrix
         PARS_history(:,j) = PARS_current(:,j)
      end do
!$    OMP END PARALLEL DO

      len_history = len_history+nchains

      !!! Main "time" loop

      do i = 2, MAXITER

         ! evolve each chain independently for nsteps (nsteps=K in ter Braak & Vrugt)
!$       OMP PARALLEL DO
         do j = 1, nchains
            do k = 1, Ksteps
               call step_chain(PARS_current(:, j), l0(j), model_likelihood, & 
            npars, PARS_history, len_history)
            end do
            ! write the chain's state after nsteps to Z
            PARS_history(:,len_history+j) = PARS_current(:,i)
         end do
!$       OMP END PARALLEL DO !!Barrier implicit ?


         ! increment length M (filled so far) of Z
         len_history = len_history+nchains
         
         ! check convergence

         ! Reorder for best chains ?  Then write to Z later.
      end do

      !!! Write out 

   end subroutine

   !> Initialize the chain's state with random values from 
   !> parameter ranges.
   !> This is not the best initialization; all later sampling is 
   !> bounded by min/max of the chains' random initial
   !> values, plus noise.
   !> Works with normalized values : returns a number between 0 and 1
   !> For each parameter
   subroutine init_random(npars, norpars)
      integer , intent(in) :: npars
      real, dimension(:), intent(out):: norpars
      integer :: i
      do i = 1, npars
         call random_number(norpars(i)) 
      end do
      ! also output the loglikelihood of the state generated

   end subroutine

   !> Evolve the state of one chain for 1 step
   subroutine step_chain(X_i, l0, model_likelihood, & 
      npars, PARS_history, len_history)
     integer, intent(in):: npars ! number of pars
    real, dimension(:), intent(inout):: X_i ! current state of the chain ; normalized values of all pars
    real, dimension(:), allocatable:: previous_vector, proposed_vector  ! internal: save previous state, proposed new state
    real, dimension(:,:), intent(in):: PARS_history  ! the matrix Z so far, to read 2 rows from
    integer, intent(in) :: len_history ! length to which Z is filled
    real, intent(inout):: l0  ! likelihood of previous accepted params
    real             :: l  ! likelihood of proposed values
    
    integer:: R1, R2  ! indices of 2 random rows 
    real:: rand
    integer:: i ! counters

    ! the function, vector of normalized par values -> loglikelihood
    interface
    function model_likelihood(param_vector) result(ML)
         implicit none
         real, dimension(:), intent(in):: param_vector  ! TODO kind double precision ?
         real:: ML
    end function model_likelihood
      end interface

      R1 = random_int(len_history)
      R2 = random_int(len_history)  
      do while(R1==R2) ! should not equal R1 ("without replacement")
         R2 = random_int(len_history) 
      end do
      previous_vector = X_i
      call step(proposed_vector, previous_vector , PARS_history(R1, :), PARS_history(R2, :))
      l = model_likelihood(proposed_vector)
      if (metropolis_choice(l, l0)) then  ! We should have this in MCMC common
         X_i = proposed_vector
         l0 = l
      end if 
   end subroutine

   ! Generate new proposed state from currect state and history
   ! ter Braak & Vrugt eq 2
   subroutine step(vout, v1, v2, v3) 
    real, dimension(:), intent(out):: vout
    real, dimension(:), intent(in):: v1, v2, v3
    ! get differential_weight, corssover_probability from module data
      vout = v1+differential_weight*(v2-v3) ! TODO + noise e
   end subroutine

   !> Accept or reject 
   !> First argument is new/proposed log(!)likelihood, 
   !> second is old log likelihood
   !> return logical 
   logical function metropolis_choice(new_loglikelihood, old_loglikelihood)  ! Really should be in or shared with MCMC
      real, intent(in):: new_loglikelihood, old_loglikelihood   
      real r  ! draw random number 0 to 1
      call random_number(r)
      ! l1/l2 > r  <=> logl1-logl2 > log(r)
      metropolis_choice = ((new_loglikelihood-old_loglikelihood) > log(r) )
   end function

   integer function random_int(N)
      integer, intent(in):: N
      real:: r
      call random_number(r)
      random_int = floor(N*r)+1
   end function

end module DEMCz_module
