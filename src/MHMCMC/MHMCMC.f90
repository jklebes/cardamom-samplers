module MHMCMC

  !!!!!!!!!!!
  ! Authorship contributions
  !
  ! This code is based on the original C verion of the University of Edinburgh
  ! CARDAMOM framework created by A. A. Bloom (now at the Jet Propulsion Laboratory).
  ! All code translation into Fortran, integration into the University of
  ! Edinburgh CARDAMOM code and subsequent modifications by:
  ! T. L. Smallman (t.l.smallman@ed.ac.uk, University of Edinburgh)
  ! J. F. Exbrayat (University of Edinburgh)
  ! See function/subroutine specific comments for exceptions and contributors
  !!!!!!!!!!!

  ! refactored jklebes 2024-25

  ! Module contains all subroutine and functions relevant specifically to the
  ! AP-MCMC method. The choice of EDC, likelihood and model are made else where and
  ! are thus contains within a seperate module

  ! Relevant source references:
  ! Haario et al., (2001) An adaptive Metropolis algorithm. Bernoulli 7.2: 223-242.
  ! Haario et al., (2006) Stat. Comput., 16:339–354, DOI 10.1007/s11222-006-9438-0, 
  ! Roberts and Rosenthal (2009), Examples of Adaptive MCMC, J. Comp. Graph. Stat. 18:349-367

   ! On normalized parameters: 
   ! "raw" values, as convenient for loglikelihood calculation and file writing, 
   ! are the default  For function arguments, internal saved data.  Parameters
   ! are converted to normalized values as needed by adaptive statistics functions.  History matrix 
   !  PARSALL and covariance matrix work refer to the normalized parameter space.
   !
   ! How to use this sampler:
   !  Create an object of type PARINFO : number and bounds of parameters 
   !  and DEMCzOPT-options, containing as many or few of the fields as needed, the rest default to the 
   !                 default values in type definiton here
   !  Create an object of type MCMC_OUTPUT to write reults to.
   !  Create a real function loglikelihood taking a vector of npars (same as in PARINFO) parameters and 
   !                 returning rel:: loglikelihood.
   !  (not implemented yet) Optionally set OMP_NUM_THREADS
   !  Call subroutine DEMCz(fct, parinfo, mcopt, mcmcout) 

use samplers_shared, only: PARINFO
use samplers_io, only: io_buffer_space, initialize_buffers, open_output_files
use OMP_LIB

implicit none

public

!> A collection of input options to the DEMCz sampler run
!> contains default values 
type MCMC_OPTIONS
integer:: MAXITER = 10000  ! overall steps, if convergence not reached
integer:: nadapt  =1000  ! steps per "local" sampling period, between adaptation steps
integer:: N_chains = 1  ! consider setting OMP env to something compatible
integer:: nwrite = 1000
integer:: nprint = 1000
real:: P_target  ! termination criteria 
! file names
character(350):: outfile = "parout.txt"
character(350):: stepfile = "stepout.txt"
character(350) ::  covfile = "covout.txt"
character(350):: covifile = "covinfoout.txt"
! Adaptive !!
! setting for adaptive AP-MCMC step size
double precision:: par_minstepsize = 0.001d0 & ! 0.0005 -> 0.001 -> 0.01 -> 0.1 -> 0.005
                              ,par_maxstepsize = 0.01d0  &
                              ,par_initstepsize = 0.005d0
double precision:: beta = 0.05d0  ! weighting for gaussian step in multivariate proposals
! Optimal scaling variable for parameter searching
double precision:: N_before_mv_target, & !
                    opt_scaling_const = 2.381204**2  ! scd = 2.381204 the optimal scaling parameter
                                ! for MCMC search, when applied to  multivariate proposal.
                                ! NOTE 1: 2.38/sqrt(npars) sometimes used when applied to the Cholesky
                                ! factor. NOTE 2: 2.381204**2 = 5.670132
double precision:: N_before_mv = 10d0
!! step
! Is current proposal multivariate or not?
logical:: multivariate_proposal = .false.
real:: fadapt  ! TODO fraction adapt-move to outsied
integer:: nout 
logical:: append
logical:: use_multivariate
logical:: restart
logical:: randparini
logical:: returnpars  ! a variable that is never used and has no effect, needs deleting in all model likelihood files
logical, dimension(:), allocatable:: fixedpars 
end type MCMC_OPTIONS


!> Collect ongoing statistics (means, covariance matrix, etc) for each chain
type MCSTATS
end type


!> Collection of info for output of the sampling run
!> Note output is mainly via file writing
type MCMC_OUTPUT
double precision:: bestll, ll
double precision, allocatable, dimension(:):: bestpars, pars
double precision:: acceptance_rate
logical:: complete 
integer:: nos_iterations
!stats collection: 
double precision:: Nparvar, Nparvar_local 
double precision, allocatable, dimension(:):: parvar, meanpar
double precision, allocatable, dimension(:,:):: covariance
logical:: cov = .false. ! Does the covariance matrix exist yet?
logical:: use_multivariate
logical:: multivariate_proposal
end type MCMC_OUTPUT

contains
  !
  !--------------------------------------------------------------------
  !
  subroutine run_parallel_mcmc( model_likelihood, PI, MCO, MCOUT_list, model_likelihood_write_in, restart_in, nchains)

    implicit none

    !/* ***********INPUTS************
    ! *
    ! * MODEL_LIKELYHOOD: A function wholly responsible for
    ! * (a) running the model given the DATA and parameters, 
    ! * (b) comparing it to observations, and
    ! * (c) returning  the (log) likelihood.
    ! * The function will be run as MODEL_LIKELIHOOD(PARS)
    ! * and returns a single log(!)likelihood value
    ! *
    ! * PARINFO: This structure contains information on
    ! * (b) initpars:        parameter starting values (optional/recommended).
    ! * (c) npars:           number of pars
    ! *
    ! * MCO: This structure contains option values for the MCMC run.
    ! * These will be set to default values if empty, restart, nchains . Options include:
    ! * (a) number of runs
    ! * (b) filename for writing file with results
    ! * (c) step adaptation frequency
    ! * (d) initial step size
    ! * */
    !
    !/* **************OUTPUTS*************
    ! *
    ! * RESULTS FILE: File includes (a) results (b) likelihood and (c) final step
    ! size
    ! *
    ! * */
    ! Write to MCMCOUT

      !! input and output structs
      ! read-only, shared beteen chains:
      type(PARINFO), intent(in):: PI
      type(MCMC_OPTIONS), intent(in):: MCO
      ! output, one for each chain
      type(MCMC_OUTPUT), dimension(:), allocatable, intent(inout):: MCOUT_list  ! array of MCOUT objects

      logical, optional, intent(in):: restart_in  ! is it a restart ? (i.e. start from data in MCOUT instead of initializing new)
      logical:: restart
      !type(MCMC_OUTPUT), optional, intent(in):: MCOUT_prev  ! TODO array ? ! process
      integer, optional:: nchains
      integer:: nchains_
      integer:: i

      ! the function to minimize  ! TODO change name to loglikelihood everywhere
      ! Completely agnostic, samples any functions pars -> loglikelihood
      interface
    subroutine model_likelihood(param_vector, n, ML)
         implicit none
         double precision, dimension(n), intent(inout):: param_vector  ! intent(in), inout for compatibility with R via C
         integer, intent(in):: n
         double precision, intent(out):: ML
    end subroutine model_likelihood
    end interface! 
    !optionally  give a second function with same shape as model_likelihood, 
      ! for writing to file; 
    procedure(model_likelihood), optional:: model_likelihood_write_in
    ! Going forwards this pointer is the alternateive likelihood function for writing to file
    procedure(model_likelihood), pointer:: model_likelihood_write
      

      ! Argument processing  !!!!!!!!!!!!!!!

      if (.not. present(restart_in)) then 
        restart = .false.
      else
        restart = restart_in
      endif

      if (.not. present(nchains)) then
        nchains_ = 1 
        ! or try to infer from OMP_THREADS or columns in MCOUT_PREV ...
      else
        nchains_ = nchains
      endif 

      ! process function arguments
      if (present(model_likelihood_write_in)) then
         ! if second function given use it for writing to file
         model_likelihood_write => model_likelihood_write_in
      else
         ! default, if no argument given: use same function for calculation and printing
         model_likelihood_write => model_likelihood
      end if

    ! TODO check the function(s) take npars arguments !

    ! Outputs
      ! each one an MCOUT struct with unitialized output data fields
      if (.not. allocated(MCOUT_list)) allocate(MCOUT_list(nchains_)) 


    ! This adaptive MCMC is "parallelized" to do N chains at the same time, but they
    ! each do a complete independent run.  The parallelization structure of this module is 
    ! trivial.  It exists mainly as template for samplers with more crossover and more complex
    ! structure in this loop.
    !$OMP parallel do
    do i = 1, nchains_  ! TODO where the user can input this
        ! saves best loglikelihood and associated parameters to MCOUT(i)
        ! MCOUT_list(i) possibly contains desired starting poisition in `pars` field, 
        ! possible entire history and stats from previous run, possible empty new MCOUT object
        call run_mcmc(model_likelihood, PI, MCO, MCOUT_list(i), model_likelihood_write, restart, i)
    end do
    !$OMP end parallel do

    !! compare chains

    !! Write files, write to output struct

  end subroutine run_parallel_mcmc

  ! Main function for a single adaptive MCMC simulation
  ! Global settings of the sampler taken from MCO 
  subroutine run_mcmc(model_likelihood, PI, MCO, MCOUT, model_likelihood_write_in, restart_in, chainid)
    use samplers_math, only: log_par2nor, log_nor2par, par2nor, nor2par
    use samplers_shared, only : init_pars_random, bounds_check, is_infinity, metropolis_choice
    use samplers_io, only: write_parameters, write_variances, write_covariance_matrix &
                          ,write_covariance_info, restart_flag, write_mcmc_output, open_output_files
    use random_uniform, ONLY: UNIF_VECTOR, initialize
    ! declare any local variables
    !! all variables in here are local to the single chain and the duration of its run
          !! input and output structs
    type(PARINFO), intent(in):: PI
    type(MCMC_OPTIONS), intent(in):: MCO
    !type(MCMC_OUTPUT), intent(in), optional:: MCOUT_prev
    type(MCMC_OUTPUT), intent(inout):: MCOUT
    !type(MCMC_OUTPUT), optional, intent(in):: prev_MCOUT  ! output of previous run, for restart  ! read from module data
    logical, intent(in), optional:: restart_in
    integer, intent(in), optional:: chainid
    logical:: restart

    type(io_buffer_space):: io_space  ! this chain has its own io buffers
    character(350):: outfile, stepfile, covfile, covifile
    double precision, dimension(PI%npars):: PARS_previous         & ! parameter values for current state
                                            ,PARS_proposed          & ! parameter values for current proposal
                                            ,BESTPARS        ! best set of parameters so far

    double precision, dimension(PI%npars, MCO%nadapt):: PARSALL  ! All accepted normalised parameters since previous step adaption
    double precision:: loglikelihood_previous, loglikelihood_proposed
        !! loglikelihood of a set of parameters 
    double precision:: output_loglikelihood
        !! loglikelihood according to  alternative loglikelihood calculation for writing to file
    double precision:: llmax
        !! best loglikelihood seen so far
    double precision:: burn_in_period & ! global TODO module data ?! for how many proposals will we adapt the covariance matrix as a minimum
                       , AM_likelihood &
                       ,outputP0      &
                       ,outputP0prior  &
                       ,Pmax, P0prior, Pprior & ! as below but for priors only
                       ,P0 & ! previously accepted observation based log-likelihood
                       ,P  & ! current observation based log-likelihood
                       ,opt_scaling & ! = opt_scaling_const/npars
                       ,beta &
                       ,par_minstepsize &
                       , P_target
    type(UNIF_VECTOR):: uniform_random_vector
    logical:: multivariate
      !! object holding array of pre-generated random values - (supposedly faster to pregenerate) - local to this 
      !! chain 
    integer:: i
    integer:: MAXITER, nchains, npars
    ! counters-local to this chain's run
    integer:: ITER, ACC, ACC_FIRST, ACCLOC
    double precision:: ACCRATE, ACCRATE_GLOBAL


    ! declare interface for the model likelihood function.
    ! A function pars vector -> loglikelihood
    ! Do not modify here, sampler is generic.
    ! Wrap the model loglikelihood function in another function elsewhere to 
    ! make it conform to this form.

    !! Collecting stats for this run, on space of npars normalized parameters.
    ! Formerly part of PI struct, but separated as belonging more to internals and output than input
    !double precision:: Nparvar 
    !double precision, dimension(PI%npars):: parvar, meanpar
    !double precision, dimension(PI%npars, PI%npars):: convariance

    interface
    subroutine model_likelihood(param_vector, n, ML)
         implicit none
         double precision, dimension(n), intent(inout):: param_vector  ! intent(in), inout for compatibility with R via C
         integer, intent(in):: n
         double precision, intent(out):: ML
    end subroutine model_likelihood
    end interface

      ! optionally give a second function with same shape as model_likelihood, 
      ! for writing to file; 
    procedure(model_likelihood), optional:: model_likelihood_write_in
    ! Going forwards this pointer is the alternateive likelihood function for writing to file
    procedure(model_likelihood), pointer:: model_likelihood_write


    logical:: accept


      if (present(model_likelihood_write_in)) then
         ! if second function given use it for writing to file
         model_likelihood_write => model_likelihood_write_in
      else
         ! default, if no argument given: use same function for calculation and printing
         model_likelihood_write => model_likelihood
      end if
    ! read settings from input struct ...
    npars = PI%npars
    nchains = MCO%N_chains
    MAXITER = MCO%MAXITER
    P_target = MCO%P_target
    MCOUT%use_multivariate = MCO%use_multivariate
    beta = MCO%beta
    par_minstepsize = MCO%par_minstepsize

    ! initialize output fields
    MCOUT%Nparvar = 0
    allocate(MCOUT%parvar(npars))
    allocate(MCOUT%meanpar(npars))
    allocate(MCOUT%covariance(npars, npars))

    ! process file names 
    outfile = MCO%outfile
    stepfile = MCO%stepfile
    covfile = MCO%covfile
    covifile = MCO%covifile
    if (MCO%n_chains > 1 .and. present(chainid) ) then
      write (outfile, '(a, i3)') MCO%outfile, chainid
      write (stepfile, '(a, i3)') MCO%stepfile, chainid
      write (covfile, '(a, i3)') MCO%covfile, chainid
      write (covifile, '(a, i3)') MCO%covifile, chainid
    endif


    !! Calculate derived  settings of the run... 
    ! Determine how long we will continue to adapt our proposal covariance
    ! matrix and use of Delayed Rejection
    burn_in_period = MCO%fADAPT*dble(MAXITER)
    ! See step() for relevant references.
    ! scd = 2.381204 the optimal scaling parameter for MCMC search, when applied
    ! to multivariate proposal.
    ! NOTE 1: 2.38/sqrt(npars) sometimes used when applied to the Cholesky factor
    ! NOTE 2: 2.381204**2 = 5.670132
    opt_scaling = MCO%opt_scaling_const/dble(PI%npars)
  
    if (.not. present(restart_in)) then 
      restart = .false.
    else
      restart = restart_in
    endif

    if (restart) then
      ! keep MCOUT
    else
      ! init MCOUT-could be a function
      MCOUT%complete = .false.
      MCOUT%nos_iterations = 0
      if (.not. allocated(MCOUT%pars)) allocate(MCOUT%pars(npars))
      if (.not. allocated(MCOUT%bestpars)) allocate(MCOUT%bestpars(npars))
    endif 

    ! Initialize basic MCMC counters
    ITER = 0

    ! Initialize further counters for adaptive
    ACC = 0
    ACC_first = 0
    ACCLOC = 0
    ACCRATE = 0d0
    ACCRATE_GLOBAL = 0d0

    ! Initialize pregenerated random numbers, if using-local to this chain
    call uniform_random_vector%initialize() 
    
    !TODO opt scaling scalin by n

    !!! prepare file writing
    ! allocate buffers (different one for each chain)
    call initialize_buffers(npars, MAXITER/MCO%nwrite, io_space)
    ! TODO potential restart handling !  outside 
    !call check_for_existing_output_files(npars, nOUT, nWRITE, sub_fraction &
    !, parname, stepname, covname, covinfoname)
    !TODO open separate output file for each chain !
    call open_output_files(MCO%outfile, MCO%stepfile, MCO%covfile, MCO%covifile)



    !!!! init params
    ! TODO implement reading in PI%fix_pars
    if (.not. restart) then
    call init_pars_random(PI, PARS_previous, PI%fix_pars, uniform_random_vector)
    ! Inform the user
    write(*,*) "Have loaded/randomly assigned PI%parini-now begin the AP-MCMC"
    ! initialize loglikelihood value of the given model with these pars
    !P = -1d0; Pprior = -1d0
    ! calculate the initial probability/log likelihood.
    ! NOTE: passing P0 -> P is needed during the EDC searching phase where we
    ! could read an EDC consistent parameter set in the first instance
    call model_likelihood(PARS_previous, npars, loglikelihood_previous); 

    write(*,*) loglikelihood_previous
    if (is_infinity(loglikelihood_previous)) then  ! TODO is this log evaluating? TODO better check for hard reject ll = inifinity P = 0
        write(*,*) "WARNING  ! loglikelihood = ",loglikelihood_previous, " - &
        & AP-MCMC will get stuck, if so please check initial conditions"
        error stop
    endif
    endif 

    ! initalize bestpars to current pars
    BESTPARS = PARS_previous
    llmax = loglikelihood_previous


    ! Begin the main AP-MCMC loop
    do while (ITER < MAXITER)! .and. Pmax < P_target)

       ! take a step in parameter space: generate proposed 
       ! new parameters PARS 
       ! should include reflectivenedd/redrawing 
       multivariate = MCOUT%use_multivariate .and. (MCOUT%Nparvar > MCO%N_before_mv_target)
       call step_pars_real(PARS_previous, PARS_proposed, PI, multivariate, MCOUT%covariance, beta, opt_scaling, par_minstepsize, &
          uniform_random_vector) 
       ! if parameter proposal in bounds check the model
       ! TODO LATER if we get a reflective gaussian kernel, no need to check
       if (bounds_check(PI, PARS_proposed)) then
           ! calculate the model likelihood
           ! TODO ideally output just one likelihood value
           call model_likelihood(PARS_proposed, npars, loglikelihood_proposed)
           accept = metropolis_choice(loglikelihood_proposed, loglikelihood_previous)
       else
           accept = .false.
           ! TODO LATER should we count the out-of-boundary in adaptation?
           ! Keep equivalent for now !  Change after checking the whole sample for equivalence to previous
       end if  ! in bound

       if (accept) then

           ! Store accepted parameter proposals (unnormalized values)
           ! keep record of all parameters accepted since step adaption
           ! (this chain)
           ! Because this history matrix is used for (normalized) statistics for adaptiveness, 
           ! store normalized version of pars
            PARSALL(1:npars, ACCLOC+1) = log_par2nor(npars, PARS_proposed, PI%parmin, PI%parmax, PI%paradj)  ! add row in history matrix
           ! store the best parameter set
           if (loglikelihood_proposed >= llmax) then
               BESTPARS = PARS_proposed; llmax = loglikelihood_proposed
           endif
           ! Keep count of the number of accepted proposals in this local period
           ACCLOC = ACCLOC+1
           ! Accepted first proposal from multivariate
           if (MCOUT%multivariate_proposal) ACC_first = ACC_first+1 !!TODO what is this counter ?
           ! TODO may have to do with whether to calc covariance matrix

           PARS_previous(1:npars) = PARS_proposed(1:npars)          ! save as previous pars
           !norPARS0(1:PI%npars) = norPARS(1:PI%npars)                    ! normalize
           loglikelihood_previous = loglikelihood_proposed;               ! save as previous loglikelihood
        else  ! TODO
            ! write to history
       endif  ! accept or reject proposed pars

       ! count iteration 
       ITER = ITER+1

       if (MCO%nwrite > 0) then
           ! TODO fct
           if (mod(ITER, MCO%nwrite) == 0) then
!              print*,"mcmc: write_mcmc_output done"

               ! calculate the likelhood for the actual uncertainties-this avoid
               ! issues with different phases of the MCMC which may use sub-samples
               ! of observations or inflated uncertainties to aid parameter
               ! searching
               call model_likelihood_write(PARS_proposed, npars, output_loglikelihood)
               ! Now write out to files
               call write_mcmc_output(MCOUT%parvar, ACCRATE, &
                                      MCOUT%covariance, &
                                      MCOUT%meanpar, MCOUT%Nparvar, &
                                      PARS_previous, output_loglikelihood, npars, ITER == MCO%nOUT, io_space)
           end if 
       end if  ! write or not to write

       ! time to adapt?
       if (mod(ITER, MCO%nadapt) == 0) then
            ! TODO fct
!           ! Debugging print statements
!           print*,"mcmc: time to adapt"
           
           !! update the acceptance counters and acceptance ratios
           ! Total accepted values
           ACC = ACC+ACCLOC
           ! Calculate global acceptance rate
           ACCRATE_GLOBAL = ACC/ITER

           ! Calculate local acceptance rate (i.e. since last adapt)
           ACCRATE = ACCLOC/dble(MCO%nadapt)

           ! Second, are we still in the adaption phase?
           ! TODO how does fortran integer division work
           if (burn_in_period > ITER .or. (ACC_first/ITER) < 0.05d0 .or. .not.MCO%use_multivariate) then

               ! Once covariance matrix has been created just update based on a
               ! single parameter set from each period.
            ! TODO ??
               if (MCOUT%cov) then
                   ACCLOC = 1
                   PARSALL(1:npars, ACCLOC) = log_par2nor(npars, PARS_previous, PI%parmin, PI%parmax, PI%paradj)
                   ! leads to call to increment_covariance_matrix with the one new row
               else if (ACCLOC > 3) then
                   PARSALL(1:npars, 2) = PARSALL(1:npars, ceiling(ACCLOC*0.5d0))
                   PARSALL(1:npars, 3) = PARSALL(1:npars, ACCLOC)
                    ACCLOC = 3
                    ! leads to attempting to make new covariance matrix from the first, middle, and last rows of local history
               endif

               ! adapt the covariance matrix for multivariate proposal
               ! TODO rename "parsall" to replect it's really a small subsample of period's history
               call update_statistics(PARSALL, npars, MCOUT, MCOUT%use_multivariate, ACCLOC, MCO%N_before_mv_target)

           end if !  have enough parameter been accepted
           ! TODO what if MCO%use_multivariate ???

           ! resets the local acceptance counter
           ACCLOC = 0

       end if  ! time to adapt?

       ! Should I be write(*,*)ing to screen or not?
       if (MCO%nPRINT > 0) then
           if (mod(ITER, MCO%nPRINT) == 0) then
               write(*,*)"Using multivariate sampling = ",MCO%use_multivariate
               write(*,*)"Total proposal = ",ITER, " out of ",MCO%nOUT
               write(*,*)"Total accepted = ",ACC
               write(*,*)"Overall acceptance rate    = ",dble(ACC)/dble(ITER)
               write(*,*)"Local   acceptance rate    = ",ACCRATE
               write(*,*)"Current obs   = ",P0, "proposed = ",P, " log-likelihood"
               write(*,*)"Current prior = ",P0prior, "proposed = ",Pprior, " log-likelihood"
               write(*,*)"Maximum likelihood = ",Pmax
               ! NOTE: that-infinity in current obs only indicates failure of EDCs
               ! but-infinity in both obs and parameter likelihood scores indicates
               ! that proposed parameters are out of bounds
           end if 
       end if  ! write(*,*) to screen or not

    end do  ! while conditions

    !!!  Finalize 

    ! write out final covariance matrix for the analysis
    if (MCO%nwrite > 0) call write_covariance_matrix(MCOUT%covariance, npars, .false.)

    ! record the best single set of parameters
    MCOUT%bestpars = BESTPARS
    MCOUT%bestll = llmax
    ! Output the current state, so that next simulations can potentially resume here
    MCOUT%pars = PARS_previous  !! (!) TODO should be done outside 
    MCOUT%ll = loglikelihood_previous
    ! record how many iterations were taken to complete
    write(*,*) MCOUT%nos_iterations, ITER
    MCOUT%nos_iterations = MCOUT%nos_iterations+ITER
    ! set flag MCMC completed
    MCOUT%complete = .true.
    ! tidy up
    !deallocate(uniform_random_vector)
    ! deallocate this chain's history-but should be fine on subroutine exit

    ! completed AP-MCMC loop
    write(*,*)"AP-MCMC loop completed"
    write(*,*)"Overall acceptance rate     = ", dble(ACC)/dble(ITER)
    write(*,*)"Final local acceptance rate = ",ACCRATE
    ! TODO function to output these two 
    write(*,*)"Best log-likelihood = ",Pmax
    !write(*,*)"Best parameters = ",MCOUT%bestpars

end subroutine

  !
  !------------------------------------------------------------------
  !
  subroutine update_statistics(PARSALL, npars, MCOUT, use_multivariate, ACCLOC, N_before_mv_target)
    use samplers_io, only: write_covariance_matrix, write_covariance_info
    use samplers_math, only: nor2par, par2nor, log_nor2par, log_par2nor, &
                              cholesky_factor, std, covariance_matrix, &
                              increment_covariance_matrix

    ! Update the multivariate propsal distribution.
    ! Ensure that this subroutine is only called if at least 1 parameter propsal
    ! has been accepted in the last adaption period.

    implicit none

    ! declare input types
    type(MCMC_OUTPUT), intent(inout):: MCOUT  ! incl statistics collection
    logical, intent (inout):: use_multivariate
    ! declare inputs variables
    integer, intent(in):: npars
    integer, intent(in):: ACCLOC
    double precision, intent(in):: PARSALL(npars, ACCLOC)  ! collection of recently accepted normalised parameter combinations
    !! two CARDAMOM quirks here :  1) only accepted steps were recorded, not repeat entries for non-accepted steps
    !!                             2) update_statistics is passed only first or first, middle, and last rows instead of whole history, extremely limnited covariance estimate

    ! declare local variables
    integer p, i, info  ! counters
    double precision, dimension(npars, npars):: cov_backup
    double precision, dimension(npars):: meanpar_backup
    double precision:: Nparvar_backup, Nparvar_local
    double precision:: N_before_mv_target

    ! if we have a covariance matrix then we want to update it, if not then we need to create one
    if (MCOUT%cov) then

        ! Increment the variance-covariance matrix with new accepted parameter sets
        ! NOTE: that this also increments the total accepted counter (PI%Nparvar)

        cov_backup = MCOUT%covariance; meanpar_backup = MCOUT%meanpar; Nparvar_backup = MCOUT%Nparvar

!        call increment_covariance_matrix(PARSALL(1:PI%npars, 1:nint(N%ACCLOC)), PI%meanpar, PI%npars &
!                                        ,PI%Nparvar, nint(N%ACCLOC), PI%covariance)
        ! Have started hardcoding a maximum number of observations to be 100.
        ! While not strictly following Haario et al., (2001) or Roberts and Rosenthal, (2009)
        ! this allows for the covariance matrix to be more responsive to its local environment.
        !! TODO discuss 
        Nparvar_local = min(N_before_mv_target, Nparvar_backup)
        call increment_covariance_matrix(PARSALL(1:npars, 1:ACCLOC), MCOUT%meanpar, npars &
                                        ,Nparvar_local, ACCLOC, MCOUT%covariance)
        ! Calculate the cholesky factor as this includes a determination of
        ! whether the covariance matrix is positive definite.
        call cholesky_factor( npars, MCOUT%covariance, info )
        ! If the updated covariance matrix is not positive definite we should
        ! reject the update in favour of the existing matrix
        ! TODO ??
        if (info == 0) then
            ! Set multivariate sampling to true
            use_multivariate = .true.
        else
            ! The current addition of a parameter leads to a matrix which is not
            ! positive definite. If we previously had a matrix which is positive
            ! definite then we should reject totally the new matrix, if not we
            ! should keep it and accumulate the information
            if (use_multivariate) then
                ! return original matrix to place
                MCOUT%covariance = cov_backup
                MCOUT%meanpar = meanpar_backup
                MCOUT%Nparvar = Nparvar_backup
            else
                ! Keep accumulating use_multivariatethe information
                use_multivariate = .false.
            end if
        endif

    else  ! PI%cov == .false.

        ! we have not yet created a covariance matrix based on accepted
        ! parameters. Assuming we have some then create one...
        if (ACCLOC > 2) then

            ! estimate covariance matrix
            call covariance_matrix(PARSALL(1:npars, 1:ACCLOC), MCOUT%meanpar, &
            & npars, ACCLOC, MCOUT%covariance)
            MCOUT%cov = .true. ; MCOUT%Nparvar = ACCLOC

            ! Calculate the cholesky factor as this includes a determination of
            ! whether the covariance matrix is positive definite.
            call cholesky_factor ( npars, MCOUT%covariance, info )
            ! If not positive definite then we should not use the multivariat
            
            ! step at this time.
            if (info /= 0) then
                ! Keep accumulating information until positive definite matrix
                ! calculated
                use_multivariate = .false.
            else


                ! Positive definite found straight away-might as well use it!
                use_multivariate = .true.
            endif

            ! write out first covariance matrix, this will be compared with the final covariance matrix
            !if (MCO%nWRITE > 0) then
                call write_covariance_matrix(MCOUT%covariance, npars, .true.)
                call write_covariance_info(MCOUT%meanpar, MCOUT%Nparvar, npars)
            !endif

        end if  ! N%ACCLOC > 2

    end if  ! PI%cov == .true.

    return

  end subroutine update_statistics

  ! Generates new proposed state from currect state in real parameter space.  
  ! Wraps step_pars
  subroutine step_pars_real(PARS0, PARS, PI, multivariate, covariance, beta, opt_scaling, &
     par_minstepsize, random_uniform_vector)
    use samplers_math, only: log_par2nor, log_nor2par
    use random_uniform, only : UNIF_VECTOR
    implicit none
    double precision, dimension(:), intent(in   ):: pars0    ! current parameters
    double precision, dimension(:), intent(out  ):: pars       ! proposal
    type(UNIF_VECTOR), intent(inout):: random_uniform_vector
    !integer, intent(in):: npars
    type(PARINFO):: PI
    !type(MCSTATS), intent(inout):: stats
    logical, intent(in):: multivariate
    double precision, dimension(:,:), intent(in):: covariance
    double precision, dimension(PI%npars)             :: pars0_norm, pars_norm
    double precision, intent(in):: beta, opt_scaling, par_minstepsize
    pars0_norm = log_par2nor(PI%npars, pars0, PI%parmin, PI%parmax, PI%paradj)
    call step_pars(pars0_norm, pars_norm, PI%npars, multivariate, covariance, beta, opt_scaling, par_minstepsize, &
         random_uniform_vector )
    pars = log_nor2par(PI%npars, pars_norm, PI%parmin, PI%parmax, PI%paradj)
  end subroutine 

  !
  !------------------------------------------------------------------
  !
  ! Applies Roberts and Rosenthal 2009-Eq 3 to generate new proposed state in (normalized)
  ! parameter space
  ! IN : Pars0 current state (normalized)
  ! IN : Stats MCOSTATS object holding info like covariance matrix of the run so far (calculated on 
  ! normalized space )
  ! OUT: PARS new proposed state (normalized)
  ! plus take beta from module data
  subroutine step_pars(PARS0, PARS, npars, multivariate, covariance, beta, opt_scaling, &
    par_minstepsize, random_uniform_vector)  ! TODO check against original !! 
    use samplers_math, only:  random_normal, random_multivariate
    use random_uniform, only : UNIF_VECTOR

    ! carries out the next step to parameters in the MCMC search

    implicit none

    ! declare input variables
    !double precision, dimension(PI%npars), intent(inout):: !norpars0 & ! normalised current parameters
                                                           !,norpars  & ! normalised proposal
    double precision, dimension(:), intent(in   ):: pars0    ! current parameters
    double precision, dimension(:), intent(out  ):: pars       ! proposal
    integer, intent(in):: npars
    type(UNIF_VECTOR), intent(inout):: random_uniform_vector
    !type(MHMCMCOPT), intent(in):: MCO
    !type(MCSTATS), intent(in):: stats
    logical, intent(in):: multivariate
    double precision, dimension(:,:), intent(in):: covariance
    double precision, intent(in):: beta, opt_scaling, par_minstepsize
    
    ! declare local variables
    integer:: p
    double precision:: rn(npars), mu(npars), rn2(npars)

    ! mean of distributions 
    mu = 0d0 

    ! Splitting step calculation based on number of parameter vectors accepted
    ! is linked to the need build a covariance matrix prior to multivariate
    ! sampling.
    ! See Roberts and Rosenthal, Examples of Adaptive MCMC, J. Comp. Graph. Stat. 18:349-367, 2009.
    
    ! Sample random normal distribution (mean = 0, sd = 1)
    ! TODO vectorize
    do p = 1, npars
        call random_normal(random_uniform_vector, rn2(p))
    end do

    if (multivariate) then !((MCOUT%use_multivariate .and. MCOUT%Nparvar > N_before_mv_target)) then

        ! Is this step a multivariate proposal or not
        ! MCOUT%multivariate_proposal = .true. ! this only affects ACC_first counter  ! TODO move to where use_multivariate, Nparvar updates

        ! Draw from multivariate random distribution
        ! NOTE: if covariance matrix provided is not positive definite
        !       a sample from normal distribution is returned
        call random_multivariate(npars, 1, covariance, mu, rn, random_uniform_vector)

        ! Estimate the step to be applied to the current parameter vector to
        ! create the new proposal. scd = a scaling parameter linking searching
        ! stepping to the number of parameters being retrieved by the analysis.
        ! See Haario et al., (2001) An adaptive Metropolis algorithm. Bernoulli 7.2: 223-242.
        ! and references therein. See also, Roberts & Rosenthal (2009) for beta scaling.
        ! TODO par_minstepsize a user input with a default value 
        pars = pars0 + (rn*opt_scaling * (1d0-beta)) + (par_minstepsize*rn2*beta)

    else 

      !MCOUT%multivariate_proposal = .false.
      pars = pars0 + (par_minstepsize*rn2)

    end if

  end subroutine step_pars
  !
  !------------------------------------------------------------------
  !
end module MHMCMC
