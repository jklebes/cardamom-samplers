module MHMCMC_MODULE

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
   !  and DEMCzOPT - options, containing as many or few of the fields as needed, the rest default to the 
   !                 default values in type definiton here
   !  Create an object of type MCMC_OUTPUT to write reults to.
   !  Create a real function loglikelihood taking a vector of n_pars (same as in PARINFO) parameters and 
   !                 returning rel::loglikelihood.
   !  (not implemented yet) Optionally set OMP_NUM_THREADS
   !  Call subroutine DEMCz(fct, parinfo, mcopt, mcmcout) 

use samplers_shared, only: PARINFO

implicit none

! make all private
private

! specify what can be seen
public:: MHMCMC

!> A collection of input options to the DEMCz sampler run
!> contains default values 
type MHMCMCOPT
integer:: MAXITER  ! overall steps, if convergence not reached
integer:: n_steps  ! steps per "local" sampling period, between adaptation steps
integer:: N_chains  ! consider setting OMP env to something compatible
integer :: n_write 
integer :: nprint
real:: P_target  ! termination criteria 
! Adaptive !!
! setting for adaptive AP-MCMC step size
double precision :: par_minstepsize = 0.001d0 & ! 0.0005 -> 0.001 -> 0.01 -> 0.1 -> 0.005
                              ,par_maxstepsize = 0.01d0  &
                              ,par_initstepsize = 0.005d0
double precision:: beta = 0.05d0  ! weighting for gaussian step in multivariate proposals
! Optimal scaling variable for parameter searching
double precision:: N_before_mv_target, & !
                    opt_scaling_const = 2.381204**2  ! scd = 2.381204 the optimal scaling parameter
                                ! for MCMC search, when applied to  multivariate proposal.
                                ! NOTE 1: 2.38/sqrt(npars) sometimes used when applied to the Cholesky
                                ! factor. NOTE 2: 2.381204**2 = 5.670132
!! step
! Is current proposal multivariate or not?
logical:: multivariate_proposal = .false.
real :: fadapt ! TODO ?? misc variables from here
real :: nout 
logical :: append
logical :: use_multivariate
logical :: restart
end type MHMCMCOPT


!> Collection of info for output of the sampling run
!> Note output is mainly via file writing
type MCMC_OUTPUT
real:: bestll
real, allocatable, dimension(:):: bestpars, pars_final
integer :: complete !TODO logical
integer :: nos_iterations
end type MCMC_OUTPUT

!> Collect ongoing statistics (means, covariance matrix, etc) for each chain
type MCSTATS
double precision :: Nparvar, Nparvar_local 
double precision, allocatable, dimension(:) :: parvar, meanpar
double precision, allocatable , dimension(:,:) :: covariance
logical :: cov = .false. ! Does the covariance matrix exist yet?
logical :: use_multivariate
logical :: multivariate_proposal
end type


contains
  !
  !--------------------------------------------------------------------
  !
  subroutine MHMCMC( model_likelihood, PI, MCO, MCOUT, model_likelihood_write_in, restart, MCOUT_prev, nchains)

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
      type(PARINFO), intent(in):: PI
      type(MHMCMCOPT), intent(in):: MCO
      type(MCMC_OUTPUT), intent(out):: MCOUT

      logical, optional, intent(inout) :: restart !is it a restart ? (i.e. start from data in MCOUT instead of initializing new)
      type(MCMC_OUTPUT), optional, intent(in):: MCOUT_prev
      integer, optional :: nchains

      ! the function to minimize ! TODO change name to loglikelihood everywhere
      ! Completely agnostic, samples any functions pars -> loglikelihood
      interface
    subroutine model_likelihood(param_vector, n, ML)
         implicit none
         double precision, dimension(:), intent(in):: param_vector
         integer, intent(in) :: n
         double precision , intent(out):: ML
    end subroutine model_likelihood
    end interface! 
    !optionally  give a second function with same shape as model_likelihood , 
      ! for writing to file;
    procedure(model_likelihood), optional :: model_likelihood_write_in
    ! Going forwards this pointer is the alternateive likelihood function for writing to file
    procedure(model_likelihood), pointer :: model_likelihood_write
      

      ! Declare some arrays
      real, allocatable, dimension(:,:):: PARS_current 
      ! and their current loglikelihood values
      real, allocatable, dimension(:):: l0 
      ! and their best likelihood values and best pars so far
      real, allocatable, dimension(:):: l_best
      real, allocatable, dimension(:,:):: PARS_best 
 
   

      ! counters
      integer:: i

      ! Argument processing  !!!!!!!!!!!!!!!

      if (.not. present(restart)) then 
        restart = .false.
      endif

      if (.not. present(nchains)) then
        nchains = 1 
        ! or try to infer from OMP_THREADS or columns in MCOUT_PREV ...
      endif 

      ! Set functions ...

      if (present(model_likelihood_write_in)) then
         ! if second function given use it for writing to file
         model_likelihood_write => model_likelihood_write_in
      else
         ! default , if no argument given: use same function for calculation and printing
         model_likelihood_write => model_likelihood
      end if

    ! TODO check the function(s) take npars arguments !


    ! Allocate arrays - the ones in commom between parallel chains


    ! This adaptive MCMC is "parallelized" to do N chains at the same time, but they
    ! each do a complete independent run.  The parallelization structure of this module is 
    ! trivial.  It exists mainly as template for samplers with more crossover and more complex
    ! structure in this loop.
    !$OMP parallel do
    do i = 1, nchains !TODO where the user can input this
        ! saves best loglikelihood and associated parameters to MCOUT(i)
        call run_mcmc(model_likelihood, PI, MCO, MCOUT, model_likelihood_write, restart)
    end do
    !$OMP end parallel do

    !! compare chains

    !! Write files , write to output struct

  end subroutine MHMCMC

  ! Main function for a single adaptive MCMC simulation
  ! Global settings of the sampler taken from MCO 
  subroutine run_mcmc(model_likelihood, PI, MCO, MCOUT, model_likelihood_write_in, restart)
    use math_functions, only: log_par2nor, log_nor2par, par2nor, nor2par
    use samplers_shared, only : init_pars_random_pregen, bounds_check, is_infinity, metropolis_choice
    use cardamom_io, only: write_parameters, write_variances, write_covariance_matrix &
                          ,write_covariance_info, restart_flag, write_mcmc_output
    use random_uniform, ONLY: UNIF_VECTOR, initialize
    ! declare any local variables
    !! all variables in here are local to the single chain and the duration of its run
          !! input and output structs
    type(PARINFO), intent(in):: PI
    type(MHMCMCOPT), intent(in):: MCO
    type(MCMC_OUTPUT), intent(out):: MCOUT
    !type(MCMC_OUTPUT), optional, intent(in):: prev_MCOUT !output of previous run, for restart ! read from module data
    logical, intent(inout), optional :: restart

    double precision, dimension(PI%n_pars) :: PARS_previous         & ! parameter values for current state
                                            ,PARS_proposed          & ! parameter values for current proposal
                                            ,BESTPARS        ! best set of parameters so far

    double precision, dimension(PI%n_pars, MCO%n_steps):: PARSALL  ! All accepted normalised parameters since previous step adaption
    double precision :: loglikelihood_previous, loglikelihood_proposed
        !! loglikelihood of a set of parameters 
    double precision :: output_loglikelihood
        !! loglikelihood according to  alternative loglikelihood calculation for writing to file
    double precision :: llmax
        !! best loglikelihood seen so far
    double precision:: burn_in_period & ! global TODO module data ?! for how many proposals will we adapt the covariance matrix as a minimum
                       , AM_likelihood &
                       ,outputP0      &
                       ,outputP0prior  &
                       ,Pmax, P0prior, Pprior & ! as below but for priors only
                       ,P0 & ! previously accepted observation based log-likelihood
                       ,P  & ! current observation based log-likelihood
                       ,opt_scaling & ! = opt_scaling_const / npars
                       ,beta &
                       ,par_minstepsize &
                       , P_target
    type(UNIF_VECTOR) :: uniform_random_vector
      !! object holding array of pre-generated random values - (supposedly faster to pregenerate) - local to this 
      !! chain 
    integer:: i
    integer :: MAXITER, nchains, npars
    ! counters - local to this chain's run
    integer :: ITER, ACC, ACC_FIRST, ACCLOC
    ! acceptance ratios derived from counters  !TODO possibly belong in stats struct
    double precision :: ACCRATE, ACCRATE_GLOBAL


    ! declare interface for the model likelihood function.
    ! A function pars vector -> loglikelihood
    ! Do not modify here, sampler is generic.
    ! Wrap the model loglikelihood function in another function elsewhere to 
    ! make it conform to this form.

    !! Collecting stats for this run, on space of n_pars normalized parameters.
    ! Formerly part of PI struct, but separated as belonging more to internals and output than input
    !double precision :: Nparvar 
    !double precision, dimension(PI%n_pars) :: parvar, meanpar
    !double precision, dimension(PI%N_pars,PI%n_pars) :: convariance
    type (MCSTATS) :: stats

    interface
    subroutine model_likelihood(param_vector, n, ML)
         implicit none
         double precision, dimension(:), intent(in):: param_vector
         integer, intent(in) :: n
         double precision , intent(out):: ML
    end subroutine model_likelihood
    end interface

      ! optionally give a second function with same shape as model_likelihood , 
      ! for writing to file;
    procedure(model_likelihood), optional :: model_likelihood_write_in
    ! Going forwards this pointer is the alternateive likelihood function for writing to file
    procedure(model_likelihood), pointer :: model_likelihood_write


    logical :: accept


      if (present(model_likelihood_write_in)) then
         ! if second function given use it for writing to file
         model_likelihood_write => model_likelihood_write_in
      else
         ! default , if no argument given: use same function for calculation and printing
         model_likelihood_write => model_likelihood
      end if

    ! read settings from input struct ...
    npars = PI%n_pars
    nchains = MCO%N_chains
    MAXITER = MCO%MAXITER
    P_target = MCO%P_target
    stats%use_multivariate = MCO%use_multivariate
    beta = MCO%beta
    par_minstepsize = MCO%par_minstepsize


    !! Calculate derived  settings of the run... 
    ! Determine how long we will continue to adapt our proposal covariance
    ! matrix and use of Delayed Rejection
    burn_in_period = MCO%fADAPT*dble(MCO%nOUT)
    ! See step() for relevant references.
    ! scd = 2.381204 the optimal scaling parameter for MCMC search, when applied
    ! to multivariate proposal.
    ! NOTE 1: 2.38/sqrt(npars) sometimes used when applied to the Cholesky factor
    ! NOTE 2: 2.381204**2 = 5.670132
    opt_scaling = MCO%opt_scaling_const/dble(PI%n_pars)
  
      if (.not. present(restart)) then 
        restart = .false.
      endif

    if (restart) then
        ! read current state, statistics, history, etc from prev_MC_OUT
    else

    ! Initialize basic MCMC counters
    ITER = 0

    ! Initialize further counters for adaptive
    ACC = 0
    ACC_first = 0
    ACCLOC = 0
    ACCRATE = 0d0
    ACCRATE_GLOBAL = 0d0

    ! Initialize pregenerated random numbers, if using - local to this chain
    call uniform_random_vector%initialize() 
    
    !TODO opt scaling scalin by n

    !!! prepare file writing
    !call prepare_file()
    ! add something here to delete previous files if wanted later
    if (.not. MCO%APPEND  .and. MCO%n_WRITE > 0) then
        write(*,*) "Oooops have requested that existing files be deleted but & 
        & not implemented yet."
    end if


    !!!! init params
    ! TODO if not restart / not flag to keep all
    if (.not. restart) then
    call init_pars_random_pregen(PI, PARS_previous, PI%fix_pars, uniform_random_vector)
    endif 
    ! Inform the user
    write(*,*) "Have loaded/randomly assigned PI%parini-now begin the AP-MCMC"
    ! initialize loglikelihood value of the given model with these pars
    !P = -1d0; Pprior = -1d0
    ! calculate the initial probability/log likelihood.
    ! NOTE: passing P0 -> P is needed during the EDC searching phase where we
    ! could read an EDC consistent parameter set in the first instance
    call model_likelihood(PARS_previous, npars, loglikelihood_previous); 
    ! P = P0; Pprior = P0prior
    !write(*,*) "Starting likelihood = ",P0, "+",P0prior
    


    if (is_infinity(loglikelihood_previous)) then  ! TODO is this log evaluating? TODO better check for hard reject ll=inifinity P=0
        write(*,*) "WARNING  ! loglikelihood = ",loglikelihood_previous, " - &
        & AP-MCMC will get stuck, if so please check initial conditions"
        stop
    endif

    ! initalize bestpars to current pars
    BESTPARS = PARS_previous
    llmax = loglikelihood_previous

    endif 


    ! Begin the main AP-MCMC loop
    do while (ITER < MCO%nOUT .and. Pmax < P_target)

       ! take a step in parameter space: generate proposed 
       ! new parameters PARS 
       ! should include reflectivenedd / redrawing 
       call step_pars_real(PARS_previous, PARS_proposed, PI, stats, beta, opt_scaling, par_minstepsize, &
          MCO%N_before_mv_target, uniform_random_vector) 

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
           if (loglikelihood_proposed >= Pmax) then
               BESTPARS = PARS_proposed; Pmax = loglikelihood_proposed
           endif
           ! Keep count of the number of accepted proposals in this local period
           ACCLOC = ACCLOC+1
           ! Accepted first proposal from multivariate
           if (stats%multivariate_proposal) ACC_first = ACC_first+1 !!TODO what is this counter ?
           ! TODO may have to do with whether to calc covariance matrix

           PARS_previous(1:npars) = PARS_proposed(1:npars)          ! save as previous pars
           !norPARS0(1:PI%npars) = norPARS(1:PI%npars)                    ! normalize
           loglikelihood_previous = loglikelihood_proposed;               ! save as previous loglikelihood
        else !TODO
            ! write to history
       endif  ! accept or reject proposed pars

       ! count iteration 
       ITER = ITER+1

       if (MCO%n_WRITE > 0) then
           ! TODO fct
           if (mod(ITER, MCO%n_WRITE) == 0) then
!              print*,"mcmc: write_mcmc_output done"

               ! calculate the likelhood for the actual uncertainties-this avoid
               ! issues with different phases of the MCMC which may use sub-samples
               ! of observations or inflated uncertainties to aid parameter
               ! searching
               call model_likelihood_write(PARS_proposed, npars, output_loglikelihood)
               ! Now write out to files
               call write_mcmc_output(stats%parvar, ACCRATE, &
                                      stats%covariance, &
                                      stats%meanpar, stats%Nparvar, &
                                      PARS_previous, output_loglikelihood, npars, ITER == MCO%nOUT)
           end if 
       end if  ! write or not to write

       ! time to adapt?
       if (mod(ITER, MCO%n_steps) == 0) then
            ! TODO fct
!           ! Debugging print statements
!           print*,"mcmc: time to adapt"
           
           !! update the acceptance counters and acceptance ratios
           ! Total accepted values
           ACC = ACC+ACCLOC
           ! Calculate global acceptance rate
           ACCRATE_GLOBAL = ACC/ITER

           ! Calculate local acceptance rate (i.e. since last adapt)
           ACCRATE = ACCLOC/dble(MCO%n_steps)

           ! Second, are we still in the adaption phase?
           ! TODO how does fortran integer division work
           if (burn_in_period > ITER .or. (ACC_first/ITER) < 0.05d0 .or. .not.MCO%use_multivariate) then

               ! Once covariance matrix has been created just update based on a
               ! single parameter set from each period.
            ! TODO ??
               if (stats%cov) then
                   ACCLOC=1
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
               call update_statistics(PARSALL, npars, stats, stats%use_multivariate, ACCLOC, MCO%N_before_mv_target)

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
               write(*,*)"Overall acceptance rate    = ",ACC/ITER
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
    if (MCO%n_WRITE > 0) call write_covariance_matrix(stats%covariance, npars, .false.)

    ! record the best single set of parameters
    MCOUT%bestpars(1:PI%n_pars) = BESTPARS(1:PI%n_pars)
    ! Output the current state, so that next simulations can potentially resume here
    MCOUT%pars_final(1:PI%n_pars) = PARS_previous(1:PI%n_pars)  !! (!) TODO should be done outside 
    ! record how many iterations were taken to complete
    MCOUT%nos_iterations = MCOUT%nos_iterations+ITER
    ! set flag MCMC completed
    MCOUT%complete = 1
    ! tidy up
    !deallocate(uniform_random_vector)
    ! deallocate this chain's history-but should be fine on subroutine exit

    ! completed AP-MCMC loop
    write(*,*)"AP-MCMC loop completed"
    write(*,*)"Overall acceptance rate     = ", ACC/ITER
    write(*,*)"Final local acceptance rate = ",ACCRATE
    ! TODO function to output these two 
    write(*,*)"Best log-likelihood = ",Pmax
    !write(*,*)"Best parameters = ",MCOUT%best_pars

end subroutine

  !
  !------------------------------------------------------------------
  !
  subroutine update_statistics(PARSALL, npars, stats, use_multivariate, ACCLOC, N_before_mv_target)
    use cardamom_io, only: write_covariance_matrix, write_covariance_info
    use math_functions, only: nor2par, par2nor, log_nor2par, log_par2nor, &
                              cholesky_factor, std, covariance_matrix, &
                              increment_covariance_matrix

    ! Update the multivariate propsal distribution.
    ! Ensure that this subroutine is only called if at least 1 parameter propsal
    ! has been accepted in the last adaption period.

    implicit none

    ! declare input types
    type(MCSTATS) , intent(inout) :: stats ! statistics collection
    logical, intent (inout) :: use_multivariate
    ! declare inputs variables
    integer, intent(in) :: npars
    integer, intent(in) :: ACCLOC
    double precision, intent(in):: PARSALL(npars, ACCLOC)  ! collection of recently accepted normalised parameter combinations
    !! two CARDAMOM quirks here :  1) only accepted steps were recorded, not repeat entries for non-accepted steps
    !!                             2) update_statistics is passed only first or first, middle, and last rows instead of whole history, extremely limnited covariance estimate

    ! declare local variables
    integer p, i, info  ! counters
    double precision, dimension(npars, npars):: cov_backup
    double precision, dimension(npars):: mean_par_backup
    double precision:: Nparvar_backup, Nparvar_local
    double precision :: N_before_mv_target

    ! if we have a covariance matrix then we want to update it, if not then we need to create one
    if (stats%cov) then

        ! Increment the variance-covariance matrix with new accepted parameter sets
        ! NOTE: that this also increments the total accepted counter (PI%Nparvar)

        cov_backup = stats%covariance; mean_par_backup = stats%meanpar; Nparvar_backup =stats%Nparvar

!        call increment_covariance_matrix(PARSALL(1:PI%npars, 1:nint(N%ACCLOC)), PI%mean_par, PI%npars &
!                                        ,PI%Nparvar, nint(N%ACCLOC), PI%covariance)
        ! Have started hardcoding a maximum number of observations to be 100.
        ! While not strictly following Haario et al., (2001) or Roberts and Rosenthal, (2009)
        ! this allows for the covariance matrix to be more responsive to its local environment.
        !! TODO discuss 
        Nparvar_local = min(N_before_mv_target, Nparvar_backup)
        call increment_covariance_matrix(PARSALL(1:npars, 1:ACCLOC), stats%meanpar, npars &
                                        ,Nparvar_local,ACCLOC, stats%covariance)
        ! Calculate the cholesky factor as this includes a determination of
        ! whether the covariance matrix is positive definite.
        call cholesky_factor( npars, stats%covariance, info )
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
                stats%covariance = cov_backup
                stats%meanpar = mean_par_backup
                stats%Nparvar = Nparvar_backup
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
            call covariance_matrix(PARSALL(1:npars, 1:ACCLOC), stats%meanpar, &
            & npars, ACCLOC, stats%covariance)
            stats%cov = .true. ; stats%Nparvar = ACCLOC

            ! Calculate the cholesky factor as this includes a determination of
            ! whether the covariance matrix is positive definite.
            call cholesky_factor ( npars, stats%covariance, info )
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
                call write_covariance_matrix(stats%covariance, npars, .true.)
                call write_covariance_info(stats%meanpar, stats%Nparvar, npars)
            !endif

        end if  ! N%ACCLOC > 2

    end if  ! PI%cov == .true.

    return

  end subroutine update_statistics

  ! Generates new proposed state from currect state in real parameter space.  
  ! Wraps step_pars
  subroutine step_pars_real(PARS0, PARS, PI, stats, beta, opt_scaling, par_minstepsize, N_before_mv_target, random_uniform_vector)
    use math_functions, only: log_par2nor, log_nor2par
    use random_uniform, only : UNIF_VECTOR
    implicit none
    double precision, dimension(:), intent(in   ):: pars0    ! current parameters
    double precision, dimension(:), intent(out  ):: pars       ! proposal
    type(UNIF_VECTOR), intent(inout) :: random_uniform_vector
    type(PARINFO) :: PI
    type(MCSTATS), intent(inout) :: stats
    double precision, dimension(PI%n_pars)             :: pars0_norm, pars_norm
    double precision, intent(in) :: beta, opt_scaling, par_minstepsize
    double precision, intent(in) :: N_before_mv_target
    pars0_norm = log_par2nor(PI%n_pars, pars0, PI%parmin, PI%parmax, PI%paradj)
    call step_pars(pars0_norm, pars_norm, PI%n_pars, stats, beta, opt_scaling, par_minstepsize, &
        N_before_mv_target, random_uniform_vector )
    pars = log_nor2par(PI%n_pars, pars_norm, PI%parmin, PI%parmax, PI%paradj)
  end subroutine 

  !
  !------------------------------------------------------------------
  !
  ! Applies Roberts and Rosenthal 2009 - Eq 3 to generate new proposed state in (normalized)
  ! parameter space
  ! IN : Pars0 current state (normalized)
  ! IN : Stats MCOSTATS object holding info like covariance matrix of the run so far (calculated on 
  ! normalized space )
  ! OUT: PARS new proposed state (normalized)
  ! plus take beta from module data
  subroutine step_pars(PARS0, PARS, npars, stats, beta, opt_scaling, par_minstepsize, N_before_mv_target, random_uniform_vector) !TODO check against original !! 
    use math_functions, only:  random_normal, random_multivariate
    use random_uniform, only : UNIF_VECTOR

    ! carries out the next step to parameters in the MCMC search

    implicit none

    ! declare input variables
    !double precision, dimension(PI%npars), intent(inout):: !norpars0 & ! normalised current parameters
                                                           !,norpars  & ! normalised proposal
    double precision, dimension(:), intent(in   ):: pars0    ! current parameters
    double precision, dimension(:), intent(out  ):: pars       ! proposal
    integer, intent(in) :: npars
    type(UNIF_VECTOR), intent(inout) :: random_uniform_vector
    !type(MHMCMCOPT) , intent(in) :: MCO
    type(MCSTATS), intent(inout) :: stats
    double precision, intent(in) :: beta, opt_scaling, par_minstepsize
    double precision, intent(in) :: N_before_mv_target
    
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
    do p=1, npars
        call random_normal(random_uniform_vector, rn2(p))
    end do

    if ((stats%use_multivariate .and. stats%Nparvar > N_before_mv_target)) then

        ! Is this step a multivariate proposal or not
        stats%multivariate_proposal = .true. ! this only affects ACC_first counter !TODO

        ! Draw from multivariate random distribution
        ! NOTE: if covariance matrix provided is not positive definite
        !       a sample form normal distribution is returned
        call random_multivariate(npars, 1, stats%covariance, mu, rn, random_uniform_vector)

        ! Estimate the step to be applied to the current parameter vector to
        ! create the new proposal. scd = a scaling parameter linking searching
        ! stepping to the number of parameters being retrieved by the analysis.
        ! See Haario et al., (2001) An adaptive Metropolis algorithm. Bernoulli 7.2: 223-242.
        ! and references therein. See also, Roberts & Rosenthal (2009) for beta scaling.
        ! TODO par_minstepsize a user input with a default value 
        pars = pars0 + (rn*opt_scaling * (1d0-beta)) + (par_minstepsize*rn2*beta)

    else 

       stats%multivariate_proposal = .false.
       pars = pars0 + (par_minstepsize*rn2)

    end if

  end subroutine step_pars
  !
  !------------------------------------------------------------------
  !
end module MHMCMC_module
