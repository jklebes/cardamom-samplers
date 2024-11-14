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

  ! refactored jklebes 2024

  ! Module contains all subroutine and functions relevant specifically to the
  ! AP-MCMC method. The choice of EDC, likelihood and model are made else where and
  ! are thus contains within a seperate module

  ! Relevant source references:
  ! Haario et al., (2001) An adaptive Metropolis algorithm. Bernoulli 7.2: 223-242.
  ! Haario et al., (2006) Stat. Comput., 16:339–354, DOI 10.1007/s11222-006-9438-0, 
  ! Roberts and Rosenthal (2009), Examples of Adaptive MCMC, J. Comp. Graph. Stat. 18:349-367

implicit none

! make all private
private

! specify what can be seen
public:: MHMCMC

! declare any module level variables needed
! Settings for MCMC run - could be common - nchains, max len
! MCOPT struct with default values
integer, parameter :: nchains =1
double precision, parameter:: N_before_mv = 10d0
! min, max pars

!! Adaptive !!
! setting for adaptive AP-MCMC step size
double precision, parameter:: par_minstepsize = 0.001d0 & ! 0.0005 -> 0.001 -> 0.01 -> 0.1 -> 0.005
                              ,par_maxstepsize = 0.01d0  &
                              ,par_initstepsize = 0.005d0

! Optimal scaling variable for parameter searching
double precision:: N_before_mv_target, & !
                    opt_scaling  ! scd = 2.381204 the optimal scaling parameter
                                ! for MCMC search, when applied to  multivariate proposal.
                                ! NOTE 1: 2.38/sqrt(npars) sometimes used when applied to the Cholesky
                                ! factor. NOTE 2: 2.381204**2 = 5.670132

!! step
double precision, parameter:: beta = 0.05d0  ! weighting for gaussian step in multivariate proposals
! Is current proposal multivariate or not?
logical:: multivariate_proposal = .false.



contains
  !
  !--------------------------------------------------------------------
  !
  subroutine MHMCMC( model_likelihood_write, model_likelihood, MCOUT)
    use cardamom_io, only: write_parameters, write_variances, write_covariance_matrix &
                          ,write_covariance_info, restart_flag, write_mcmc_output

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
    ! * These will be set to default values if empty. Options include:
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

    !type(adaptmcopt), intent(in):: MCOPT ! collection of settings for the adaptive MCMC run
    !type(parinfo) , intent(in):: PI      ! collection of parameter info incl number, type, bounds
    !type(mcout), intent(out):: MCOUT     ! struct for output: best loglikelihood, best parameters !TODO from multiple chains?  get array of them.

    ! declare interface for the model likelihood function.
    ! A function pars vector -> loglikelihood
    ! Do not modify here, sampler is generic.
    ! Wrap the model loglikelihood function in another function elsewhere to 
    ! make it conform to this form.
    interface
      subroutine model_likelihood(param_vector, ML_obs_out, ML_prior_out)
        use MCMCOPT, only: PI
        use CARBON_MODEL_MOD, only: carbon_model
           implicit none
           ! declare input variables
           double precision, dimension(PI%npars), intent(in):: param_vector
           ! output
           double precision, intent(out):: ML_obs_out, ML_prior_out
      end subroutine model_likelihood
    end interface

    ! an alternative function, usually strictly a log-likelihood where the actual 
    ! calculation is performed with some variant, 
    ! only used for writing the true log-likelihood to output file
    interface
      subroutine model_likelihood_write(param_vector, ML_obs_out, ML_prior_out)
        use MCMCOPT, only: PI
        use CARBON_MODEL_MOD, only: carbon_model
           implicit none
           ! declare input variables
           double precision, dimension(PI%npars), intent(in):: param_vector
           ! output
           double precision, intent(out):: ML_obs_out, ML_prior_out
      end subroutine model_likelihood_write
    end interface

    !! Calculate shared settings of the run, following from input options.  
    ! Keep them at module data or in a struct ?
    ! Determine how long we will continue to adapt our proposal covariance
    ! matrix and use of Delayed Rejection
    burn_in_period = MCO%fADAPT*dble(MCO%nOUT)
    N_before_mv_target = N_before_mv*dble(PI%npars)


    ! This adaptive MCMC is "parallelized" to do N chains at the same time, but they
    ! each do a complete independent run.  The parallelization structure of this module is 
    ! trivial.  It exists mainly as template for samplers with more crossover and more complex
    ! structure in this loop.
    !$OMP parallel do
    do i=1, nchains
        ! saves best loglikelihood and associated parameters to MCOUT(i)
        run_mcmc(model_likelihood_write, model_likelihood, PI, MCO, MCOUT(i))
    end do
    !$OMP end parallel do

    !! compare chains

  end subroutine MHMCMC

  subroutine run_mcmc(model_likelihood_write, model_likelihood, PI, MCO, MCOUT)
    use math_functions, only: random_uniform, log_par2nor, log_nor2par, par2nor, nor2par
    ! declare any local variables
    !! all variables in here are local to the single chain and the duraction of its run
    type ( counters ):: N
    double precision, dimension(PI%npars):: norPARS0      & ! normalised parameter values for current state
                                            ,norPARS       & ! normalised parameter values for current proposal
                                            ,PARS0         & ! parameter values for current state
                                            ,PARS          & ! parameter values for current proposal
                                            ,BESTPARS        ! best set of parameters so far

    double precision, dimension(PI%npars, MCO%nADAPT):: PARSALL  ! All accepted normalised parameters since previous step adaption
    double precision:: burn_in_period & ! for how many proposals will we adapt the covariance matrix as a minimum
                       , AM_likelihood &
                       ,outputP0      &
                       ,outputP0prior  &
                       ,Pmax, P0prior, Pprior & ! as below but for priors only
                       ,P0 & ! previously accepted observation based log-likelihood
                       ,P    ! current observation based log-likelihood
    integer:: i
    ! declare interface for the model likelihood function.
    ! A function pars vector -> loglikelihood
    ! Do not modify here, sampler is generic.
    ! Wrap the model loglikelihood function in another function elsewhere to 
    ! make it conform to this form.
    interface
      subroutine model_likelihood(param_vector, ML_obs_out, ML_prior_out)
        use MCMCOPT, only: PI
        use CARBON_MODEL_MOD, only: carbon_model
           implicit none
           ! declare input variables
           double precision, dimension(PI%npars), intent(in):: param_vector
           ! output
           double precision, intent(out):: ML_obs_out, ML_prior_out
      end subroutine model_likelihood
    end interface

    ! an alternative function, usually strictly a log-likelihood where the actual 
    ! calculation is performed with some variant, 
    ! only used for writing the true log-likelihood to output file
    interface
      subroutine model_likelihood_write(param_vector, ML_obs_out, ML_prior_out)
        use MCMCOPT, only: PI
        use CARBON_MODEL_MOD, only: carbon_model
           implicit none
           ! declare input variables
           double precision, dimension(PI%npars), intent(in):: param_vector
           ! output
           double precision, intent(out):: ML_obs_out, ML_prior_out
      end subroutine model_likelihood_write
    end interface

    ! Initialize basic MCMC counters
    N%ITER = 0d0

    ! Initialize further counters for adaptive
    ! TODO should some of these be integer
    N%ACC = 0d0; N%ACC_first = 0d0;
    N%ACCLOC = 0d0; N%ACCRATE = 0d0; N%ACCRATE_GLOBAL = 0d0

    !!! prepare file writing
    !call prepare_file()
    ! add something here to delete previous files if wanted later
    if (MCO%APPEND == 0 .and. MCO%nWRITE > 0) then
        write(*,*) "Oooops have requested that existing files be deleted but & 
        & not implemented yet."
    end if

    ! allocate the history matrix PARSALL

    !!!! init params
    ! start random sampling if MCO%randparini set
    ! TODO function
    call init_pars()
    

    ! Begin the main AP-MCMC loop
    ! do for i in nchains
    do while (N%ITER < MCO%nOUT .and. Pmax < P_target)

       ! take a step in parameter space - generate proposed 
       ! new parameters PARS and norPARS  
       call step(N, PARS0, PARS, norPARS0, norPARS) 

       ! if parameter proposal in bounds check the model
       ! TODO if we get a reflective gaussian kernel, no need to check
       if (minval(norPARS) > 0d0 .and. maxval(norPARS) < 1d0) then
           ! calculate the model likelihood
           ! TODO ideally output just one likelihood value
           call model_likelihood(PARS, P, Pprior)
           accept = metropolis_decision( (P+Pprior) , (P0+P0prior) ) !make sure it's OMP loop private
       else
           accept = .false.
       end if  ! in bound

       if (accept) then

           ! Store accepted parameter proposals
           ! keep record of all parameters accepted since step adaption
           ! (this chain)
           PARSALL(1:PI%npars, (nint(N%ACCLOC)+1)) = norPARS(1:PI%npars) !add row in history matrix
           ! store the best parameter set
           if ((P+Pprior) >= Pmax) then
               BESTPARS = PARS; Pmax = P+Pprior
           endif
           ! Keep count of the number of accepted proposals in this local period
           N%ACCLOC = N%ACCLOC+1d0 !TODO int
           ! Accepted first proposal from multivariate
           if (multivariate_proposal) N%ACC_first = N%ACC_first+1d0 !!TODO what is this counter ?

           PARS0(1:PI%npars) = PARS(1:PI%npars)                          ! save as previous pars
           norPARS0(1:PI%npars) = norPARS(1:PI%npars)                    ! normalize
           P0 = P; P0prior = Pprior                                      ! save as previous loglikelihood

       endif  ! accept or reject condition

       ! count iteration 
       N%ITER = N%ITER+1

       if (MCO%nWRITE > 0) then
           ! TODO fct
           if (mod(nint(N%ITER), MCO%nWRITE) == 0) then
!              print*,"mcmc: write_mcmc_output done"

               ! calculate the likelhood for the actual uncertainties-this avoid
               ! issues with different phases of the MCMC which may use sub-samples
               ! of observations or inflated uncertainties to aid parameter
               ! searching
               call model_likelihood_write(PARS0, outputP0, outputP0prior)
               ! Now write out to files
               call write_mcmc_output(PI%parvar, N%ACCRATE, &
                                      PI%covariance, &
                                      PI%mean_par, PI%Nparvar, &
                                      PARS0, (outputP0+outputP0prior), PI%npars, N%ITER == MCO%nOUT)
           end if 
       end if  ! write or not to write

       ! time to adapt?
       if (mod(nint(N%ITER), MCO%nADAPT) == 0) then
            ! TODO fct
!           ! Debugging print statements
!           print*,"mcmc: time to adapt"

           ! Total accepted values
           N%ACC = N%ACC+N%ACCLOC
           ! Calculate global acceptance rate
           N%ACCRATE_GLOBAL = N%ACC/N%ITER

           ! Calculate local acceptance rate (i.e. since last adapt)
           N%ACCRATE = N%ACCLOC/dble(MCO%nADAPT)

           ! Second, are we still in the adaption phase?
           if (burn_in_period > N%ITER .or. (N%ACC_first/N%ITER) < 0.05d0 .or. .not.PI%use_multivariate) then

               ! Once covariance matrix has been created just update based on a
               ! single parameter set from each period.
               if (PI%cov) then
                   N%ACCLOC = 1d0; PARSALL(1:PI%npars, nint(N%ACCLOC)) = norPARS0(1:PI%npars)
               else if (.not.PI%cov .and. N%ACCLOC > 3d0) then
                   PARSALL(1:PI%npars, 2) = PARSALL(1:PI%npars, ceiling(N%ACCLOC*0.5d0))
                   PARSALL(1:PI%npars, 3) = PARSALL(1:PI%npars, nint(N%ACCLOC))
                   N%ACCLOC = 3d0
               endif

               ! adapt the covariance matrix for multivariate proposal
               call adapt_step_size(PARSALL, N)

           end if !  have enough parameter been accepted

           ! resets to local counter
           N%ACCLOC = 0d0

       end if  ! time to adapt?

       ! Should I be write(*,*)ing to screen or not?
       if (MCO%nPRINT > 0) then
           if (mod(nint(N%ITER), MCO%nPRINT) == 0) then
               write(*,*)"Using multivariate sampling = ",PI%use_multivariate
               write(*,*)"Total proposal = ",N%ITER, " out of ",MCO%nOUT
               write(*,*)"Total accepted = ",N%ACC
               write(*,*)"Overall acceptance rate    = ",N%ACC/N%ITER
               write(*,*)"Local   acceptance rate    = ",N%ACCRATE
               write(*,*)"Current obs   = ",P0, "proposed = ",P, " log-likelihood"
               write(*,*)"Current prior = ",P0prior, "proposed = ",Pprior, " log-likelihood"
               write(*,*)"Maximum likelihood = ",Pmax
               ! NOTE: that-infinity in current obs only indicates failure of EDCs
               ! but-infinity in both obs and parameter likelihood scores indicates
               ! that proposed parameters are out of bounds
           end if 
       end if  ! write(*,*) to screen or not

    end do  ! while conditions

    ! write out final covariance matrix for the analysis
    if (MCO%nWRITE > 0) call write_covariance_matrix(PI%covariance, PI%npars, .false.)

    ! record the best single set of parameters
    MCOUT%best_pars(1:PI%npars) = BESTPARS(1:PI%npars)
    ! set the initial parameter set the final one accepted
    PI%parini(1:PI%npars) = PARS0(1:PI%npars)  !! (!) TODO should be done outside 
    ! record how many iterations were taken to complete
    MCOUT%nos_iterations = MCOUT%nos_iterations+N%ITER
    ! set flag MCMC completed
    MCOUT%complete = 1
    ! tidy up
    deallocate(uniform_random_vector)
    ! deallocate this chain's history - but should be fine on subroutine exit

    ! completed AP-MCMC loop
    write(*,*)"AP-MCMC loop completed"
    write(*,*)"Overall acceptance rate     = ",N%ACC/N%ITER
    write(*,*)"Final local acceptance rate = ",N%ACCRATE
    ! TODO function to output these two 
    write(*,*)"Best log-likelihood = ",Pmax
    !write(*,*)"Best parameters = ",MCOUT%best_pars
end subroutine()

subroutine init_pars()
    ! initial values
    P = -1d0; Pprior = -1d0

    PI%parfix = 0d0
    do i = 1, PI%npars
       ! parfix = 1 stay at prior value, parfix = 0 randomly search
       if (MCO%fixedpars .and. PI%parini(i) /= -9999d0) PI%parfix(i) = 1d0
       ! only assign random parameters if (a) randparini == .true. or (b) PI$parini(n) == -9999)
       if (MCO%randparini .and. PI%parfix(i) == 0d0 .and. .not.restart_flag) then
!           call nor2par(1, uniform_random_vector(uniform), PI%parmin(i), PI%parmax(i), PI%parini(i))
           call log_nor2par(1, uniform_random_vector(uniform), PI%parmin(i), PI%parmax(i), PI%paradj(i), PI%parini(i))
           uniform = uniform+1
       end if
       ! write(*,*) parameter values to screen
       if (MCO%print_pars) write(*,*) "p",i, "=",PI%parini(i)
    end do  ! for PI%npar loop

    ! Inform the user
    write(*,*) "Have loaded/randomly assigned PI%parini-now begin the AP-MCMC"

    ! Initialise the initial and best pars vectors
    PARS0(1:PI%npars) = PI%parini(1:PI%npars)
    BESTPARS(1:PI%npars) = PI%parini(1:PI%npars)

    ! Track normalised initial proposal
    do i = 1, PI%npars
!       call par2nor(1, PARS0(i), PI%parmin(i), PI%parmax(i), norPARS0(i))
       call log_par2nor(1, PARS0(i), PI%parmin(i), PI%parmax(i), PI%paradj(i), norPARS0(i))
    end do

    ! calculate the initial probability/log likelihood.
    ! NOTE: passing P0 -> P is needed during the EDC searching phase where we
    ! could read an EDC consistent parameter set in the first instance
    call model_likelihood(PI%parini, P0, P0prior); P = P0; Pprior = P0prior

    write(*,*) "Starting likelihood = ",P0, "+",P0prior
    Pmax = P0+P0prior

    ! checks whether the EDCs (combined with P0 not P0prior) have been met in the initial parameter set
    infini = 0d0
    if (P0 == log(infini)) then !TODO is this log evaluating?
        write(*,*) "WARNING  ! P0 = ",P0, " - AP-MCMC will get stuck, if so please check initial conditions"
        stop
    endif
  end subroutine 


  !
  !------------------------------------------------------------------
  !
  subroutine update_covariance(PARSALL, N)
    use cardamom_io, only: write_covariance_matrix, write_covariance_info
    use MCMCOPT, only: MCO, PI, COUNTERS
    use math_functions, only: nor2par, par2nor, log_nor2par, log_par2nor, &
                              cholesky_factor, std, covariance_matrix, &
                              increment_covariance_matrix

    ! Update the multivariate propsal distribution.
    ! Ensure that this subroutine is only called if at least 1 parameter propsal
    ! has been accepted in the last adaption period.

    implicit none

    ! declare input types
    type ( counters ), intent(inout):: N


    ! declare inputs variables
    double precision, intent(in):: PARSALL(PI%npars, MCO%nADAPT)  ! collection of recently accepted normalised parameter combinations

    ! declare local variables
    integer p, i, info  ! counters
    double precision, dimension(PI%npars, PI%npars):: cholesky, cov_backup
    double precision, dimension(PI%npars):: mean_par_backup
    double precision:: Nparvar_backup, Nparvar_local

    ! if we have a covariance matrix then we want to update it, if not then we need to create one
    if (PI%cov) then

        ! Increment the variance-covariance matrix with new accepted parameter sets
        ! NOTE: that this also increments the total accepted counter (PI%Nparvar)
        cov_backup = PI%covariance; mean_par_backup = PI%mean_par; Nparvar_backup = PI%Nparvar
!        call increment_covariance_matrix(PARSALL(1:PI%npars, 1:nint(N%ACCLOC)), PI%mean_par, PI%npars &
!                                        ,PI%Nparvar, nint(N%ACCLOC), PI%covariance)
        ! Have started hardcoding a maximum number of observations to be 100.
        ! While not strictly following Haario et al., (2001) or Roberts and Rosenthal, (2009)
        ! this allows for the covariance matrix to be more responsive to its local environment.
        !! TODO discuss 
        Nparvar_local = min(N_before_mv_target, Nparvar_backup)
        call increment_covariance_matrix(PARSALL(1:PI%npars, 1:nint(N%ACCLOC)), PI%mean_par, PI%npars &
                                        ,Nparvar_local, nint(N%ACCLOC), PI%covariance)
        ! Calculate the cholesky factor as this includes a determination of
        ! whether the covariance matrix is positive definite.
        cholesky = PI%covariance
        call cholesky_factor( PI%npars, cholesky, info )
        ! If the updated covariance matrix is not positive definite we should
        ! reject the update in favour of the existing matrix
        ! TODO ??
        if (info == 0) then
            ! Set multivariate sampling to true
            PI%use_multivariate = .true.
        else
            ! The current addition of a parameter leads to a matrix which is not
            ! positive definite. If we previously had a matrix which is positive
            ! definite then we should reject totally the new matrix, if not we
            ! should keep it and accumulate the information
            if (PI%use_multivariate) then
                ! return original matrix to place
                PI%covariance = cov_backup
                PI%mean_par = mean_par_backup
                PI%Nparvar = Nparvar_backup
            else
                ! Keep accumulating the information
                PI%use_multivariate = .false.
            end if
        endif

    else  ! PI%cov == .true.

        ! we have not yet created a covariance matrix based on accepted
        ! parameters. Assuming we have some then create one...
        if (N%ACCLOC > 2d0) then

            ! estimate covariance matrix
            call covariance_matrix(PARSALL(1:PI%npars, 1:nint(N%ACCLOC)), PI%mean_par, PI%npars, nint(N%ACCLOC), PI%covariance)
            PI%cov = .true. ; PI%Nparvar = N%ACCLOC

            ! Calculate the cholesky factor as this includes a determination of
            ! whether the covariance matrix is positive definite.
            cholesky = PI%covariance
            call cholesky_factor ( PI%npars, cholesky, info )
            ! If not positive definite then we should not use the multivariat
            
            ! step at this time.
            if (info /= 0) then
                ! Keep accumulating information until positive definite matrix
                ! calculated
                PI%use_multivariate = .false.
            else
                ! Positive definite found straight away-might as well use it!
                PI%use_multivariate = .true.
            endif

            ! write out first covariance matrix, this will be compared with the final covariance matrix
            if (MCO%nWRITE > 0) then
                call write_covariance_matrix(PI%covariance, PI%npars, .true.)
                call write_covariance_info(PI%mean_par, PI%Nparvar, PI%npars)
            endif

        end if  ! N%ACCLOC > 2

    end if  ! PI%cov == .true.

    ! adjust step size by local variance
    !do p = 1, PI%npars
       ! calculate standards deviation (variability) for the local
       ! window extracted from the variance
       ! TODO used? Only for writing, which is not needed as cov is also written
       !PI%parvar(p) = PI%covariance(p, p)
    !end do  ! p

    return

  end subroutine update_covariance
  !
  !------------------------------------------------------------------
  !
  subroutine step(N, pars0, pars, norpars0, norpars)
    use math_functions, only: par2nor, nor2par, log_par2nor, log_nor2par, &
                              random_normal, random_multivariate, random_uniform
    use MCMCOPT, only: PI, COUNTERS

    ! carries out the next step to parameters in the MCMC search

    implicit none

    ! declare input variables
    type ( counters ), intent(inout):: N
    double precision, dimension(PI%npars), intent(inout):: norpars0 & ! normalised current parameters
                                                           ,norpars  & ! normalised proposal
                                                           ,pars0    & ! current parameters
                                                           ,pars       ! proposal

    ! declare local variables
    integer:: p
    double precision:: rn(PI%npars), mu(PI%npars), rn2(PI%npars)

    ! mean of distributions 
    mu = 0d0 

    ! Splitting step calculation based on number of parameter vectors accepted
    ! is linked to the need build a covariance matrix prior to multivariate
    ! sampling.
    ! See Roberts and Rosenthal, Examples of Adaptive MCMC, J. Comp. Graph. Stat. 18:349-367, 2009.
    
    ! Sample random normal distribution (mean = 0, sd = 1)
       ! TODO vectorize
       do p = 1, PI%npars
          call random_normal(rn2(p))
       end do !

    if ((PI%use_multivariate .and. PI%Nparvar > N_before_mv_target)) then

        ! Is this step a multivariate proposal or not
        ! multivariate_proposal = .true.

        ! Draw from multivariate random distribution
        ! NOTE: if covariance matrix provided is not positive definite
        !       a sample form normal distribution is returned
        call random_multivariate(PI%npars, 1, PI%covariance, mu, rn)

        ! Estimate the step to be applied to the current parameter vector to
        ! create the new proposal. scd = a scaling parameter linking searching
        ! stepping to the number of parameters being retrieved by the analysis.
        ! See Haario et al., (2001) An adaptive Metropolis algorithm. Bernoulli 7.2: 223-242.
        ! and references therein. See also, Roberts & Rosenthal (2009) for beta scaling.
        norpars = norpars0 + (rn*opt_scaling * (1d0-beta)) + (par_minstepsize*rn2*beta)

    else  ! nint(PI%Nparvar) > N_before_mv*PI%npars

       ! multivariate_proposal = .false.
       norpars = norpars0 + (par_minstepsize*rn2)

    end if  ! nint(PI%Nparvar) > 10*PI%npars

    ! reverse normalisation on the new parameter step
    ! TODO not here ?
    do p = 1, PI%npars
!       call nor2par(1, norpars(p), PI%parmin(p), PI%parmax(p), pars(p))
       call log_nor2par(1, norpars(p), PI%parmin(p), PI%parmax(p), PI%paradj(p), pars(p))
    end do

  end subroutine step
  !
  !------------------------------------------------------------------
  !
end module MHMCMC_module
