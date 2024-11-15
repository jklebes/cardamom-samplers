module DEMCz_module

   !!!!!!!!!!!
   ! MC differential evolution z algorithm
   ! jklebes 2024
   ! Implementing ter Braak & Vrugt 2008
   !!!!
   implicit none

   public 

   public :: DEMCz

   type DEMCzOPT
      integer :: MAXITER ! overall steps
      integer :: n_steps !related to time for adapt
      integer :: N_chains ! consider setting OMP env to something compatible
      real :: f = 0.8 !TODO names, default values
      real ::  gamma = 0.8 ! parameters of DE algorithm
      real :: P_target ! termination criteria
   end type DEMCzOPT

   type PARINFO
      integer :: n_pars
      real, allocatable, dimension(:) :: parmin, parmax
   end type PARINFO

   type MCMC_OUTPUT
      real :: bestll
      real, allocatable, dimension(:) :: bestpars
   end type MCMC_OUTPUT

contains
   !
   !--------------------------------------------------------------------
   !
! PI and MCO structs come from model MODEL_LIKELIHOOD files, at the moment
   subroutine DEMCz(model_likelihood_write, model_likelihood, PI, MCO, MCOUT)
    implicit none
      !------------------------------------------------------------------
      !
      ! interface for function
       interface
      function model_likelihood_write(param_vector) result(ML) 
           implicit none
           ! declare input variables
           real, dimension(:), intent(in):: param_vector
           ! output
           real :: ML
      end function model_likelihood_write
        end interface
       interface
      function model_likelihood(param_vector) result(ML)
           implicit none
           ! declare input variablesinteger
           real, dimension(:), intent(in):: param_vector
           ! output
           real :: ML
      end function model_likelihood
        end interface
      

      !! inputs
      type(PARINFO), intent(in) :: PI
      type(DEMCzOPT), intent(in) :: MCO
      type(MCMC_OUTPUT), intent(out) :: MCOUT

      real, allocatable, dimension(:,:) :: PARS_current ! Matrix X , d x N
      ! and their current loglikelihood values
      real, allocatable, dimension(:) :: l0
      ! and their best likelihood values 
      real, allocatable, dimension(:) :: l_best
      real, allocatable, dimension(:,:) :: PARS_best 

      real, allocatable, dimension(:,:) :: PARS_history ! Matrix Z , d x final value of M 
      
      integer :: npars 
      integer :: nchains, nsteps, len_history, MAXITER
      integer :: i,j

      npars = PI%n_pars
      nchains = MCO%N_chains
      nsteps = MCO%N_chains
      MAXITER = MCO%MAXITER
      !TODO check the function takes npars arguments 

      ! TODO put these arrays the right way around 
      allocate(PARS_current(npars, Nchains)) !TODO or each parallel worker could hold its own array, private
      ! but we need these to persist between parallel regions
      allocate(l0(nchains))
      allocate(l_best(nchains))
      allocate(PARS_best(npars,nchains))
      allocate(PARS_history(npars, MCO%MAXITER))

      ! number of entries in history matrix so far
      len_history = 0

!$    OMP PARALLEL DO
      do i=1, nchains
         ! choose initial values
         call init_random(PI, PARS_current(:,i), l0(i))
         ! set pars from nor

         ! potential burnin steps
         ! write first values to history matrix

      end do
!$    OMP END PARALLEL DO
      len_history = len_history+n_chains

      do i = 2, MAX_ITER
!$       OMP PARALLEL DO
         do j=1, n_chains
            call step_chain(PARS_current(i, :), l0(i), PARS_history, model_likelihood, & 
            npars, len_history, nsteps)
         end do
         ! write to Z
!$       OMP END PARALLEL DO !!Barrier implicit ?
      len_history = len_history+n_chains
         ! check convergence
         ! Reorder for best chains ?  Then write to Z later.
      end do

   end subroutine

   subroutine init_random(PI, pars, ll)
      type(PARINFO), intent(in) :: PI
      real, dimension(:), intent(out) :: pars
      real, intent(out) :: ll
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
