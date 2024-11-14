module DEMCz_module

   !!!!!!!!!!!
   ! MC differential evolution z algorithm
   ! jklebes 2024
   ! Implementing ter Braak & Vrugt 2008
   !!!!
   implicit none

   private 
   ! Hold settings of the DEMCz run here, they are common to all
   ! (possibly OMP parallel) chains and are set in DEMCz from 
   ! TODO optional arguments or defaults
   ! TODO package in DEMCzOPT struct
   real :: f = 0.8
   real :: gamma = 0.8
   integer :: k = 100 !TODO
   integer :: N_chains = 12 !TODO
   real :: burn_in_period
   ! termination criteria:
   integer :: MAX_ITER = 500 
   real :: P_target 

   public :: DEMCz
!

contains
   !
   !--------------------------------------------------------------------
   !
   subroutine DEMCz(model_likelihood_write, model_likelihood, MCOUT)
    use MCMCOPT, only: MCOUT, PI !we can assume this has been set to hold number, bounds of pars
    ! Not the ideal way to pass settings in and results out 
    implicit none
      !------------------------------------------------------------------
      !
      ! interface for function
       interface
      subroutine model_likelihood_write(param_vector, ML)
        use MCMCOPT, only: PI
           implicit none
           ! declare input variables
           real, dimension(PI%npars), intent(in):: param_vector
           ! output
           real, intent(out):: ML
      end subroutine model_likelihood_write
       interface
      subroutine model_likelihood(param_vector, ML)
        use MCMCOPT, only: PI
           implicit none
           ! declare input variables
           real, dimension(PI%npars), intent(in):: param_vector
           ! output
           real, intent(out):: ML
      end subroutine model_likelihood
        end interface
      

      real, allocatable, dimension(:,:) :: PARS_current ! Matrix X , d x N
      real, allocatable, dimension(:,:) :: PARS_history ! Matrix Z , d x final value of M 
      
      integer :: npars 

      ! TODO math constants
      real :: loglikelihood0 = infini 

      integer :: i,j

      npars = PI%npars
      !TODO check the function takes npars arguments 

      ! TODO put these arrays the right way around 
      allocate(PARS_current(npars, N_chains)) !TODO or each parallel worker could hold its own array, private
      ! but we need these to persist between parallel regions

      allocate(PARS_history(npars, N_steps))

!$    OMP PARALLEL DO
      do i=1, n_chains
         ! choose initial values
         call init_random(PARS_current(:,i))
         ! set pars from nor

         ! potential burnin steps
      end do
!$    OMP END PARALLEL DO


      do i = 1, MAX_ITER
!$       OMP PARALLEL DO
         do j=1, n_chains
            call step_chain(PARS_current(i, :), l0, PARS_history, model_likelihood, npars, len_history)
         end do
         ! write to Z
!$       OMP END PARALLEL DO !!Barrier implicit ?
         ! check convergence
         ! Reorder for best chains ?  Then write to Z later.
      end do

   end subroutine

   !> Evolve the state of one chain by k steps
   subroutine step_chain(X_i, l0, PARS_history, model_likelihood, npars, len_history)
      USE DEMCzOPT, only: k
    real, dimension(:), intent(inout) :: X_i
    real, dimension(:), allocatable :: vector , prop_vector !internal: save previous state, proposed new state
    real, dimension(:,:), intent(in) :: PARS_history !the matrix Z so far, to read 2 rows from 
    real, intent(inout) :: l0 ! likelihood of previous accepted params
    real             :: l !likelihood of proposed values
    integer, intent(in) :: npars, len_history 
    integer :: R1, R2 !indices of 2 random rows (may be same)
    real :: rand
    integer :: i, k

    interface
    subroutine model_likelihood(param_vector, ML)
      use MCMCOPT, only: PI
         implicit none
         real, dimension(PI%npars), intent(in):: param_vector !TODO kind double precision ?
         real, intent(out):: ML
    end subroutine model_likelihood
      end interface

      do i=1, k 
      R1 = random_int(len_history)
      R2 = random_int(len_history) ! without replacement - possibly equal R1
      vector = X_i
      call step(prop_vector, vector, PARS_history(R1,:), PARS_history(R2,:))
      call model_likelihood(prop_vector, l)
      if (metropolis_Choice(l, l0)) then !We should have this in MCMC common
         X_i = prop_vector
         l0=l
      end if 
      end do
   end subroutine

   subroutine step(vout, v1, v2, v3) 
    real, dimension(:), intent(out) :: vout
    real, dimension(:), intent(in) :: v1, v2, v3
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
