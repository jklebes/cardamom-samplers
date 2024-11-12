module DEMCz_module

   !!!!!!!!!!!
   ! MC differential evolution z algorithm
   ! jklebes 2024
   ! Implementing ter Braak & Vrugt 2008
   !!!!
   implicit none

   public !expose all for testing
   real :: f = 0.8
   real :: gamma = 0.8
   integer :: k = 100 !TODO

!

contains
   !
   !--------------------------------------------------------------------
   !
   subroutine DEMCz(model_likelihood, MAX_ITER)
    use MCMCOPT, only: PI
    implicit none
      !------------------------------------------------------------------
      !
      ! interface for function


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
      ! function e : always gaussian
      !

      real, allocatable, dimension(:,:) :: PARS_current ! Matrix X , d x N
      real, allocatable, dimension(:,:) :: PARS_history ! Matrix Z , d x final value of M

      integer, intent(in) :: MAX_ITER
      integer :: npars 
      integer :: N_chains = 12 !TOTO
      integer :: N_Steps = 10000 ! TODO
      integer :: len_history = 12 !TODO

      real :: l0 = 0.0 ! TODO initalize parameters, likelihoods

      integer :: i,j
      npars = PI%npars
      ! TODO put these arrays the right way around 
      allocate(PARS_current(npars, N_chains)) !TODO or each parallel worker could hold its own array, private
      ! but we need these to persist between parallel regions

      allocate(PARS_history(npars, N_steps))

!$    OMP PARALLEL DO
      do i=1, n_chains
         ! choose initial values
         ! potential burnin steps
      end do
!$    OMP END PARALLEL DO


      do i = 1, MAX_ITER
!$       OMP PARALLEL DO
         do j=1, n_chains
            call step_chain(PARS_current(i, :), l0, PARS_history, model_likelihood, npars, len_history)
         end do
         ! write to Z
!$       OMP END PARALLEL DO !!Barrier implicit    
         ! check convergence
         ! Reorder for best chains ?  Then write to Z later.
      end do

   end subroutine

   subroutine step_chain(X_i, l0, PARS_history, model_likelihood, npars, len_history)
    real, dimension(:), intent(inout) :: X_i
    real, dimension(:), allocatable :: vector , prop_vector !internal: save previous state, proposed new state
    real, dimension(:,:), intent(in) :: PARS_history !the matrix Z so far, to read 2 rows from 
    real, intent(inout) :: l0 ! likelihood of previous accepted params
    real             :: l !likelihood of proposed values
    integer, intent(in) :: npars, len_history 
    integer :: R1, R2 !indices of 2 random rows (may be same)
    real :: rand
    integer :: i

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
      R2 = random_int(len_history) ! and not equal R1
      vector = X_i
      call mutate(prop_vector, vector, PARS_history(R1,:), PARS_history(R2,:))
      call model_likelihood(prop_vector, l)
      if (metropolis_Choice(l, l0)) then !We should have this in MCMC common
         X_i = prop_vector
         l0=l
      end if 
      end do
   end subroutine

   subroutine mutate(vout, v1, v2, v3) 
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
