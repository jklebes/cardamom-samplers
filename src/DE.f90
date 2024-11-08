module DEMCz_module

  !!!!!!!!!!!
  ! MC differential evolution z algorithm 
  ! jklebes 2024
  ! Implementing ter Braak & Vrugt 2008 
  !!!!
  implicit none

private

public :: DEMCz ! main function

!

contains
  !
  !--------------------------------------------------------------------
  !
  subroutine DEMCz(function)
    
  !------------------------------------------------------------------
  !
  real, allocatable, dimension(:,:) :: PARS_current ! Matrix X , d x N
  real, allocatable, dimension(:,:) :: PARS_history ! Matrix Z , d x final value of M

  allocate(PARS_current(len_pars, N_population))
  allocate(PARS_history(len_pars, N_steps))

  ! potential burnin steps

  ! choose initial values

  do i = 1, MAX_ITER 
    call step_population(PARS_current, PARS_history)
  end do

  end subroutine

  subroutine step_population(PARS_current, PARS_history)
    do i=1, N_population !TODO in case it doesnt divide and we're short
        PARS_history(:, i) = step_vector(PARS_current(:,i), PARS_history)
    end do
  end subroutine 

  subroutine step_vector(vector)
    R1 = random()
    R2 = Random() ! and not equal R1
    do i=1, len_pars
        if ((i == R1) .or. (i==R2)) then
            prop(i) = mutate(prop(i))
        else
            rand = random() ! 0 to 1
            if (rand <= f) then
                prop(i) = mutate(prop(i))
            endif 
        endif 
    end do
    call metropolisChoice(vector, prop_vector) !We should have this in MCMC common
  end subroutine

  subroutine mutate(v1, v2, v3, gamma) !TODO function
    ! update a single value to v1 = v1 + gamma(v2 -v3) + noise
  end subroutine

end module DEMCz_module
