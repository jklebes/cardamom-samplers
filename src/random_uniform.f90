module random_uniform
   !! An object-oriented approach to a pre-filled array of random uniform values (0,1)
   !! Can no longer be kept centrally, because each chain needs its own
   !! Checking for need to re-fill and re-filling in a getter seems neater than checking 
   !! and potentially refilling in each place it's used

  public UNIF_VECTOR, get_random_uniform, next_random_uniform
    type UNIF_VECTOR
        integer :: length
        double precision, dimension(:), allocatable :: u
        integer :: index
        contains 
        procedure :: initialize
        procedure :: get_random_uniform 
        procedure :: next_random_uniform 
    end type

    contains

  subroutine initialize(this)
    class(UNIF_VECTOR) :: this
    call fill_random_uniform(this%u, this%length)
    this%index = 1
  end subroutine


  subroutine get_random_uniform(this, n, x)
    !! getter from array of pre-generated random values, 
    !! handling re-filling of the array when needed
    class(UNIF_VECTOR) :: this
    integer, intent(in)  :: n 
        !! number of random values to get
    double precision, dimension(:), allocatable , intent(out) :: x 
        !! array of n values out

    if (this%index + n  > this%length) then ! refill if running out of random values
        call fill_random_uniform(this%u, this%length)
    endif
    ! TODO not handled, will get stuck in infinite loop:  n > this%length   

    x(1:n) = this%u(this%index : this%index + n)

    ! update the index pointer
    this%index = this%index+ n

  end subroutine

  double precision function next_random_uniform(this) result(x)
    !! getter for one (scalar) pre-generated random value, 
    !! handling re-filling of the array when needed
    class(UNIF_VECTOR) :: this
        !! random value out

    if (this%index + 1  > this%length) then ! refill if running out of random values
        call fill_random_uniform(this%u, this%length)
    endif

    x = this%u(this%index)

    ! update the index pointer
    this%index = this%index+ 1

  end function

  subroutine fill_random_uniform(u, n)

    ! Generate an array of n double precision values between 0 and 1.
    ! from Seminumerical Algorithms by D E Knuth, 3rd edition (1997)
    !       including the MODIFICATIONS made in the 9th printing (2002)
    ! https://www-cs-faculty.stanford.edu/~knuth/programs.html#rng
    ! ********* see the book for explanations and caveats! *********
    ! Author: Steve Kifowit
    ! http://ourworld.compuserve.com/homepages/steve_kifowit
    ! with modifications by Alan Miller to rnarry and rnstrt based upon
    ! Knuth's code.
    ! Code converted using TO_F90 by Alan Miller
    ! Date: 2000-09-10, last update 16 January 2003
    ! Modified for integration into CARDAMOM by T. Luke Smallman (t.l.smallman@ed.ac.uk)
    ! 03/05/2019

    integer, intent(in)  :: n  ! number of random values wanted
    double precision, intent(out):: u(n)  ! output vector

    ! Local array
    integer, allocatable, dimension(:)  :: aa

    ! allocate memory
    allocate(aa(n))

    call rnarry(aa, n)
    u(1:n) = scale( dble(aa), -30)

    ! tidy
    deallocate(aa)

    return

  end subroutine fill_random_uniform

end module