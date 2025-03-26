module random_uniform
   !! An object-oriented approach to a pre-filled array of random uniform values (0, 1)
   !! Can no longer be kept centrally, because each chain needs its own
   !! Checking for need to re-fill and re-filling in a getter seems neater than checking 
   !! and potentially refilling in each place it's used
  
   ! rand(), narray(), rnstrt() related
  integer, parameter  :: kk = 100, ll = 37, mm = 2**30, tt = 70, kkk = kk+kk-1
  integer, save       :: ranx(kk)


  public UNIF_VECTOR, get_random_uniform, next_random_uniform
    type UNIF_VECTOR
        integer :: length = 1000 !TODO get default from orig cardamom
        double precision, dimension(:), allocatable :: u
        integer:: index
        contains 
        procedure:: initialize
        procedure:: get_random_uniform
        procedure:: next_random_uniform
    end type

    contains

  subroutine initialize(this)
    class(UNIF_VECTOR):: this
    integer :: seed = 100
    call rnstrt(seed) ! TODO no longer necessary to do from outside take out of main
    ! TODO pass seed down
    if (.not.allocated(this%u)) then
      allocate(this%u(this%length))
    end if
    call fill_random_uniform(this%u, this%length)
    this%index = 1
  end subroutine


  function get_random_uniform(this, n) result(x)
    !! getter from array of pre-generated random values, 
    !! handling re-filling of the array when needed
    class(UNIF_VECTOR):: this
    integer, intent(in)  :: n 
        !! number of random values to get
    double precision, dimension(:), allocatable:: x 
        !! array of n values out

    if (this%index+n  > this%length) then  ! refill if running out of random values
        call fill_random_uniform(this%u, this%length)
        this%index = 1
    endif
    ! TODO not handled, will get stuck in infinite loop:  n > this%length   

    x  = this%u(this%index : this%index+n)

    ! update the index pointer
    this%index = this%index+n

  end function

  double precision function next_random_uniform(this) result(x)
    !! getter for one (scalar) pre-generated random value, 
    !! handling re-filling of the array when needed
    class(UNIF_VECTOR):: this
        !! random value out

    if (this%index+1  > this%length) then  ! refill if running out of random values
        call fill_random_uniform(this%u, this%length)
        this%index = 1
    endif
    x = this%u(this%index)

    ! update the index pointer
    this%index = this%index+1

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
  !
  !--------------------------------------------------------------------
  !
  subroutine rnarry(aa, n)

    ! Generate an array of n integers between 0 and 2^30-1.
    ! Part of process to an array of n double precision values between 0 and 1.
    ! from Seminumerical Algorithms by D E Knuth, 3rd edition (1997)
    !       including the MODIFICATIONS made in the 9th printing (2002)
    ! ********* see the book for explanations and caveats! *********
    ! Author: Steve Kifowit
    ! http://ourworld.compuserve.com/homepages/steve_kifowit
    ! with modifications by Alan Miller to rnarry and rnstrt based upon
    ! Knuth's code.
    ! Code converted using TO_F90 by Alan Miller
    ! Date: 2000-09-10, last update 16 January 2003
    ! Modified for integration into CARDAMOM by T. Luke Smallman (t.l.smallman@ed.ac.uk)
    ! 03/05/2019

    integer, intent(in)   :: n
    integer, intent(out)  :: aa(n)

    ! Local variables
    integer  :: j

    aa(1:kk) = ranx(1:kk)
    do j = kk+1, n
       aa(j) = aa(j-kk) - aa(j-ll)
       if (aa(j) < 0) aa(j) = aa(j) + mm
    end do
    do j = 1, ll
       ranx(j) = aa(n+j-kk) - aa(n+j-ll)
       if (ranx(j) < 0) ranx(j) = ranx(j) + mm
    end do
    do j = ll+1, kk
       ranx(j) = aa(n+j-kk) - ranx(j-ll)
       if (ranx(j) < 0) ranx(j) = ranx(j) + mm
    end do

    return

  end subroutine rnarry
  !
  !--------------------------------------------------------------------
  !
  subroutine rnstrt(seed)

    ! Initialize integer array ranx using the input seed.
    ! Part of process to an array of n double precision values between 0 and 1.
    ! from Seminumerical Algorithms by D E Knuth, 3rd edition (1997)
    !       including the MODIFICATIONS made in the 9th printing (2002)
    ! ********* see the book for explanations and caveats! *********
    ! Author: Steve Kifowit
    ! http://ourworld.compuserve.com/homepages/steve_kifowit
    ! with modifications by Alan Miller to rnarry and rnstrt based upon
    ! Knuth's code.
    ! Code converted using TO_F90 by Alan Miller
    ! Date: 2000-09-10, last update 16 January 2003
    ! Modified for integration into CARDAMOM by T. Luke Smallman (t.l.smallman@ed.ac.uk)
    ! 03/05/2019

    integer, intent(in)  :: seed

    ! Local variables
    integer  :: x(kkk), j, ss, sseed, t

    if (seed < 0) then
        sseed = mm-1 - mod(-1-seed, mm)
    else
        sseed = mod(seed, mm)
    end if
    ss = sseed-mod(sseed, 2) + 2
    do j = 1, kk
       x(j) = ss
       ss = ishft(ss, 1)
       if (ss >= mm) ss = ss-mm+2
    end do
    x(kk+1:kkk) = 0
    x(2) = x(2)+1
    ss = sseed
    t = tt-1
10  do j = kk, 2, -1
       x(j+j-1) = x(j)
    end do
    do j = kkk, kk+1, -1
       x(j-(kk-ll)) = x(j-(kk-ll)) - x(j)
       if (x(j-(kk-ll)) < 0) x(j-(kk-ll)) = x(j-(kk-ll)) + mm
       x(j-kk) = x(j-kk) - x(j)
       if (x(j-kk) < 0) x(j-kk) = x(j-kk) + mm
    end do
    if (mod(ss, 2) == 1) then
        do j = kk, 1, -1
           x(j+1) = x(j)
        end do
        x(1) = x(kk+1)
        x(ll+1) = x(ll+1) - x(kk+1)
        if (x(ll+1) < 0) x(ll+1) = x(ll+1) + mm
    end if
    if (ss /= 0) THEN
        ss = ishft(ss, -1)
    else
        t = t-1
    end if
    if (t > 0) GO TO 10

    do j = 1, ll
       ranx(j+kk-ll) = x(j)
    end do
    do j = ll+1, kk
       ranx(j-ll) = x(j)
    end do

    do j = 1, 10
       call rnarry(x, kkk)
    end do

    return
  end subroutine rnstrt
end module
