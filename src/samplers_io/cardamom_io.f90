
module cardamom_io

  !!!!!!!!!!!
  !
  ! Samplers_io functions split from cardamom_io.f90 
  ! 
  ! Authorship contributions
  !
  ! This code is based on the original C verion of the University of Edinburgh
  ! CARDAMOM framework created by A. A. Bloom (now at the Jet Propulsion Laboratory).
  ! All code translation into Fortran, integration into the University of
  ! Edinburgh CARDAMOM code and subsequent modifications by:
  ! T. L. Smallman (t.l.smallman@ed.ac.uk, University of Edinburgh)
  ! J. F. Exbrayat (University of Edinburgh)
  ! See function / subroutine specific comments for exceptions and contributors
  !!!!!!!!!!!

  ! Module contains subroutines and variables needed to output parameter,
  ! likelihood and step size information from the MHMCMC.

  implicit none

  ! declare private
  private

  ! allow access to specific functions
  public :: write_mcmc_output               &
           ,write_parameters                &
           ,write_variances                 &
           ,write_covariance_matrix         &
           ,write_covariance_info           &
           ,update_for_restart_simulation   &
           ,check_for_existing_output_files &
           ,open_output_files               &
           ,close_output_files              

  ! allow access to needed variable
  public :: restart_flag

  ! declare module level variables
  integer :: pfile_unit = 10, sfile_unit = 11, cfile_unit = 12, cifile_unit = 13, ifile_unit = 14
  ! default assumption is that this is not a restart fun
  logical :: restart_flag = .false.

  ! parameters
  integer, parameter :: real_bytes = 8 ! number of bytes in real variable, 8 bytes is to make double precision

  ! TODO taken from cardamom_structures.f90, where to park this in the future?
  type io_buffer_space
    integer :: io_buffer, io_buffer_count
    double precision, allocatable, dimension(:,:) :: &
                                    variance_buffer, &
                                   mean_pars_buffer, &
                                        pars_buffer

    double precision, allocatable, dimension(:) :: &
                                   nsample_buffer, &
                               accept_rate_buffer, &
                                      prob_buffer
  end type
  type(io_buffer_space) :: io_space

  save

  contains
  
  !
  !------------------------------------------------------------------
  !
  subroutine check_for_existing_output_files(npars,nOUT,nWRITE,sub_fraction &
                                            ,parname,stepname,covname,covinfoname)

    ! subroutine checks whether both the parameter and step files exist for this
    ! job. If they do we will assume that this is a restart job that we want to
    ! finish off. Important for large jobs or running on machines with may crash
    ! / have runtime limits
    implicit none
    ! declare input variables
    integer, intent(in) :: npars, nOUT, nWRITE
    double precision, intent(in) :: sub_fraction
    character(350), intent(in) :: parname, stepname, covname, covinfoname
    ! local variables
    logical :: par_exists, step_exists, cov_exists, covinfo_exists
    double precision :: dummy
    integer :: num_lines, status

    ! Check that all files exist
    inquire(file=trim(parname),     exist=par_exists)
    inquire(file=trim(stepname),    exist=step_exists)
    inquire(file=trim(covname),     exist=cov_exists)
    inquire(file=trim(covinfoname), exist=covinfo_exists)

    ! now determine the correct response
    if (par_exists .and. step_exists .and. cov_exists .and. covinfo_exists) then

        ! All files exist therefore this might be a restart run.
        ! lets see if there is anything in the files that we might use
        ! count the number of remaining lines in the file..
        ! open the relevant output files
        call open_output_files(parname,stepname,covname,covinfoname)
        status = 0 ; num_lines = 0
        do
          read(pfile_unit,iostat=status) dummy
          if ( status .ne. 0 ) exit
          num_lines = num_lines + 1
        enddo
        ! Re-use dummy to calculate the target file size to be considered for
        ! restart
        dummy = ((dble(nOUT)/dble(nWRITE)) * sub_fraction) * dble(npars+1)
        if (num_lines > dummy) then
            ! Then there is something in the file we we can use it
            restart_flag = .true.
            print*,"...have found parameter file = ",trim(parname)
            print*,"...have found step file = ",trim(stepname)
            print*,"...have found cov file = ",trim(covname)
            print*,"...have found cov_info file = ",trim(covinfoname)
        else
            ! The file exists but is empty / no enough so treat it as a fresh start
            restart_flag = .false.
            print*,"Output files are present, however they are too small for a restart"
        endif
        ! Either way we open the file up later on so now we need to close them
        call close_output_files

    else ! par_exists .and. step_exists

        ! Then or of these files exists and the other does not so it is
        ! ambiguous whether or not this is a restart
        print*,"One or more of the analysis files cannot be found."
        print*,"CARDAMOM must start from scratch... "
        restart_flag = .false.

    endif ! par_exists .and. step_exists

  end subroutine check_for_existing_output_files
  !
  !------------------------------------------------------------------
  !
  subroutine close_output_files

    ! where you open a file you've got to make sure that you close them too. It
    ! just tidy

    implicit none

    ! close the files we have in memory
    close(pfile_unit)
    close(sfile_unit)
    close(cfile_unit)
    close(cifile_unit)

  end subroutine close_output_files
  
  
  subroutine open_output_files(parname,stepname,covname,covinfoname)

    ! Subroutine opens the needed output files and destroys any previously
    ! existing files with the same name, just in case mind!
    ! NOTE: that is unless I have not remove the 'UNKNOWN' status in which case
    ! then the files are appended to

    implicit none

    ! declare input variables
    character(350), intent(in) :: parname, stepname, covname,covinfoname

    ! declare local variables
    integer :: ios, reclen
    double precision :: a = 1d0

    ! open files now
    ! most of these will require new information to be appended to the end at
    ! all times - therefore we use the unformatted stream access
    open(pfile_unit,file=trim(parname),form="UNFORMATTED",access="stream",status="UNKNOWN",iostat=ios)
    if (ios /= 0) print*,"error ",ios," opening file",trim(parname)
    open(sfile_unit,file=trim(stepname),form="UNFORMATTED",access="stream",status="UNKNOWN",iostat=ios)
    if (ios /= 0) print*,"error ",ios," opening file",trim(stepname)
    open(cifile_unit,file=trim(covinfoname),form="UNFORMATTED",access="stream",status="UNKNOWN",iostat=ios)
    if (ios /= 0) print*,"error ",ios," opening file",trim(covinfoname)
    ! for the covariance matrix we have a fixed size containing two matrices,
    ! the initial and the current output - therefore we use
    inquire(iolength = reclen) a !; print*,reclen
    open(cfile_unit,file=trim(covname),form="UNFORMATTED",access="direct",recl=reclen,iostat=ios)
    if (ios /= 0) print*,"error ",ios," opening file",trim(covname)

    return

  end subroutine open_output_files

  !
  !------------------------------------------------------------------
  !

  subroutine read_options(solutions_wanted,freq_print,freq_write,outfile)
    use MCMCOPT, only: MCO, MCOUT

    ! loads required options about the MHMCMC either form hardcoded sources or
    ! from variables which were read form the command line

    implicit none

    ! declare input variables
    character(350), intent(in) :: outfile
    integer, intent(in) :: solutions_wanted, freq_print, freq_write

    ! defining hardcoded MCMC options
    MCO%append = 1
    MCO%nADAPT = 1000 ! TLS: 500 -> 1000 -> 5000 -> 10000
    MCO%fADAPT = 0.5d0
    MCO%randparini = .false.
    MCO%returnpars = .false.
    MCO%fixedpars  = .true. ! TLS: changed from .false. for testing 16/12/2019

    ! command line options

    ! how many accepted parameters to accept before completion
    if (solutions_wanted > 0) then
        MCO%nOUT = solutions_wanted
    else
        MCO%nOUT = 1000
        write(*,*)"Default MCO%nOUT value used"
    endif
    ! how frequently to print information to screen
    if (freq_print >= 0) then
        MCO%nPRINT = freq_print
    else
        MCO%nPRINT = 1000
        write(*,*)"Default MCO%nPRINT value used"
    end if
    ! how frequently to write information to file
    if (freq_write >= 0) then
        MCO%nWRITE = freq_write
    else
        MCO%nWRITE = 1000
        write(*,*)"Default MCO%nWRITE value used"
    end if

    ! Assume that sub-sampling process, if completed, will use 10 % of the
    ! simulation time therefore we want to adjust the output frequency to
    ! correct for this
    MCO%nOUT = max(1,MCO%nOUT - MCOUT%nos_iterations)

    ! construct file names
    write(MCO%outfile,fmt='(A)')trim(outfile)//"PARS"
    write(MCO%stepfile,fmt='(A)')trim(outfile)//"STEP"
    write(MCO%covfile,fmt='(A)')trim(outfile)//"COV"
    write(MCO%covifile,fmt='(A)')trim(outfile)//"COVINFO"

  end subroutine read_options
  !
  !-------------------------------------------------------------------
  !
  subroutine update_for_restart_simulation()
    use MCMCOPT, only: PI, MCOUT, MCO
    ! TODO DATAin, an object of type DATAtype saving soring simulation data, is used only for n_pars.  
    ! find alternative source of PI%N_pars
    use math_functions, only: std, covariance_matrix, inverse_matrix, par2nor

    ! subroutine is responsible for loading previous parameter and step size
    ! information into the current

    implicit none

    ! local variables
    integer :: a, b, c, i, j, num_lines, status, n_pars
    double precision :: dummy
    double precision,dimension(:,:), allocatable :: tmp

    ! the parameter and step files should have already been openned so
    ! read the parameter and step files to get to the end

    ! rewind to the beginning
    rewind(pfile_unit) ; rewind(sfile_unit) ; rewind(cifile_unit)

    !
    ! As this subroutine will only be called once reading each file will occur
    ! separately to improve simplicity.
    !

    !
    ! Parameter file - stored as non-normalised values
    !

    n_pars= PI%npars

    ! count the number of lines in the file..
    status = 0 ; num_lines = 0
    do
      read(pfile_unit,iostat=status) dummy
      if ( status .ne. 0 ) exit
      num_lines = num_lines + 1
    enddo
    ! Determine the number of complete parameter vectors stored. Note that the +
    ! 1 is due to the log-likelihood score being saved as well.
    num_lines = num_lines/(n_pars+1)

    ! Allocate memory to our temperary variable and the normalised parameter
    ! vector equivalent.
    allocate(tmp(num_lines,(n_pars+1)))
    ! rewind so that we can read the contents now correctly
    rewind(pfile_unit)

    ! Read the data for real
    do i = 1, num_lines
       do j = 1, (n_pars+1)
          read(pfile_unit) tmp(i,j)
       end do ! j for parameter
    end do ! i for combinations

    ! Determine the total number of iterations processed so far
    MCOUT%nos_iterations = num_lines * MCO%nWRITE
    ! Extract the final parameter set and load into the initial parameter vector
    ! for the analysis. NOTE: This parameter set will be normalised on entry
    ! into the MHMCMC subroutine
    PI%parini(1:n_pars) = tmp(num_lines,1:n_pars)

    ! free up variable for new file
    deallocate(tmp)

    !
    ! Variance file - stores output of the current parameter variance
    !

    ! count the number of remaining lines in the file..
    status = 0 ; num_lines = 0
    do
      read(sfile_unit,iostat=status) dummy
       if ( status .ne. 0. ) exit
       num_lines = num_lines + 1
    enddo

    ! Determine the number of actual stepsize vectors present.
    ! The + 1 is due to the local acceptance rate being provided too.
    num_lines = num_lines/(n_pars+1)
    ! allocate memory
    allocate(tmp(num_lines,(n_pars+1)))
    ! rewind, for actual reading
    rewind(sfile_unit)

    ! now read the data for real
    do i = 1, num_lines
       do j = 1, (n_pars+1)
          read(sfile_unit) tmp(i,j)
       end do ! j for parameter
    end do ! i for combinations
    ! Save the current acceptance_rate
    MCOUT%acceptance_rate = MCOUT%nos_iterations * MCO%nWRITE

    ! tidy up for the next file
    deallocate(tmp)

    !
    ! Covariance matrix file
    !

    ! The covariance matrix may contain either 1 or 2 complete matrices. We will
    ! just want to the latest one.

    ! count the number of remaining lines in the file..
    status = 0 ; num_lines = 1
    do
       read(cfile_unit,iostat=status,rec = num_lines) dummy
       if ( status .ne. 0. ) exit
       num_lines = num_lines + 1
    enddo

    ! Determine whether there is 1 or more matrice here
    if ((num_lines/n_pars)/n_pars == 1) then
        ! the size of the file is consistent with a single matrix having been
        ! saved
        a = 1
    else if ((num_lines/n_pars)/n_pars == 2) then
        !
        a = 2
    else
        ! something has gone wrong - best stop
        print*,"Error reading COV file"
        print*,"DATAin%nopars = ",n_pars,"COV length = ",num_lines *n_pars
        stop
    endif

    ! now read the data for real
    c = 1
    do b = 1, a
       do i = 1, n_pars
          do j = 1, n_pars
             read(cfile_unit, rec = c) PI%covariance(i,j)
             c = c + 1
          end do ! j for parameter
       end do ! i for combinations
    end do

    ! extract current variance information
    do i = 1, PI%npars
       PI%parvar(i) = PI%covariance(i,i)
    end do
    ! estimate status of the inverse covariance matrix
    call inverse_matrix( PI%npars, PI%covariance, PI%iC )

    !
    ! Covariance information file
    !

    ! The number of parameters on which the covariance matrix is based must be
    ! known to allow for correct updating. Similarly the mean normalised
    ! parameter values are also needed

    ! count the number of remaining lines in the file..
    status = 0 ; num_lines = 0
    do
      read(cifile_unit,iostat=status) dummy
       if ( status .ne. 0. ) exit
       num_lines = num_lines + 1
    enddo

    ! how many parameter vectors have been output. Note the + 1 is accounting
    ! for the number of samples underlying the mean
    num_lines = num_lines / (PI%npars + 1)
    ! allocate memory
    allocate(tmp(num_lines,(n_pars+1)))
    ! rewind, for actual reading
    rewind(cifile_unit)

    ! now read the data for real
    do i = 1, num_lines
       do j = 1, (n_pars+1)
          read(cifile_unit) tmp(i,j)
       end do ! j for parameter
    end do ! i for combinations
    ! Store the most recent step size, which corresponds with the saved
    ! parmeters (above) and covariance matrix (below)
    PI%mean_par = tmp(num_lines,1:n_pars)
    PI%Nparvar = tmp(num_lines,n_pars+1)

    return

  end subroutine update_for_restart_simulation
  !
  !------------------------------------------------------------------
  !
  subroutine write_covariance_matrix(covariance,npars,initial_cov)

    ! subroutine writes MCMC accepted parameters and step values to binary files

    implicit none

    ! arguments
    logical, intent(in) :: initial_cov
    integer, intent(in) :: npars
    double precision, dimension(npars,npars), intent(in) :: covariance

    ! declare local variables
    integer :: i,j,irec

    ! If we have already written the initial covariance matrix we want to keep
    ! over-writing the current matrix. We do this to avoid large files form
    ! writing out multiple covariance matrices
    if (.not.initial_cov) then
        irec = npars * npars
    else
        irec = 0
    end if

    ! write out the file. Its binary format has already been determined at the
    ! openning of the file

    do i = 1, npars
       do j = 1, npars
          irec = irec + 1
          write(cfile_unit, rec = irec) covariance(i,j)
       end do
    end do

    return

  end subroutine write_covariance_matrix
  !
  !------------------------------------------------------------------
  !
  subroutine write_covariance_info(mean_pars,nsample,npars)

    ! subroutine writes MCMC accepted parameters and step values to binary files

    implicit none

    ! arguments
    integer, intent(in) :: npars
    double precision, intent(in) :: nsample
    double precision, dimension(npars), intent(in) :: mean_pars

    ! declare local variables
    integer :: i,j

    ! write out the file. Its binary format has already been determined at the
    ! openning of the file

    do i = 1, npars
       write(cifile_unit) mean_pars(i)
    end do

    write(cifile_unit) nsample

    return

  end subroutine write_covariance_info
  !
  !------------------------------------------------------------------
  !
  subroutine write_variances(variance,npars,accept_rate)

    ! subroutine writes parameter variance for corresponding parameter values

    implicit none

    ! declare input variables
    integer, intent(in) :: npars
    double precision, dimension(npars), intent(in) :: variance
    double precision, intent(in) :: accept_rate ! local acceptance rate

    ! declare local variables
    integer :: n

    ! write out the file. Its binary format has already been determined at the
    ! openning of the file

    do n = 1, npars
       write(sfile_unit) variance(n)
    end do

    ! we will need to know the current acceptance rate for restarts
    write(sfile_unit) accept_rate

    return

  end subroutine write_variances
  !
  !------------------------------------------------------------------
  !
  subroutine write_parameters(pars,prob,npars)

    ! subroutine writes parameter values to binary file`

    implicit none

    ! declare input variables
    integer, intent(in) :: npars
    double precision, dimension(npars), intent(in) :: pars
    double precision, intent(in) :: prob

    ! declare local variables
    integer :: n

    ! write out the file. Its binary format has already been determined at the
    ! openning of the file

    do n = 1, npars
       write(pfile_unit) pars(n)
    end do

    ! now add the probability
    write(pfile_unit) prob

    ! close will occur at the end of the MCMC

    ! return back
    return

  end subroutine write_parameters
  !
  !--------------------------------------------------------------------
  !
  subroutine write_mcmc_output(variance, accept_rate, &
                               covariance, mean_pars, nsample, &
                               pars, prob, npars, dump_now)

    ! Arguments
    integer, intent(in) :: npars
    double precision, dimension(npars,npars), intent(in) :: covariance
    double precision, dimension(npars), intent(in) :: mean_pars, &
                                                       variance, &
                                                           pars
    double precision, intent(in) :: nsample, accept_rate, prob
    logical, intent(in) :: dump_now

    ! Local variables
    integer :: i

!    ! Debugging print statements
!    print*,"write_mcmc_output:"

    ! Increment buffer
    io_space%io_buffer_count = io_space%io_buffer_count + 1

    ! Store information in buffer for later writing
    io_space%variance_buffer(1:npars,io_space%io_buffer_count) = variance
    io_space%mean_pars_buffer(1:npars,io_space%io_buffer_count) = mean_pars
    io_space%pars_buffer(1:npars,io_space%io_buffer_count) = pars
    io_space%prob_buffer(io_space%io_buffer_count) = prob
    io_space%nsample_buffer(io_space%io_buffer_count) = nsample
    io_space%accept_rate_buffer(io_space%io_buffer_count) = accept_rate

    ! Are we storing information in buffer or writing to file?
    if (io_space%io_buffer_count == io_space%io_buffer .or. dump_now) then

        ! Then we are writing out to file
        ! Only write the most current covariance matrix as this would be an overwrite anyway
        call write_covariance_matrix(covariance,npars,.false.)
        ! Everything else loop through the buffered output to write out
        do i = 1, io_space%io_buffer_count
           call write_covariance_info(io_space%mean_pars_buffer(:,i),io_space%nsample_buffer(i),npars)
           call write_variances(io_space%variance_buffer(:,i),npars,io_space%accept_rate_buffer(i))
           call write_parameters(io_space%pars_buffer(:,i),io_space%prob_buffer(i),npars)
        end do

        ! Reset buffer increment
        io_space%io_buffer_count = 0

    endif

!    ! Debugging print statements
!    print*,"write_mcmc_output:done"

  end subroutine write_mcmc_output
  !
  !--------------------------------------------------------------------
  !
end module cardamom_io
