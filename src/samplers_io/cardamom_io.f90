
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
  ! See function/subroutine specific comments for exceptions and contributors
  !!!!!!!!!!!

  ! Module contains subroutines and variables needed to output parameter, 
  ! likelihood and step size information from the MHMCMC.

  implicit none

  ! declare private
  private

  ! allow access to specific functions
  public:: write_mcmc_output               &
           ,write_parameters                &
           ,write_variances                 &
           ,write_covariance_matrix         &
           ,write_covariance_info           &
           ,check_for_existing_output_files &
           ,open_output_files               &
           ,close_output_files              

  ! allow access to needed variable
  public:: restart_flag

  ! declare module level variables
  integer:: pfile_unit = 10, sfile_unit = 11, cfile_unit = 12, cifile_unit = 13, ifile_unit = 14
  ! default assumption is that this is not a restart fun
  logical:: restart_flag = .false.

  ! parameters
  integer, parameter:: real_bytes = 8  ! number of bytes in real variable, 8 bytes is to make double precision

  type io_buffer_space
    integer:: io_buffer, io_buffer_count
    double precision, allocatable, dimension(:,:) :: &
                                    variance_buffer, &
                                   mean_pars_buffer, &
                                        pars_buffer

    double precision, allocatable, dimension(:) :: &
                                   nsample_buffer, &
                               accept_rate_buffer, &
                                      prob_buffer
  end type

  save

  contains
  
  !
  !------------------------------------------------------------------
  !
  subroutine check_for_existing_output_files(npars, nOUT, nWRITE, sub_fraction &
                                            ,parname, stepname, covname, covinfoname)

    ! subroutine checks whether both the parameter and step files exist for this
    ! job. If they do we will assume that this is a restart job that we want to
    ! finish off. Important for large jobs or running on machines with may crash
    ! / have runtime limits
    implicit none
    ! declare input variables
    integer, intent(in):: npars, nOUT, nWRITE
    double precision, intent(in):: sub_fraction
    character(350), intent(in):: parname, stepname, covname, covinfoname
    ! local variables
    logical:: par_exists, step_exists, cov_exists, covinfo_exists
    double precision:: dummy
    integer:: num_lines, status

    ! Check that all files exist
    inquire(file = trim(parname),     exist = par_exists)
    inquire(file = trim(stepname),    exist = step_exists)
    inquire(file = trim(covname),     exist = cov_exists)
    inquire(file = trim(covinfoname), exist = covinfo_exists)

    ! now determine the correct response
    if (par_exists .and. step_exists .and. cov_exists .and. covinfo_exists) then

        ! All files exist therefore this might be a restart run.
        ! lets see if there is anything in the files that we might use
        ! count the number of remaining lines in the file..
        ! open the relevant output files
        call open_output_files(parname, stepname, covname, covinfoname)
        status = 0; num_lines = 0
        do
          read(pfile_unit, iostat = status) dummy
          if ( status .ne. 0 ) exit
          num_lines = num_lines+1
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
            ! The file exists but is empty/no enough so treat it as a fresh start
            restart_flag = .false.
            print*,"Output files are present, however they are too small for a restart"
        endif
        ! Either way we open the file up later on so now we need to close them
        call close_output_files

    else  ! par_exists .and. step_exists

        ! Then or of these files exists and the other does not so it is
        ! ambiguous whether or not this is a restart
        print*,"One or more of the analysis files cannot be found."
        print*,"CARDAMOM must start from scratch... "
        restart_flag = .false.

    endif  ! par_exists .and. step_exists

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
  
  
  subroutine open_output_files(parname, stepname, covname, covinfoname)

    ! Subroutine opens the needed output files and destroys any previously
    ! existing files with the same name, just in case mind!
    ! NOTE: that is unless I have not remove the 'UNKNOWN' status in which case
    ! then the files are appended to

    implicit none

    ! declare input variables
    character(350), intent(in):: parname, stepname, covname, covinfoname

    ! declare local variables
    integer:: ios, reclen
    double precision:: a = 1d0

    ! open files now
    ! most of these will require new information to be appended to the end at
    ! all times-therefore we use the unformatted stream access
    open(pfile_unit, file = trim(parname), form="UNFORMATTED",access="stream",status="UNKNOWN",iostat = ios)
    if (ios /= 0) print*,"error ",ios, " opening file",trim(parname)
    open(sfile_unit, file = trim(stepname), form="UNFORMATTED",access="stream",status="UNKNOWN",iostat = ios)
    if (ios /= 0) print*,"error ",ios, " opening file",trim(stepname)
    open(cifile_unit, file = trim(covinfoname), form="UNFORMATTED",access="stream",status="UNKNOWN",iostat = ios)
    if (ios /= 0) print*,"error ",ios, " opening file",trim(covinfoname)
    ! for the covariance matrix we have a fixed size containing two matrices, 
    ! the initial and the current output-therefore we use
    inquire(iolength = reclen) a !; print*,reclen
    open(cfile_unit, file = trim(covname), form="UNFORMATTED",access="direct",recl = reclen, iostat = ios)
    if (ios /= 0) print*,"error ",ios, " opening file",trim(covname)


  end subroutine open_output_files

  subroutine initialize_buffers(npars, n_write_events, io_space)
    integer, intent(in) :: npars, n_write_events
    type(io_buffer_space), intent(inout) :: io_space
    ! TODO oops each chain needs its own buffer obeject and file streams
    ! also initialize buffers
    ! TODO move because not fitting function name
   ! Initialise counters used to track the output of parameter sets
    io_space%io_buffer_count = 0 
    io_space%io_buffer = min(1000, max(10, n_write_events / 10))
   
    ! Allocate variables used in io buffering, 
    ! these could probably be moved to a more sensible place within cardamom_io.f90
    allocate(io_space%variance_buffer(npars, io_space%io_buffer), &
             io_space%mean_pars_buffer(npars, io_space%io_buffer), &
             io_space%pars_buffer(npars, io_space%io_buffer), &
             io_space%prob_buffer(io_space%io_buffer), &
             io_space%nsample_buffer(io_space%io_buffer), &
             io_space%accept_rate_buffer(io_space%io_buffer))
 

    return

  end subroutine 

 
  !
  !------------------------------------------------------------------
  !
  subroutine write_covariance_matrix(covariance, npars, initial_cov)

    ! subroutine writes MCMC accepted parameters and step values to binary files

    implicit none

    ! arguments
    logical, intent(in):: initial_cov
    integer, intent(in):: npars
    double precision, dimension(npars, npars), intent(in):: covariance

    ! declare local variables
    integer:: i, j, irec

    ! If we have already written the initial covariance matrix we want to keep
    ! over-writing the current matrix. We do this to avoid large files form
    ! writing out multiple covariance matrices
    if (.not.initial_cov) then
        irec = npars*npars
    else
        irec = 0
    end if

    ! write out the file. Its binary format has already been determined at the
    ! openning of the file

    do i = 1, npars
       do j = 1, npars
          irec = irec+1
          write(cfile_unit, rec = irec) covariance(i, j)
       end do
    end do

    return

  end subroutine write_covariance_matrix
  !
  !------------------------------------------------------------------
  !
  subroutine write_covariance_info(mean_pars, nsample, npars)

    ! subroutine writes MCMC accepted parameters and step values to binary files

    implicit none

    ! arguments
    integer, intent(in):: npars
    double precision, intent(in):: nsample
    double precision, dimension(npars), intent(in):: mean_pars

    ! declare local variables
    integer:: i, j

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
  subroutine write_variances(variance, npars, accept_rate)

    ! subroutine writes parameter variance for corresponding parameter values

    implicit none

    ! declare input variables
    integer, intent(in):: npars
    double precision, dimension(npars), intent(in):: variance
    double precision, intent(in):: accept_rate  ! local acceptance rate

    ! declare local variables
    integer:: n

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
  subroutine write_parameters(pars, prob, npars)

    ! subroutine writes parameter values to binary file`

    implicit none

    ! declare input variables
    integer, intent(in):: npars
    double precision, dimension(npars), intent(in):: pars
    double precision, intent(in):: prob

    ! declare local variables
    integer:: n

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
                               pars, prob, npars, dump_now, io_space)

    ! Arguments
    integer, intent(in):: npars
    double precision, dimension(npars, npars), intent(in):: covariance
    double precision, dimension(npars), intent(in):: mean_pars, &
                                                       variance, &
                                                           pars
    double precision, intent(in):: nsample, accept_rate, prob
    logical, intent(in):: dump_now
    type(io_buffer_space), intent(inout) :: io_space

    ! Local variables
    integer:: i

!    ! Debugging print statements
!    print*,"write_mcmc_output:"

    ! Increment buffer
    io_space%io_buffer_count = io_space%io_buffer_count+1

    ! Store information in buffer for later writing
    io_space%variance_buffer(1:npars, io_space%io_buffer_count) = variance
    io_space%mean_pars_buffer(1:npars, io_space%io_buffer_count) = mean_pars
    io_space%pars_buffer(1:npars, io_space%io_buffer_count) = pars
    io_space%prob_buffer(io_space%io_buffer_count) = prob
    io_space%nsample_buffer(io_space%io_buffer_count) = nsample
    io_space%accept_rate_buffer(io_space%io_buffer_count) = accept_rate

    ! Are we storing information in buffer or writing to file?
    if (io_space%io_buffer_count == io_space%io_buffer .or. dump_now) then

        ! Then we are writing out to file
        ! Only write the most current covariance matrix as this would be an overwrite anyway
        call write_covariance_matrix(covariance, npars, .false.)
        ! Everything else loop through the buffered output to write out
        do i = 1, io_space%io_buffer_count
           call write_covariance_info(io_space%mean_pars_buffer(:,i), io_space%nsample_buffer(i), npars)
           call write_variances(io_space%variance_buffer(:,i), npars, io_space%accept_rate_buffer(i))
           call write_parameters(io_space%pars_buffer(:,i), io_space%prob_buffer(i), npars)
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
