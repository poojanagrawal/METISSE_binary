subroutine zcnsts(z,zpars)
    use track_support
    use z_support
    use remnant_support

    real(dp), intent(in) :: z
    real(dp), intent(out) :: zpars(20)

    integer :: i,ierr
    logical :: debug
    character(len=strlen):: path
    
    ierr = 0
    debug = .false.
    
    if (direct_call .and. (.not. defined(z)))then
        print*,"Error: initial_Z is not defined "
        STOP
    endif
    
    initial_Z = z
    
    !reading defaults option first
    call read_defaults(ierr); if (ierr/=0) STOP

    !read inputs from evolve_metisse.in
    call read_input(ierr); if (ierr/=0) STOP
    
    !read metallicity related variables
    call get_metallcity_file_from_Z(z,ierr); if (ierr/=0) STOP
    
    
    !reading format file
    call read_format(format_file,ierr); if (ierr/=0) STOP

    !getting filenames
    call get_files_from_path(INPUT_FILES_DIR,ierr); if (ierr/=0) STOP

    if (verbose) print*,"Number of input tracks: ", num_tracks

    call read_key_eeps()
    if (debug) print*, "key eeps", key_eeps

    if (read_eep_files) then
        if (debug) print*,"reading eep files"
        do i=1,num_tracks
            call read_eep(s(i))
            if(debug) write(*,'(a50,f8.2,99i8)') trim(s(i)% filename), s(i)% initial_mass, s(i)% ncol
        end do
    else
        !read and store column names in temp_cols from the the file if header location is not provided
        if (header_location<=0) then
            if (debug) print*,"Reading column names from file"

            call process_columns(column_name_file,temp_cols,ierr)
            
            if(ierr/=0) then
                print*,"Failed while trying to read column_name_file"
                print*,"Check if header location and column_name_file are correct "
                STOP
            endif

            if (size(temp_cols) /= total_cols) then
                print*,'Number of columns in the column_name_file does not matches with the total_cols'
                print*,'Check if column_name_file and total_cols are correct'
                STOP
            endif
        end if

        do i=1,num_tracks
            call read_input_file(s(i))
            if(debug) write(*,'(a50,f8.2,99i8)') trim(s(i)% filename), s(i)% initial_mass, s(i)% ncol
        end do

    endif

    if(debug) print*, s(1)% cols% name, s(1)% tr(:,1)
    
    do i = 1,size(s)
        s(i)% has_mass_loss = check_mass_loss(s(i))
    end do

    !TODO: check for monotonicity of initial masses
    if (debug) print*,s% initial_mass

    !first calculate it SSE way as a default
    call calculate_sse_zpars(z,zpars)

    !then reset z parameters where available
    !and determine cutoff masses
    call set_zparameters(zpars)

    if (direct_call) then
    ! sets remnant schmeme from SSE_input_controls
        call set_remnant_scheme()
    else
    ! reads
        call assign_commons()
    endif

end subroutine zcnsts

