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
    call read_input(ierr);if (ierr/=0) STOP
    
    !read metallicity related variables
    call read_metallicity_file(metallicity_file,ierr);if (ierr/=0) STOP
    
     !get folder name if read_files_from_Z
    if (read_files_from_Z) then
        if (Z_folder_list == '') then
            print*,"Error: Z_folder_list not defined for read_files_from_Z"
            STOP
        else
            call get_folder_from_Z(INPUT_FILES_DIR,initial_Z,path)
        endif
    else
        path = INPUT_FILES_DIR
    endif
    
    !reading format file
    call read_format(format_file,ierr)

    !getting filenames
    call get_files_from_path(path)

    if (verbose) print*,"Number of input tracks: ", num_tracks

    !determine key columns 
    if (key_columns_file /= '') call process_columns(key_columns_file,key_cols,ierr)
    if (debug) print*,"Using key columns: ", key_cols % name

    !get column numbers for core related quantities
!    call get_core_columns()

    if (.not. read_eep_files .and. header_location<=0)then
            if (debug) print*,"Reading column names from file"

            call process_columns(column_name_file,temp_cols,ierr)
            if(ierr/=0) then
                print*, ""
                print*,"Check if header location and column_name_file are correct "
                STOP
            endif

            if (size(temp_cols) /= total_cols) then
                print*,'Erorr reading number of columns'
                print*,'Check if column_name_file and total_cols are correct'
                STOP
            endif
    end if

    call read_key_eeps()
    if (debug) print*, "key eeps", key_eeps

        if (read_eep_files) then
            if (debug) print*,"reading eep files"
            do i=1,num_tracks
                call read_eep(s(i))
                if(debug) write(*,'(a50,f8.2,99i8)') trim(s(i)% filename), s(i)% initial_mass, s(i)% ncol
            end do
        else
            do i=1,num_tracks
                call read_input_file(s(i))
                if(debug) write(*,'(a50,f8.2,99i8)') trim(s(i)% filename), s(i)% initial_mass, s(i)% ncol
            end do

        endif

        if (size(key_cols)<1) then
            allocate(key_cols(s(1)% ncol))
            key_cols% name = s(1)% cols% name
        end if

    if(debug) print*, s(1)% cols% name, s(1)% tr(:,1)

    !locate columns of mass, age etc.
    call locate_column_numbers(s,key_cols)
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
        call set_remnant_scheme()
    else
        call assign_commons()
    endif

end subroutine zcnsts

