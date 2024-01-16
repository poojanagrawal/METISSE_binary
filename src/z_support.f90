module z_support
    use track_support
    implicit none

    character(LEN=strlen) :: INPUT_FILES_DIR
    logical :: read_eep_files, read_all_columns

    integer :: max_metallicity_files = 50
    character(LEN=strlen) :: metallicity_file,format_file, extra_columns_file
    character(LEN=strlen), allocatable :: metallicity_file_list(:)
    character(LEN=col_width) :: extra_columns(100)

    real(dp) :: Z_files, Z_accuracy_limit

    !format_specifications
    character(LEN=5):: file_extension
    integer :: header_location, eep_location
    character(LEN=strlen) :: column_name_file
    type(column), allocatable :: key_cols(:), temp_cols(:)
    integer :: total_cols

    character:: extra_char

    real(dp) :: mass, max_age, min_mass, max_mass

    character(LEN=strlen), allocatable :: track_list(:)

    !used to set star_type_from_history
    ! central limits for high- / intermediate-mass stars, set these from input eep_controls nml
    real(dp) :: center_gamma_limit = 1d2
    real(dp) :: center_carbon_limit = 1d-4
    real(dp) :: log_center_T_limit = 9d0
    real(dp) :: high_mass_limit = 1d1 !Msun
    real(dp) :: he_core_mass_limit = 2.2

!    real(dp) ::  Mup_core,Mec_core
    integer, allocatable :: m_cutoff(:)

    type critical_mass
        integer :: loc
        real(dp) :: mass
    end type critical_mass

    type(critical_mass) :: Mcrit(9)

    namelist /SSE_input_controls/ initial_Z, max_age,read_mass_from_file,&
                        input_mass_file, number_of_tracks, max_mass, min_mass, &
                        WD_mass_scheme,use_initial_final_mass_relation, allow_electron_capture, &
                        BHNS_mass_scheme, max_NS_mass,pts_1, pts_2, pts_3, write_track_to_file

    namelist /METISSE_input_controls/ metallicity_file_list, Z_accuracy_limit, &
                        mass_accuracy_limit, construct_wd_track, verbose, &
                        write_eep_file,write_error_to_file
            
    namelist /metallicity_controls/ INPUT_FILES_DIR, Z_files,format_file, extra_columns_file, &
                        read_all_columns, extra_columns,Mhook, Mhef, Mfgb, Mup, Mec, Mextra, Z_H, Z_He
                        
    namelist /format_controls/ file_extension, read_eep_files,total_cols,&
                        extra_char, header_location, column_name_file, &
                        PreMS_EEP, ZAMS_EEP, IAMS_EEP, TAMS_EEP, BGB_EEP, cHeIgnition_EEP, &
                        cHeBurn_EEP, TA_cHeB_EEP, TPAGB_EEP, cCBurn_EEP, post_AGB_EEP, WD_EEP, &
                        Initial_EEP, Final_EEP, Extra_EEP1 ,Extra_EEP2, Extra_EEP3, &
                        fix_track, low_mass_final_eep, high_mass_final_eep, lookup_index, &
                        age_colname, mass_colname, log_L_colname ,log_T_colname, &
                        log_R_colname, he_core_mass, c_core_mass, &
                        log_Tc, c12_mass_frac, o16_mass_frac,he4_mass_frac, &
                        Lum_colname, Teff_colname, Radius_colname,&
                        log_mdot_colname, mdot_colname, he_core_radius, co_core_radius,&
                        mass_conv_envelope, radius_conv_envelope, moment_of_inertia

    contains

    subroutine read_defaults(ierr)
        integer, intent(out) :: ierr

        ierr = 0
        if (front_end <0) then
            print*, 'Error: front_end is not initialized for METISSE'
            ierr = 1
            return
        endif
        
        allocate(metallicity_file_list(max_metallicity_files))

        include 'defaults/evolve_metisse_defaults.inc'
        
        ! Note that unlike other variables, metallicity_file_list = ''
        ! in the namelist only sets the value of the first element and not the whole array
        ! that's why the default for metallicity is specified here
        ! Todo: this might be removable now
        metallicity_file_list = ''
        extra_columns = ''
        
        !initialize metallicity related variables
        include 'defaults/metallicity_defaults.inc'

        !initialize file format specs
        include 'defaults/format_defaults.inc'
    
    end subroutine read_defaults
    
    
    subroutine read_metisse_input(ierr)
        integer :: io
        integer, intent(out) :: ierr
        character(len=strlen) :: infile


        ierr = 0
        io = alloc_iounit(ierr)
        !reading user input
        
        infile = trim(METISSE_DIR)// '/evolve_metisse.in'
        open(io,FILE=infile,action="read",iostat=ierr)
            if (ierr /= 0) then
               print*, 'Error: Failed to open evolve_metisse.in'
               call free_iounit(io)
               return
            end if
            if (front_end == main) then
                read(unit = io, nml = SSE_input_controls)
                if (.not. defined(initial_Z ))then
                    print*,"Error: initial_Z is not defined"
                    ierr = 1
                    return
                endif
            endif
            read(unit = io, nml = METISSE_input_controls)
        close(io)
        call free_iounit(io)
        
    end subroutine read_metisse_input
    
    
    subroutine get_metisse_input(path_to_tracks)

    character(LEN=strlen) :: path_to_tracks
    character(LEN=strlen), allocatable :: temp_list(:)

    integer :: ierr, i

        ierr = 0
        ! use inputs from COSMIC
        call get_files_from_path(path_to_tracks,'_metallicity.in',temp_list,ierr)
        
        if (.not. allocated(temp_list)) then
            print*, 'Could not find metallicity file(s) in ',trim(path_to_tracks)
            ierr = 1
            return
        else
            if (allocated(metallicity_file_list)) deallocate(metallicity_file_list)
            allocate(metallicity_file_list(size(temp_list)))
            do i = 1, size(temp_list)
                metallicity_file_list(i) = temp_list(i)
            end do
        endif
        
    end subroutine get_metisse_input
    
    
    subroutine get_metallcity_file_from_Z(Z_req,ierr)
        real(dp), intent(in) :: Z_req
        integer, intent(out) :: ierr

        integer :: i,c
        logical:: found_z,debug
        
        debug = .false.

        if (verbose) write(*,'(a,f7.3)') 'Input Z is', Z_req
        ierr = 0
!        metallicity_file = ''
        c = 0
        found_z = .false.
        do i = 1, size(metallicity_file_list)
            if (len_trim(metallicity_file_list(i))>0) c=c+1
        end do
        
        if (c<1) then
            print*, "Error: metallicity_file_list is not defined"
            ierr = 1
            return
        endif
        
        do i = 1, size(metallicity_file_list)
            if (len_trim(metallicity_file_list(i))>0) then
                if (debug) print*, 'Reading : ', trim(metallicity_file_list(i))
                call read_metallicity_file(metallicity_file_list(i),ierr)
                if (ierr/=0) cycle
                if (.not. defined (Z_files)) then
                    print*, 'Warning: Z_files not defined in "'//trim(metallicity_file_list(i))//'"'
                else
                    if (debug) print*, 'Z_files is', Z_files
                    if (relative_diff(Z_files,Z_req) < Z_accuracy_limit) then
                        if (debug) print*, 'Z_files matches with input Z'
                        found_z = .true.
                        exit
                    endif
                endif

            endif
        end do
        if ((found_z .eqv. .false.) .and. (ierr==0)) then
            print*, 'Error: metallicity value =', Z_req, 'not found amongst given Z_files'
            print*, 'Check metallicity_file_list and value of Z_files for each file'
            print*, 'If needed, Z_accuracy_limit can be relaxed (set to a greater value).'
            ierr = 1
            return
        endif
     
    end subroutine get_metallcity_file_from_Z
    
    subroutine read_metallicity_file(filename,ierr)
        character(LEN=strlen), intent(in) :: filename

        integer :: io
        integer, intent(out) :: ierr

        ierr = 0
        
        ! reset the defaults (even if already set)
        include 'defaults/metallicity_defaults.inc'
        
        io = alloc_iounit(ierr)
        open(io, file=filename, action='read', iostat=ierr)
            if (ierr /= 0) then
               print*, 'Error: failed to open metallicity_file: "'//trim(filename)//'"'
               call free_iounit(io)
               return
            end if
            read(io, nml = metallicity_controls)
        close(io)
        call free_iounit(io)
        
    end subroutine read_metallicity_file
    
    subroutine read_format(filename,ierr)
        character(LEN=strlen), intent(in) :: filename
        integer :: io
        integer, intent(out) :: ierr
        
        ierr = 0
        io = alloc_iounit(ierr)
        !read file format specs
        open(unit=io,file=trim(filename),action='read',iostat=ierr)
            if(ierr/=0)then
                print*,'Erorr: failed to open format file: "'//trim(filename)//'"'
                print*,'check if format file is correct'
                return
            endif
            read(unit = io, nml = format_controls)
        close(io)
        call free_iounit(io)

    end subroutine read_format

    
    subroutine get_files_from_path(path,extension,file_list,ierr)
    character(LEN=strlen), intent(in) :: path
    character(LEN=*), intent(in) :: extension
    character(LEN=strlen), allocatable :: file_list(:)

    integer, intent (out) ::  ierr

    character(LEN=strlen) :: str,find_cmd
    integer :: n,i, io
    
        ierr = 0
        
        find_cmd = 'find '//trim(path)//'/*'//trim(extension)//' -maxdepth 1 > .file_name.txt'
        call system(find_cmd,ierr)
        
        if (ierr/=0) return

        io = alloc_iounit(ierr)
        open(io,FILE='.file_name.txt',action="read")

        !count the number of tracks
        n = 0
        do while(.true.)
            read(io,*,iostat=ierr)
            if(ierr/=0) exit
            n = n+1
        end do

        allocate(file_list(n))
        rewind(io)
        ierr = 0
        do i = 1,n
            read(io,'(a)',iostat=ierr)str
            if (ierr/=0) exit
            file_list(i) = trim(str)
        end do
        
        close(io)
        call free_iounit(io)
        
    end subroutine get_files_from_path

    subroutine read_eep(x)      !from iso/make_track.f90
    type(eep_track), intent(inout) :: x
    real(dp), allocatable :: temp_tr(:,:)
    integer :: ierr, io, i,j

    logical :: read_phase
    character(LEN=8) :: phase_info
    character(LEN=strlen) :: eepfile
    character(LEN=10) :: type_label
    
    logical :: debug

    ierr = 0
    debug = .false.
    read_phase = .false.

    eepfile = trim(x% filename)

    io = alloc_iounit(ierr)
    open(io,file=trim(eepfile),status='old',action='read',iostat=ierr)

    !check if the file was opened successfully; if not, then fail
    if(ierr/=0) then
       x% ignore=.true.
       write(*,*) 'PROBLEM OPENING EEP FILE: ', trim(eepfile)
       close(io)
       call free_iounit(io)
       return
    endif

    read(io,'(25x,a8)') !x% version_string
    read(io,'(25x,i8)') !x% MESA_revision_number
    read(io,*) !comment line
    read(io,*) !comment line
    read(io,'(2x,f6.4,1p1e13.5,0p3f9.2)') x% initial_Y, x% initial_Z, x% Fe_div_H, x% alpha_div_Fe, x% v_div_vcrit
    read(io,*) !comment line
    read(io,*) !comment line
    read(io,'(2x,1p1e16.10,3i8,a8,2x,a10)') x% initial_mass, x% ntrack, x% neep, total_cols, phase_info, type_label

    if (debug) print*,'reading',eepfile,x% neep
    call set_star_type_from_label(type_label,x)

    if(index(phase_info,'YES')/=0) then
       x% has_phase = .true.
       allocate(x% phase(x% ntrack))
       total_cols = total_cols - 1
    else
       x% has_phase = .false.
       total_cols = total_cols
    endif

    allocate(temp_tr(total_cols, x% ntrack), temp_cols(total_cols))
    allocate(x% eep(x% neep))

    read(io,'(8x,299i8)') x% eep
    read(io,*) ! comment line
    read(io,*) ! column numbers

    !to exclude pms -read from eep(2)
    read(io,'(1x,299a32)') temp_cols% name
    do j = x% eep(1), x% ntrack
        read(io,'(1x,299(1pes32.16e3))') temp_tr(:,j)
    enddo

    close(io)
    call free_iounit(io)
    
    !determine column of mass, age etc.
    if (i_mass<1) call locate_column_numbers(temp_cols, total_cols)
    
    if (read_all_columns) then
        x% ncol = total_cols
        allocate(x% tr(x% ncol, x% ntrack),x% cols(x% ncol))
        x% cols% name = temp_cols% name
        x% tr = temp_tr
    else
        !determine key columns
        if (.not. allocated(key_cols)) call get_key_columns(temp_cols, total_cols)
        
        x% ncol = size(key_cols)
        if (debug) print*,'i', total_cols,x% ncol,x% initial_mass
        allocate(x% tr(x% ncol, x% ntrack), x% cols(x% ncol))
        do i = 1, x% ncol
            if (debug) print*, 'key column ',i,':',key_cols(i)% name,key_cols(i)% loc
            x% cols(i)% name = key_cols(i)% name
            x% tr(i,:) = temp_tr(key_cols(i)% loc,:)
        end do
    endif
    deallocate(temp_tr, temp_cols)
  end subroutine read_eep

    !adapted from read_history_file of iso_eep_support
    subroutine read_input_file(x)
        type(eep_track), intent(inout) :: x
        character(LEN=8192) :: line
        integer :: i, io, j,ierr
        real(dp), allocatable :: temp_tr(:,:)
        logical :: debug

        ierr = 0
        debug = .false.

        if (debug) print*,"in read_input_file",x% filename
        
        io = alloc_iounit(ierr)
        open(unit=io,file=trim(x% filename),status='old',action='read')
        !read lines of header as comments

        if (header_location >0)then
            do i = 1,header_location-1
                read(io,*) !header
            end do
            allocate(temp_cols(total_cols))
            !get column names
            read(io,'(a)') line
            do i =1, total_cols
                j = scan(line," ")
                temp_cols(i)% name = line(1:j)
                if (trim(temp_cols(i)% name)==extra_char) then
                    line = adjustl(line(j:))
                    j = scan(line," ")
                    temp_cols(i)% name = line(1:j)
                endif
                line = adjustl(line(j:))
            end do
        endif
        !figure out how many data lines
        j=0
        do while(.true.)
            read(io,*,iostat=ierr)
            if(ierr/=0) exit
            j=j+1
        enddo

        x% ntrack = j

        rewind(io)

        if (header_location >0)then
        !ignore file header, already read it once
            do i=1,header_location
               read(io,*) !header
            enddo
        endif

        allocate(temp_tr(total_cols, x% ntrack))

        do j=1, x% ntrack
            read(io,'(a)') line
            call split(line, temp_tr(:,j), total_cols)
        enddo

        close(io)
        call free_iounit(io)

        !determine column of mass, age etc.
        if (i_mass<1) call locate_column_numbers(temp_cols, total_cols)

        
        if (read_all_columns) then
            x% ncol = total_cols
            allocate(x% tr(x% ncol, x% ntrack),x% cols(x% ncol))
            x% cols% name = temp_cols% name
            x% tr = temp_tr
        else
            !determine key columns
            if (.not. allocated(key_cols)) call get_key_columns(temp_cols, total_cols)
            x% ncol = size(key_cols)
            if (debug) print*,'i', total_cols,x% ncol
            allocate(x% tr(x% ncol, x% ntrack), x% cols(x% ncol))
            do i = 1, x% ncol
                if (debug) print*, 'key column ',i,':',key_cols(i)% name,key_cols(i)% loc
                x% cols(i)% name = key_cols(i)% name
                x% tr(i,:) = temp_tr(key_cols(i)% loc,:)
            end do
        endif

        deallocate(temp_tr)
        if(header_location >0) deallocate(temp_cols)

        x% neep = count(key_eeps .le. x% ntrack,1)
        allocate(x% eep(x% neep))
        x% eep = pack(key_eeps,mask = key_eeps .le. x% ntrack)
        !print*,x% eep
        x% initial_mass = x% tr(i_mass,1)
        x% initial_Z = initial_Z

        if (debug) print*,x% initial_mass, x% initial_Z, x% ncol
    end subroutine read_input_file

    !from C.Flynn's driver routine

    subroutine split(line,values,ncol)
    character(LEN=*) :: line
    real(dp) :: values(:)
    integer:: i,ncol, iblankpos
    line = adjustl(line)
        do i =1, ncol
            !print*,i,trim(line)
            iblankpos = scan(line," ")
            if (trim(line)/= '') read(line(1:iblankpos),*) values(i)
            !print*, values(i)
            line = adjustl(line(iblankpos:))
        end do
    end subroutine split

!locating essential columns here
    subroutine locate_column_numbers(cols,ncol)
        type(column), intent(in) :: cols(:)
        integer, intent(in) :: ncol
        logical :: essential
        essential = .true.
        
        ! i_age is the extra age column for recording age values of new tracks
        ! it is used for interpolating in surface quantities after any explicit mass gain/loss
        ! It is the same as i_age2 if the input tracks already include mass loss due to winds/no mass loss
        i_age = ncol+1
        
        i_age2 = locate_column(cols, age_colname, essential)
        i_mass = locate_column(cols, mass_colname, essential)
        
        if (log_L_colname /= '') then
            !find the log luminosity column
            i_logL = locate_column(cols, log_L_colname, essential)

        else
            !find the luminosity column and convert it into log
            i_lum = locate_column(cols, Lum_colname, essential)
            call make_logcolumn(s, i_logL)
        endif

        if (log_R_colname/= '') then
            i_logR = locate_column(cols, log_R_colname, essential)
        else
            i_logR = locate_column(cols, Radius_colname, essential)
            call make_logcolumn(s, i_logR)
        endif

        i_he_core = locate_column(cols, he_core_mass, essential)
        i_co_core = locate_column(cols, c_core_mass, essential)
        
        essential  = .false.
        
        !optional columns

        !TODO: - make log_T optional, Teff will get calculated in the code
        if (log_T_colname/= '') then
            i_logTe = locate_column(cols, log_T_colname, essential)
        else
            i_logTe = locate_column(cols, Teff_colname, essential)
            call make_logcolumn(s, i_logTe)
        endif

        i_RHe_core = -1
        i_RCO_core = -1

        if (he_core_radius/= '') i_RHe_core = locate_column(cols, he_core_radius)
        if (co_core_radius/= '') i_RCO_core = locate_column(cols, co_core_radius)
            

        i_mcenv = -1
        if (mass_conv_envelope/= '') i_mcenv = locate_column(cols, mass_conv_envelope)

        i_Rcenv = -1
        if (radius_conv_envelope/= '') i_rcenv = locate_column(cols, radius_conv_envelope)

!        print*, 'mcenv, rcenv columns',i_mcenv, i_Rcenv

        i_MoI = -1
        if (moment_of_inertia/= '') i_MoI = locate_column(cols, moment_of_inertia)
            
        i_he4 = locate_column(cols, he4_mass_frac)
        i_c12 = locate_column(cols, c12_mass_frac)
        i_o16 = locate_column(cols, o16_mass_frac)
        i_Tc = locate_column(cols, log_Tc)


!        if (log_mdot_colname/= '') then
!            i_mdot = locate_column(cols, log_mdot_colname)
!            call make_pow10column(s,i_mdot,"Mdot")
!        elseif (mdot_colname/= '') then
!             i_mdot =  locate_column(cols, mdot_colname)          !star_mdot
!        !            call make_logcolumn(s, i_mdot)
!        else
!             i_mdot = -1
!        endif

!        if (i_RHe_core > 0) number_of_core_columns = number_of_core_columns+1
!        if (i_RCO_core > 0) number_of_core_columns = number_of_core_columns+1

        number_of_core_columns = 5
        allocate(core_cols(number_of_core_columns))
        core_cols = -1
        core_cols(1) = i_logL
        core_cols(2) = i_he_core
        core_cols(3) = i_co_core

        if (i_RHe_core > 0) core_cols(4) = i_RHe_core
        if (i_RCO_core > 0) core_cols(5) = i_RCO_core   
    end subroutine locate_column_numbers


    integer function locate_column(cols,colname,essential)
        character(LEN=col_width), intent(in) :: colname
        type(column), intent(in) :: cols(:)

        logical, intent(in),optional :: essential
        logical :: essential1

        integer :: i

        !unless explicitly specified
        !assume that the column is not essential
        essential1 = .false.
        if (present(essential)) essential1 = essential

        !now find the column
        locate_column = -1
        if (trim(colname)=='') return
        do i=1,size(cols)
           if(adjustl(adjustr(cols(i)% name))==trim(colname)) then
              locate_column = i
              return
           endif
        enddo

        !check whether the column has been successfully located
        if(locate_column<0) then
            write(0,*) 'Could not find column: ', trim(colname)
            !STOP the code if cannot locate one of the essential columns
            if(essential1) STOP
        endif
        
    end function locate_column
      
    subroutine make_logcolumn(s, itemp)
        type(eep_track) :: s(:)
        integer :: itemp,k
        do k = 1, size(s)
            s(k)% tr(itemp,:) = log10(s(k)% tr(itemp,:))
            s(k)% cols(itemp)% name = "log("//trim(s(k)% cols(itemp)% name)//")"
        end do
    end subroutine make_logcolumn
    
    subroutine make_pow10column(s, itemp,newname)
        type(eep_track) :: s(:)
        integer :: itemp,k
        character(LEN=col_width), intent(in), optional :: newname
        do k = 1, size(s)
            s(k)% tr(itemp,:) = 10.d0**(s(k)% tr(itemp,:))
            if (present(newname)) s(k)% cols(itemp)% name = trim(newname)
        end do
    end subroutine make_pow10column

    subroutine get_key_columns(cols,ncol)
        type(column), intent(in) :: cols(:)
        integer, intent(in) :: ncol
        type(column) :: temp(ncol)
        type(column), allocatable :: temp_extra_columns(:)

        integer :: i,j,n,c,ierr
     
!        if (debug) print*, 'assigning key columns'

        ! Essential columns get reassigned here to match to reduced array format
        temp% loc = -1
        temp% name  =  ''
        
        temp(1)% loc = i_age2
        temp(1)% name = age_colname
        i_age2 = 1
        
        temp(2)% loc = i_mass
        temp(2)% name = mass_colname
        i_mass = 2
        
        temp(3)% loc = i_logL
        temp(3)% name = log_L_colname
        i_logL = 3
        
        temp(4)% loc = i_logR
        temp(4)% name = log_R_colname
        i_logR = 4
        
        temp(5)% loc = i_he_core
        temp(5)% name = he_core_mass
        i_he_core = 5
        
        temp(6)% loc = i_co_core
        temp(6)% name = c_core_mass
        i_co_core = 6
        
        
        n=7
        if (i_logTe >0) then
            temp(n)% loc = i_logTe
            temp(n)% name = log_T_colname
            i_logTe = n
            n=n+1
        endif
        
        if (i_RHe_core >0) then
            temp(n)% loc = i_RHe_core
            temp(n)% name = he_core_radius
            i_RHe_core = n
            n=n+1
        endif
        
        if (i_RCO_core >0) then
            temp(n)% loc = i_RCO_core
            temp(n)% name = co_core_radius
            i_RCO_core = n
            n=n+1
        endif
        
        if (i_mcenv>0) then
            temp(n)% loc = i_mcenv
            temp(n)% name = mass_conv_envelope
            i_mcenv = n
            n=n+1
        endif
        
        if (i_Rcenv>0) then
            temp(n)% loc = i_Rcenv
            temp(n)% name = radius_conv_envelope
            i_Rcenv = n
            n=n+1
        endif
        
        if (i_MoI>0) then
            temp(n)% loc = i_MoI
            temp(n)% name = moment_of_inertia
            i_MoI = n
            n=n+1
        endif
        
        if (i_Tc >0) then
            temp(n)% loc = i_Tc
            temp(n)% name = log_Tc
            i_Tc = n
            n=n+1
        endif
        
        if (i_he4 >0) then
            temp(n)% loc = i_he4
            temp(n)% name = he4_mass_frac
            i_he4 = n
            n=n+1
        endif
        
        if (i_c12 >0) then
            temp(n)% loc = i_c12
            temp(n)% name = c12_mass_frac
            i_c12= n
            n=n+1
        endif
        
        if (i_o16>0) then
            temp(n)% loc = i_o16
            temp(n)% name = o16_mass_frac
            i_o16 = n
            n=n+1
        endif
            
        c = 0
        do i = 1, size(extra_columns)
            if (len_trim(extra_columns(i))>0) then
                j = locate_column(cols,extra_columns(i))
                if (j>0 )then
                    temp(n)% loc = j
                    temp(n)% name = extra_columns(i)
                    n = n+1
                endif
                c=c+1
            endif
        end do
        
        ierr = 0

        if (c<1 .and. extra_columns_file /= '') call process_columns(extra_columns_file,temp_extra_columns,ierr)
        
        if(ierr/=0) print*, "Failed while trying to read extra_columns_file"

        if (allocated(temp_extra_columns)) then
            do i = 1, size(temp_extra_columns)
                if (len_trim(temp_extra_columns(i)% name)<1) exit
                temp(n)% loc = locate_column(cols,temp_extra_columns(i)% name)
                temp(n)% name = temp_extra_columns(i)% name
                n = n+1
            end do
            deallocate(temp_extra_columns)
        endif
    
        allocate(key_cols(n-1))
        key_cols% name = temp(1:n-1)% name
        key_cols% loc = temp(1:n-1)% loc

        i_age = n
    end subroutine get_key_columns


    !reading column names from file - from iso_eep_support.f90
    subroutine process_columns(filename,cols,ierr)
        character(LEN=strlen), intent(in) :: filename
        integer, intent(out) :: ierr
        integer :: i, ncols(2), nchar, column_length, pass
        character(LEN=strlen) :: line, column_name
        logical :: is_int,debug
        type(column), allocatable, intent(out) :: cols(:)
        integer :: ncol,io

        debug =.false.
        ierr = 0
        io = alloc_iounit(ierr)
        open(io,file=trim(filename),action='read',status='old',iostat=ierr)
        if(ierr/=0) then
           write(*,*) 'failed to open columns list file: ', trim(filename)
           call free_iounit(io)
           return
        endif
        ncols=0
        do pass=1,2
           if(pass==2) allocate(cols(ncols(1)))
           inner_loop: do while(.true.)
              is_int = .false.
              read(io,'(a)',iostat=ierr) line
              if(ierr/=0) exit inner_loop

              !remove any nasty tabs
              do while(index(line,char(9))>0)
                 i=index(line,char(9))
                 line(i:i)=' '
              enddo

              nchar=len_trim(line)
              if(nchar==0) cycle inner_loop ! ignore blank line

              line=adjustl(line)
              i=index(line,'!')-1
              if(i<0) i=len_trim(line)

              if(i==0) then       !comment line
                 if(debug) write(*,*) ' comment: ', trim(line)
                 cycle inner_loop
              else if(index(line(1:i),'number')>0) then
                 if(debug) write(*,*) '****** ', trim(line)
                 if(debug) write(*,*) '****** index of number > 0 => integer'
                 is_int = .true.
              else if(index(line(1:i),'num_')==1)then
                 if(debug) write(*,*) '****** ', trim(line)
                 if(debug) write(*,*) '****** index of num_ == 1 => integer'
                 is_int = .true.
              endif

              column_name = line
              ncols(pass)=ncols(pass)+1
              if(i==0) then
                 column_length=len_trim(column_name)
              else
                 column_length=len_trim(column_name(1:i))
              endif
              do i=1,column_length
                 if(column_name(i:i)==' ') column_name(i:i)='_'
              enddo
              !if(debug) write(*,'(2i5,a32,i5)') pass, ncols(pass),trim(column_name(1:column_length)), column_length
              if(pass==2) then
                 cols(ncols(pass))% name = trim(column_name(1:column_length))
                 if(is_int) then
                    cols(ncols(pass))% type = column_int
                 else
                    cols(ncols(pass))% type = column_dbl
                 endif
                 cols(ncols(pass))% loc = ncols(pass)
              endif
           end do inner_loop
           if(pass==1) rewind(io)
           if(pass==2) close(io)
        end do
        if(ncols(1)==ncols(2)) then
           ierr=0
           ncol=ncols(1)
        endif
        call free_iounit(io)
        if(debug) write(*,*) 'process_columns: ncol = ', ncol

      end subroutine process_columns


    subroutine read_key_eeps()
    integer :: temp(15), neep,ieep
        temp = -1

        ieep = 1
        temp(ieep) = PreMS_EEP
        if(add_eep(temp,ieep)) ieep=ieep+1

        temp(ieep) = ZAMS_EEP
        if(add_eep(temp,ieep)) ieep=ieep+1

        temp(ieep) = IAMS_EEP
        if(add_eep(temp,ieep)) ieep=ieep+1

        temp(ieep) = TAMS_EEP
        if(add_eep(temp,ieep)) ieep=ieep+1

        temp(ieep) = BGB_EEP
        if(add_eep(temp,ieep)) ieep=ieep+1

        temp(ieep) = cHeIgnition_EEP
        if(add_eep(temp,ieep)) ieep=ieep+1

        temp(ieep) = cHeBurn_EEP
        if(add_eep(temp,ieep)) ieep=ieep+1

        temp(ieep) = TA_cHeB_EEP
        if(add_eep(temp,ieep)) ieep=ieep+1

        temp(ieep) = cCBurn_EEP
        if(add_eep(temp,ieep)) ieep=ieep+1

        temp(ieep) = TPAGB_EEP
        if(add_eep(temp,ieep)) ieep=ieep+1

        temp(ieep) = post_AGB_EEP
        if(add_eep(temp,ieep)) ieep=ieep+1

        temp(ieep) = WD_EEP
        if(add_eep(temp,ieep)) ieep=ieep+1

        temp(ieep) = Extra_EEP1
        if(add_eep(temp,ieep)) ieep=ieep+1

        temp(ieep) = Extra_EEP2
        if(add_eep(temp,ieep)) ieep=ieep+1

        temp(ieep) = Extra_EEP3
        if(.not. add_eep(temp,ieep)) temp(ieep) = -1

    neep = count(temp > 0,1)
    allocate(key_eeps(neep))
    key_eeps = pack(temp,temp > 0)
    
    !define initial and final eep if not already defined
    if (Initial_EEP <0 .or. Initial_EEP< minval(key_eeps))  Initial_EEP = ZAMS_EEP
    if (Final_EEP < 0 .or. Final_EEP > maxval(key_eeps))  Final_EEP = maxval(key_eeps)
    
    end subroutine

    logical function add_eep(temp, i)
    integer :: last, i, temp(:)
        last = 0
        add_eep = .false.
        if (i>1 ) last = temp(i-1)
        if (temp(i) >0 .and. temp(i)<= Final_EEP .and. temp(i)/= last) then
            add_eep = .true.
        endif
    end function

    
  subroutine set_star_type_from_history(x)
    type(eep_track), intent(inout) :: x
    integer :: n

    !set the WDCS primary EEP if center_gamma < center_gamma_limit
    !center_gamma_limit = 19

    !set the CarbonBurn primary EEP if center_c12 < center_carbon_limit
    center_carbon_limit = 0.05

    !set star_type to high_mass_star if max(log_center_T) > this
    log_center_T_limit = 8.7 !changing from 8.5

    !set star_type to high mass star if M_init >= this
    high_mass_limit = 10.0 !Msun

    !from Pols et al. 1998- set star_type to high mass star if He core mass>= this
    he_core_mass_limit = 2.2 !Msun

    n = x% ntrack
    
    !- can be based on co mass in the end<1.4 maybe
    
    !only reach center_gamma_limit if the star evolves to a WD
    !if( x% tr(i_gamma,n) > center_gamma_limit) then
       !x% star_type = star_low_mass
       !return
    !endif

    x% star_type = star_low_mass

    !last gasp test for high-mass stars is the initial mass...
    if(x% initial_mass >= high_mass_limit) then
       x% star_type = star_high_mass
       return
    endif
    
    !i_he_core
    if (x% tr(i_he_core,n)>= he_core_mass_limit) then
        x% star_type = star_high_mass
        return
    endif
    
    !simple test for high-mass stars is that central C is depleted
    if (i_c12 >0) then
        if(maxval(x% tr(i_c12,:)) > 0.4d0 .and. x% tr(i_c12,n) < center_carbon_limit)then
           x% star_type = star_high_mass
           return
        endif
    endif
    
    !alternative test for high-mass stars is that they reach a
    !central temperature threshhold
    if (i_Tc >0) then
        if(x% tr(i_Tc,n) > log_center_T_limit)then
            x% star_type = star_high_mass
            return
        endif
    endif

  end subroutine set_star_type_from_history

    logical function check_mass_loss(x)
    type(eep_track), intent(in) :: x
    real(dp) :: dm

    dm = (x% tr(i_mass,1)-x% tr(i_mass,x% ntrack))/x% initial_mass
    if (dm<tiny) then
        check_mass_loss = .false.
    else
        check_mass_loss = .true.
    endif
    end function

    subroutine set_star_type_from_label(label,s)
        character(LEN=10), intent(in) :: label
        type(eep_track), intent(inout) :: s
        integer :: n,i
        n = size(star_label)
        do i=1,n
            if(label==star_label(i)) s% star_type = i
        enddo
    end subroutine set_star_type_from_label

    subroutine set_zparameters(zpars)
        real(dp), intent(out) :: zpars(20)
        real(dp) :: old_co_frac,co_fraction,change_frac
        real(dp) :: smass,Teff,last_val,he_diff
        real(dp), allocatable :: T_centre(:)
        integer :: len_track, i, min_index
        integer:: j_bagb, j_tagb, i_start
        real(dp), allocatable :: mass_list(:)

        logical:: debug

        debug = .false.
        
        old_co_frac = 0.d0
        Mup_core = 0.d0
        Mec_core = 0.d0

!        allocate(Mcrit(9))
        Mcrit% mass= -1.d0
        Mcrit% loc = 0

        Mcrit(1)% mass = s(1)% initial_mass
        Mcrit(1)% loc = 1

        Mcrit(2)% mass = very_low_mass_limit
        Mcrit(3)% mass = Mhook
        Mcrit(4)% mass = Mhef
        Mcrit(5)% mass = Mfgb
        Mcrit(6)% mass = Mup
        Mcrit(7)% mass = Mec
        Mcrit(8)% mass = Mextra

        Mcrit(9)% mass = s(num_tracks)% initial_mass
        Mcrit(9)% loc = num_tracks+1 !TODO: explain why+1?

        if (verbose) write(*,'(a,f7.1)') 'Minimum initial mass', Mcrit(1)% mass
        if (verbose) write(*,'(a,f7.1)') 'Maximum initial mass', Mcrit(9)% mass
        
        if (.not. defined(Mcrit(7)% mass)) then
          do i = 1,size(s)
            call set_star_type_from_history(s(i))
!            print*, s(i)% initial_mass, s(i)% star_type
          end do
        endif

        allocate(mass_list(num_tracks))
        mass_list = s% initial_mass

        !if already defined, do index search here otherwise search below
        do i = 2, size(Mcrit)-1
            if (.not. defined(Mcrit(i)% mass)) cycle
            call index_search (num_tracks, mass_list, Mcrit(i)% mass, min_index)
            Mcrit(i)% mass = s(min_index)% initial_mass
            Mcrit(i)% loc = min_index
!            if (debug) print*, i, Mcrit(i)% mass
        end do

        i_start = max(Mcrit(1)% loc, Mcrit(2)% loc)
        
        do i = i_start, num_tracks
            smass = s(i)% initial_mass
            !print*,smass, s(i)% star_type
            len_track = s(i)% ntrack
            if (smass<=3.0) then
                !determining Mhook
                !where the maximum of the central temperature (Tc) between
                !the IAMS and the TAMS EEPs is greater than the Tc at
                !the TAMS EEP i.e., Tc,max>Tc,TAMS

                if (.not. defined(Mcrit(3)% mass))then
                    if (len_track >= TAMS_EEP) then
                    !T_centre = s(i)%tr(i_Tc,IAMS_EEP:TAMS_EEP)
                    allocate(T_centre,source=s(i)% tr(i_Tc,IAMS_EEP:TAMS_EEP))

                    last_val = T_centre(size(T_centre))
                    if (maxval(T_centre)>last_val) then
                        Mcrit(3)% mass = smass
                        Mcrit(3)% loc = i
                        if (debug) print*,"Mhook",smass,i
                    endif
                    deallocate(T_centre)

                    endif
                endif

                !determining Mhef
                !the minimum temperature for core helium burning is about 100 million Kelvin.
                !In stars that undergo the helium flash, a slight expansion of the core
                ! following the flash causes the central temperature to decrease
                !a little before increasing again with stable helium burning.
                
                if (.not. defined(Mcrit(4)% mass))then
                if (len_track>=TA_cHeB_EEP) then
                    allocate(T_centre,source=s(i)% tr(i_Tc,cHeIgnition_EEP:TA_cHeB_EEP-1))
                    !T_centre = s(i)% tr(i_Tc,cHeIgnition_EEP:TA_cHeB_EEP-1)
                    if (minval(T_centre)>7.4) then
                        Mcrit(4)% mass = smass
                        Mcrit(4)% loc = i
                        if (debug) print*,"Mhef",smass,i
                    endif
                    deallocate(T_centre)
                endif
                endif

            else        !if (smass>=3.0) then
                !determining Mup
                !where the absolute fractional change in the
                !central carbon-oxygen mass fraction exceeds
                !0.01 at the end of the AGB (TPAGB EEP)
                if ((.not. defined(Mcrit(6)% mass)) .and. smass<8.0) then
                    if (i_c12>0 .and. i_o16 >0) then
                        j_tagb = min(cCBurn_EEP,TPAGB_EEP)      !end of agb
                        j_tagb = min(len_track,j_tagb)
                        co_fraction = s(i)% tr(i_c12,j_tagb)+s(i)% tr(i_o16,j_tagb)
                        if (old_co_frac>0.0) then
                            change_frac = abs(co_fraction-old_co_frac)
                            change_frac = change_frac/old_co_frac
                            if (change_frac>0.01) then
                                ! this is the mass at which C/O ignition occur
                                ! we need the mass preceeding it
                                ! (Mup = M below which C/O ignition doesn't occur)
                                Mcrit(6)% loc = i-1
                                Mcrit(6)% mass = s(Mcrit(6)% loc)% initial_mass
                                if (debug) print*,"Mup",Mcrit(6)% mass, Mcrit(6)% loc
                            endif
                        endif
                    old_co_frac = co_fraction
                    endif
                endif
            endif

            !determining Mfgb- all masses
            if (.not. defined(Mcrit(5)% mass))then
                if (smass<=20.0 .and. len_track>=cHeIgnition_EEP) then
                    Teff = s(i)% tr(i_logTe,cHeIgnition_EEP-1)       !temp at the end of HG/FGB
                    he_diff = abs(s(i)% tr(i_he4, cHeIgnition_EEP-1)-s(i)% tr(i_he4, TAMS_EEP))
!                    print*,"bgb",smass,Teff, he_diff
                    if (Teff> T_bgb_limit .or. he_diff >0.01) then !
                        Mcrit(5)% mass = smass
                        Mcrit(5)% loc = i
                        if (debug) print*,"Mfgb",smass,i
                    endif
                endif
            endif

            !determining Mec
            if (.not. defined(Mcrit(7)% mass))then
                if (s(i)% star_type == star_high_mass) then
                Mcrit(7)% mass = smass
                Mcrit(7)% loc = i
                if (debug) print*,"Mec",smass,i
                endif
            endif

        end do

        !If the tracks are beyond the zpars limits, above procedure
        !picks up the first or second track, which can lead to errors later,
        !hence those values need to be reverted

        do i = 2,size(Mcrit)-1
            if (Mcrit(i)% mass <= Mcrit(1)% mass) then
                Mcrit(i)% mass= -1.d0
                Mcrit(i)% loc = 0
            endif
        end do

        Mcrit(7)% loc = max(Mcrit(7)% loc,1)
        j_bagb = min(s(Mcrit(7)% loc)% ntrack,TA_cHeB_EEP)
        Mec_core = s(Mcrit(7)% loc)% tr(i_he_core,j_bagb)
        
        !if cannot locate Mup or located it beyond Mec (which is incorrect),

        if (Mcrit(6)% loc < 1 .or. Mcrit(6)% loc >=  Mcrit(7)% loc) then
            if (debug) print*, 'Mcrit(6)/Mup not found or is incorrect (exceeds Mcrit(7)/Mec)'
            if (debug) print*, 'Using value closest to Mec-1.8 (the SSE way)'
            !modify Mup by SSE's way
            Mcrit(6)% mass = Mcrit(7)% mass - 1.8d0
            call index_search (num_tracks, mass_list, Mcrit(6)% mass, Mcrit(6)% loc)
            !make sure the new location for Mup does not exceed Mec
            Mcrit(6)% loc = max(1,min(Mcrit(6)% loc,Mcrit(7)% loc-1))
            Mcrit(6)% mass = s(Mcrit(6)% loc)% initial_mass
            if (debug) print*,"new Mup",Mcrit(6)% mass, Mcrit(6)% loc
        endif


        j_bagb = min(s(Mcrit(6)% loc)% ntrack,TA_cHeB_EEP)
        Mup_core = s(Mcrit(6)% loc)% tr(i_he_core,j_bagb)

        if (debug) print*,"Mup_core =", Mup_core
        if (debug) print*,"Mec_core =", Mec_core


        call sort_mcutoff()
        if (debug) print*, "m_cutoffs: ", m_cutoff
    
        !now redefine zpars where applicable
        do i = 3,7
        if (defined (Mcrit(i)% mass)) zpars (i-2) = Mcrit(i)% mass
        end do

        !Redefine these for use later in the code
        Mhook = zpars(1)
        Mhef = zpars(2)
        Mfgb = zpars(3)
        Mup = zpars(4)
        Mec = zpars(5)
    
    
        if (defined(Z_H)) zpars(11) = Z_H
        if (defined(Z_He)) zpars(12) = Z_He
        Z04 = zpars(14)

        if (debug) print*, 'zpars:',  zpars(1:5)

        deallocate(mass_list)
!        call sort(Mcrit% loc, m_cutoff)
    end subroutine set_zparameters

    subroutine sort_mcutoff()
     !subroutine to sort Mcutoffs, removing the ones who are at less than 2 distance from the last one
    !making sure there are at least 2 tracks between subsequent mcutoff
    
!        integer:: mloc(:)
        integer, allocatable :: mloc(:)
        integer :: val,n
        integer :: i, k, loc,a
        
        n = size(Mcrit)
        allocate (m_cutoff(n),mloc(n))
        mloc = Mcrit% loc
        m_cutoff = 0
        m_cutoff(1) = 1
        k=1
        
        do i = 1, n
            val = minval(mloc(i:n))
            a = minloc(mloc(i:n),dim=1)
            loc = (i - 1) + a
            mloc(loc) = mloc(i)
            mloc(i) = val
            if ((val - m_cutoff(k))>=2) then
                k=k+1
                m_cutoff(k) = val
            endif
        end do
        
        !   deallocate(mloc)
        !        allocate(mloc(size(m_cutoff)))
        !        mloc = m_cutoff
        
        m_cutoff = pack(m_cutoff,mask = m_cutoff .ne. 0)
        deallocate(mloc)
    end subroutine sort_mcutoff

    subroutine sort(mloc, list)
        integer, intent(in) :: mloc(:)
        integer, allocatable, intent(out) :: list(:)
        integer :: val,n
        integer :: i, a, loc

        n = size(mloc)
        allocate(list(n))
        list=pack(mloc,mask = mloc .ne. 0)
        
        !sort array
        do i = 0, n-1
            val = minval(list(i:n-1))
            a = minloc(list(i:n-1),dim=1)
            loc = (i - 1) + a
            list(loc) = list(i)
            list(i) = val
        end do
    end subroutine sort

    !ZPARS
    !finds critical masses and their locations
    !1; M below which hook doesn't appear on MS, Mhook. 3
    !2; M above which He ignition occurs non-degenerately, Mhef. 4
    !3; M above which He ignition occurs on the HG, Mfgb. 5
    !4; M below which C/O ignition doesn't occur, Mup. 6
    !5; M above which C ignites in the centre, Mec. 7

    subroutine calculate_sse_zpars(z,zpars)

    real(dp),intent(in) :: z
    real(dp),intent(out) :: zpars(14)
    real(dp) :: lzs,dlzs,lz,lzd

        lzs = log10(z/0.02d0)
        dlzs = 1.d0/(z*log(10.d0))
        lz = log10(z)
        lzd = lzs + 1.d0

        zpars(1) = 1.0185d0 + lzs*(0.16015d0 + lzs*0.0892d0)
        zpars(2) = 1.995d0 + lzs*(0.25d0 + lzs*0.087d0)
        zpars(3) = 16.5d0*z**0.06d0/(1.d0 + (1.0d-04/z)**1.27d0)
        zpars(4) = MAX(6.11044d0 + 1.02167d0*lzs, 5.d0)
        zpars(5) = zpars(4) + 1.8d0
        zpars(6) = 5.37d0 + lzs*0.135d0
        !* set the hydrogen and helium abundances
        zpars(11) = 0.76d0 - 3.d0*z
        zpars(12) = 0.24d0 + 2.d0*z
        !* set constant for low-mass CHeB stars
        zpars(14) = z**0.4d0

    end subroutine
    
    elemental function relative_diff(z1,z2) result(y)
        real(dp),intent(in) :: z1,z2
        real(dp) :: y
        !the formula can catch the difference in small values of metallicity
        !while allowing for differences due to precision errors
        y = abs(z1-z2)/MIN(z1,z2)
    end function

end module z_support
