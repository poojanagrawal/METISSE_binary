subroutine initialize_front_end(front_end_name,path_to_metisse)
    use track_support
    character (len=*), intent(in) :: front_end_name
    character(len=*), intent(in), optional :: path_to_metisse

    if (front_end_name == 'main') then
        ! METISSE's main code as described in Agrawal et al. 2020
        ! Can be used to evolve single stars and/or debugging purposes.
        ! Note: currently not functional.
        ! Needs to be updated for latest updates.
        front_end = main
        METISSE_DIR = '.'
        !using relative path here for convenience, should use absolute path
    elseif (front_end_name == 'SSE' .or. front_end_name == 'BSE') then
        ! SSE (Single Star Evolution) from Hurley et al. 2000
        ! BSE (Binary Star Evolution) from Hurley et al. 2002
        front_end = BSE
        METISSE_DIR = '.'
        !using relative path here for convenience, should use absolute path

    elseif (front_end_name == 'COSMIC') then
        ! COSMIC (Compact Object Synthesis and Monte Carlo Investigation Code)
        ! Binary evolution code from Breivik et al 2020
        front_end = COSMIC
        METISSE_DIR = path_to_metisse
        
    else
        print*, "Error: Unrecongnized front_end_name for METISSE"
        print*, "Choose from 'main', 'SSE', 'BSE', 'COSMIC' "
        STOP
    endif
    
end subroutine initialize_front_end
