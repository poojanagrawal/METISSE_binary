subroutine evolv1(mass,max_age,ierr)

    !subroutine to evolve a single star and write output to file
    !calls star.f90 for interpolation in mass and calls hrdiag.f90 for interpolation in age
    !timesteps are determined through deltat.f90

    use track_support
    use remnant_support, only:time_He_MS

    implicit none

    !variable declaration
    real(dp), intent(in):: mass,max_age
    integer, intent(out):: ierr
    integer :: id

    !type(track), target :: t
    integer :: str,lines,old_phase
    real(dp):: tphys,timestep
    character(len=strlen) :: output_file,eep_filename
    logical :: output, t_end
    real(dp) :: dms, M_env
   type(track), pointer :: t => NULL()
    
    t=> tarr(id)
   
    ierr = 0
    call star(id,mass)

    if (t% complete .eqv..false.) ierr = 1
    str=int(mass*100)

  !write the mass interpolated track if write_eep_file is true
    if (write_eep_file) then
        write(eep_filename,"(a,a,i5.5,a)") trim(METISSE_DIR),"/output/",str,"M_full.track.eep"
        call write_eep_track(eep_filename,t)
    end if

    lines = 0
    tphys = 0.0
    timestep = 0.0
    old_phase = -1
    t_end  = .false.
    t% pars% phase = t% phase(initial_eep)

    !SSE like output file if write_track_to_file is true
    
!TODO: update the output file location is no longer valid:remove it
    output = write_track_to_file
    if (output) write (output_file,"(a,a,i5.5,a)") trim(METISSE_DIR), "/output/evolve_", str, "M.dat"
    if (output) open (120,FILE=trim(output_file),action="write")

    if (output) write(120,'(9a15,2a10)') "time", "age", "mass","core_mass","He_core" &
                ,"CO_core","log_L","log_Teff","log_radius", "phase","e"

    do while(.true.)
        if (tphys>=max_age) then
            tphys = max_age
            t_end  = .true.
            !flag for telling the loop that this is the last timestep, so save
            !it is different to end_of_file defined in hrdiag
        end if
        
        !evolve the star- calculate stellar parameters at tphys
        t% pars% age = tphys - t% pars% epoch
        call hrdiag(id,tphys)
        
         if (t% pars% phase>6 .or. t% post_agb) then
            t% pars% log_L = log10(t% pars% luminosity)
            t% pars% Teff = 1000*((1130.d0*t% pars% luminosity/(t% pars% radius**2))**0.25)
            t% pars% log_Teff = log10(t% pars% Teff)
            t% pars% log_R =  log10(t% pars% radius)
        endif
       
        !write output if flags are true
        if (old_phase /=t% pars% phase .or. t_end) then
            if (verbose) write(*,'(a10,f10.1,3a10,f7.3)') "Time ", tphys, "Phase ", phase_label(t% pars% phase+1), &
                                                                                "Mass ", t% pars% mass
            if (output) call write_dat_track(tphys, t% pars)
        else if (lines<3000 .and. output) then
            call write_dat_track(tphys,t% pars)
             !if (mod(i,1)==0 .or. dtp intervals)then  !***TODO
        endif
    
        lines=lines+1
        old_phase = t% pars% phase
        if (t_end) exit

        !calculate next time step
        call deltat(id,kw,age,tm,tn,tscls,dt,dtr)
        timestep = min(dt,dtr)
         
         !Calculate mass loss and modify timestep if needed  -- only for He stars.
        if (t% pars% phase>=He_MS .and. t% pars% phase<=He_GB) then
        
            ! Calculate mass loss for the timestep.
            dms = mlwind(t% pars,initial_z)
           dms = dms *1.0d+06*timestep
           M_env = pars% mass - pars% core_mass
            if(dms>=M_env)then
              timestep = (M_env/dms)*timestep
                !print*,"timestep modified1",timestep,dms
              dms = M_env
           endif
        !Limit to 1% mass loss.
        if(dms.gt.0.01d0*pars% mass)then
            timestep = 0.01d0*timestep*pars% mass/dms
            !print*,"timestep modified2",timestep,dms

            dms = 0.01d0*pars% mass
        endif
        pars% mass = pars% mass - dms
        !print*,"in mlwind, dms",dms, pars% mass,timestep,pars% phase
        return
            if (t% pars% phase == He_MS) then
                t% zams_mass = t% pars% mass
                t% pars% epoch= tphys-(t% pars% age *time_He_MS(t% pars% mass)/t% ms_time)
            endif
        endif
        
        tphys = tphys+timestep
    end do
    
    if (output) close(120)

   call dealloc_track(id)
    nullify(t)
    if (verbose) write(*,*) "-------------------------------------------------------------------------"
    return
    end subroutine evolv1

      
