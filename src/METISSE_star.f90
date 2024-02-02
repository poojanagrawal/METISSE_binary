subroutine METISSE_star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars,dtm,id)

    use track_support
    use interp_support
    use sse_support
    implicit none

    integer, intent(in), optional :: id
    real(dp) :: mass,mt,tm,tn,tscls(20),lums(10),GB(10),zpars(20)
    integer :: kw

!    real(dp), allocatable:: hecorelist(:),ccorelist(:),Lum_list(:),age_list(:)
    real(dp), allocatable:: rlist(:)

    real(dp) :: times_old(11), nuc_old, delta,dtm, delta1,delta_wind, quant

    integer :: idd, ierr,  age_col
    logical :: debug, exclude_core,consvR
    type(track), pointer :: t

    idd = 1
    if(present(id)) idd = id
    t => tarr(idd)
        
    ierr = 0
    
    debug = .false.
!    if ((id == 1) .and. (kw==8))debug = .true.

    if (debug) print*, '-----------STAR---------------'
    if (debug) print*,"in star", mass,mt,kw,id,t% pars% phase

    consvR = .false.
    exclude_core = .false.
    delta = 0.d0
    t% zams_mass = mass

    if(kw>= He_MS .and. kw<=He_GB .and. use_sse_NHe) t% star_type = sse_he_star
    if(kw>= HeWD) t% star_type = remnant       !to be double sure, in case star becomes a remnant outside hrdiag
    
    select case(t% star_type)
    case(unknown)
        ! Initial call- just do normal interpolation
        t% initial_mass = mass
        
        if (debug) print*, "initial interpolate mass", t% initial_mass,t% zams_mass,t% pars% mass,mt,kw
        
        t% pars% delta = 0.d0
        t% pars% extra = 0
        t% pars% bhspin = 0.d0
        t% is_he_track = .false.
        t% ierr = 0
        call interpolate_mass(id,exclude_core)
        call calculate_timescales(t)

        !Todo: explain what age 2 and Times_new are
        t% times_new = t% times
        t% tr(i_age,:) = t% tr(i_age2,:)
        mt = t% tr(i_mass,ZAMS_EEP)
        t% ms_old = t% times(MS)
        t% pars% age_old = t% pars% age
        t% initial_mass_old = t% initial_mass
        
        !call write_eep_track(t,t% initial_mass)

    case(rejuvenated) !gntage
        t% pars% mass = mt
        t% pars% delta = 0.d0
        t% pars% phase = kw
        t% times_new = -1.d0
        if(kw>= He_MS .and. kw<=He_GB) then
            t% is_he_track = .true.
            age_col = i_he_age
        else
            t% is_he_track = .false.
            age_col = i_age
        endif
    
        call get_initial_mass_for_new_track(t,id)
        if (debug)print*, 'rejuvanted giant, rewrite with new track'
        call interpolate_mass(id,exclude_core)
        
        ! Calculate timescales and assign SSE phases (Hurley et al.2000)
        call calculate_timescales(t)
        t% times_new = t% times
        t% tr(age_col,:) = t% tr(i_age2,:)
        mass = t% initial_mass
        t% ms_old = t% times(MS)
        t% pars% age_old = t% pars% age
        t% initial_mass_old = t% initial_mass
        t% zams_mass = mass
    case(remnant)
        tm = 1.0d+10
        tscls(1) = t% MS_time
        tn = 1.0d+10
        if (debug) print*, 'remnant'
    case(sse_he_star)
         call calculate_SSE_He_timescales(t)
         call calculate_SSE_He_star(t,tscls,lums,GB,tm,tn)
    case(sub_stellar:star_high_mass)
        if (debug) print*,'nuc burn star'
        if (t% pars% core_mass.ge.mt) then
            !Check whether star lost its envelope during binary interaction
            !to avoid unneccesssary call to interpolation routine
!            if (debug)write(UNIT=err_unit,fmt=*)
            if (debug) print*,'star has lost envelope, exiting star',t% pars% core_mass,mt,kw
        else
            ! Check if mass has changed since last time star was called.
            ! For tracks that already have wind mass loss,
            ! exclude mass loss due to winds
            if(kw>= He_MS .and. kw<=He_GB) then
                t% is_he_track = .true.
                age_col = i_he_age
                write(UNIT=err_unit,fmt=*) 'switching to he tracks'
            else
                t% is_he_track = .false.
                age_col = i_age
            endif
            
            if (t% has_mass_loss .and. (abs(mt-t% pars% mass))>1.0d-06) then
                delta_wind = (t% pars% dms*dtm*1.0d+06)
            else
                delta_wind = 0.d0
            endif

            delta = (t% pars% mass-mt) - delta_wind
    !       if (id==1)print*, 'delta org', delta,t% pars% mass,mt,delta_wind
            t% pars% delta = t% pars% delta+ delta
            delta1 = 1.0d-04*mt
            
            if (debug) print*, 'delta is', delta,t% pars% delta,delta1,t% pars% mass
            if (delta.ge.10) then
                write(UNIT=err_unit,fmt=*)'large delta',delta,t% pars% mass,id
    !           call stop_code
            endif
            
            if ((abs(t% pars% delta).ge.delta1).and.(t% post_agb.eqv..false.)) then
                if(debug)print*,"mass loss in interpolate mass called for",t% initial_mass,mt,t% pars% mass,id
                    
                if (kw>MS .and. kw /=He_MS) exclude_core = .true.

    !           if (dtm<0.d0) print*,'dtm<0',t% pars% age,dtm

                quant = (t% pars% age*t% MS_old/t% ms_time)+dtm
                if (dtm<0.d0 .and. (quant .le.t% pars% age_old).and. kw<=MS) then
                    if (debug)print*,'reverse intp',quant,t% pars% age_old
                    t% initial_mass = t% initial_mass_old
                else
                    t% ms_old = t% times(MS)
                    t% pars% age_old = t% pars% age
                    t% initial_mass_old = t% initial_mass
                    call get_initial_mass_for_new_track(t,idd)
                endif
                    
                if (debug)print*, 'initial mass for the new track',t% initial_mass
                
                ! reset delta
                t% pars% delta = 0.d0
                t% times_new = -1.d0
                
                if (exclude_core) then
                    !store core properties for post-main sequence evolution
                    if (debug) print*,'post-main-sequence star'

                    times_old = t% times
                    nuc_old = t% nuc_time
                    
                    if ((t% pars% mcenv/t% pars% mass).ge.0.2) then
!                     .and.(t% pars% env_frac.ge.0.2)
                        allocate(rlist(t% ntrack))
                        rlist = t% tr(i_logR,:)
                        consvR= .true.
                    endif
    !
                    call interpolate_mass(id,exclude_core)
                   !write the mass interpolated track if write_eep_file is true
                    if (kw>=1 .and. kw<=4 .and. .false.) call write_eep_track(t,mt)
                    
                    call calculate_timescales(t)
                    t% times_new = t% times
                    t% times = times_old
                    t% nuc_time = nuc_old

                    !TEST: R remains unchanged if Mcenv is significant
                    if (consvR) then
                        t% tr(i_logR,:) = rlist(1:t% ntrack)
                        deallocate(rlist)
                    endif
                else
                    ! kw=0,1,7: main-sequence star, rewrite all columns with new track
                    if (debug)print*, 'main-sequence star, rewrite with new track'
                    call interpolate_mass(id,exclude_core)
                    
                    ! Calculate timescales and assign SSE phases (Hurley et al.2000)
                    call calculate_timescales(t)
                    t% times_new = t% times
                    t% tr(age_col,:) = t% tr(i_age2,:)
                endif
            endif
            t% pars% mass = mt !+ delta_wind
            ! for certain phases star is called twice to calculate epoch
            ! above is to prevent multiple mass calculations in the same step
            ! until hrdiag is called and t% pars% mass gets re-evaluated there
            endif
    end select
    
    if (kw<= TPAGB) call calculate_SSE_parameters(t,zpars,tscls,lums,GB,tm,tn)
    
    t% MS_time = tm
    t% nuc_time = tn

!    if (debug) print*, "in star end", mt,delta,kw,tm,tscls(1),tn,t%initial_mass
    nullify(t)
    return
end subroutine METISSE_star

