subroutine METISSE_star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars,dtm,id)

    use track_support
    use interp_support
    use sse_support
    implicit none

    integer, intent(in), optional :: id
    real(dp) :: mass,mt,tm,tn,tscls(20),lums(10),GB(10),zpars(20)
    integer :: kw

    real(dp), allocatable:: hecorelist(:),ccorelist(:),Lum_list(:),age_list(:)
    real(dp), allocatable:: rlist(:)

    real(dp) :: times_old(11), nuc_old, delta,dtm! ,tnext,mnext

    integer :: idd, nt ,ierr
    logical :: debug, interpolate_all,consvR
    type(track), pointer :: t

    idd = 1
    if(present(id)) idd = id
    t => tarr(idd)
        
    ierr=0
    consvR = .false.
    debug = .false.
!    if ((id == 1) .and. (kw>=6))debug = .true.

    if (debug) print*, '-----------STAR---------------'
    if (debug) print*, "in star", mass,mt,kw,id


    delta = 0.d0
    t% zams_mass = mass

    select case(kw)
        case(low_mass_MS:EAGB)
            if (t% pars% mass < 0) then
                ! Initial call- just do normal interpolation
                if (debug) print*, "normal interpolate mass", t% initial_mass,t% zams_mass,t% pars% mass,mt,kw
                t% initial_mass = mass
                call interpolate_mass(t% initial_mass,idd)
                call calculate_phases_and_times(t)

                !Todo: explain what age 2 and Times_new are
                t% times_new = t% times
                t% tr(i_age,:) = t% tr(i_age2,:)
                t% ms_old = t% times(MS)
                mt = t% tr(i_mass,ZAMS_EEP)
        !            call write_eep_track(t,t% initial_mass)
            else
                !first check if mass has changed since last time star was called
                ! for tracks that already have wind mass loss, only check
                ! for binary mass loss
                if (t% has_mass_loss) then
                    delta = (t% pars% mass-mt)-(t% pars% dms*dtm *1.0d+06)
                else
                    delta = (t% pars% mass-mt)
                endif
                delta = delta -t% pars% delta
                if (debug) print*,'delta in star',kw, t% pars% mass, mt,delta

                !next check whether star lost its envelope during binary interaction
                !to avoid unneccesssary call to interpolation routine
                if (t% pars% core_mass.ge.mt) then
                    if (debug) print*, 'star has lost envelope, exiting star',t% pars% core_mass,mt
                elseif (abs(delta) .gt. 1.0d-6) then
                    if (debug) print*, "mass loss in interpolate mass called for", &
                                                    t% initial_mass,delta,t% tr(i_mass,1),id
                    t% times_new = -1.0
                    nt = t% ntrack
                    call get_initial_mass_for_new_track(idd, delta,interpolate_all)
                    if (debug)print*, 'initial mass for the new track',t% initial_mass
                    t% ms_old = t% times(MS)

                    if (interpolate_all) then
                        ! kw=0,1: main-sequence star, rewrite all columns with new track
                        if (debug)print*, 'main-sequence star, rewrite with new track'
                        call interpolate_mass(t% initial_mass,idd)
                        if (t% ntrack<nt) print*, '***WARNING: track length reduced***',t% initial_mass,nt,t% ntrack
                        
                        ! Calculate timescales and assign SSE phases (Hurley et al.2000)
                        call calculate_phases_and_times(t)
                        t% times_new = t% times
                        t% tr(i_age,:) = t% tr(i_age2,:)
                        
                    else
                        !store core properties for post-main sequence evolution
                        if (debug)print*, 'post-main-sequence star'
                        allocate(hecorelist(nt), ccorelist(nt),Lum_list(nt),age_list (nt))
                        hecorelist = t% tr(i_he_core,:)
                        ccorelist =  t% tr(i_co_core,:)
                        Lum_list = t% tr(i_logL,:)

                        times_old = t% times
                        nuc_old = t% nuc_time
                        age_list = t% tr(i_age,:)
                        rlist = t% tr(i_logR,:)
                        if (t% pars% cenv_frac.ge.0.2) consvR= .true.
!                            .and.(t% pars% env_frac.ge.0.2)
                        call interpolate_mass(t% initial_mass,idd)
                       !write the mass interpolated track if write_eep_file is true
                        if (kw>=1 .and. kw<=4 .and. .false.) then
                            call write_eep_track(t,mt)
                        end if
                        
                        if (t% ntrack<nt) print*, '**WARNING: track length reduced**',t% initial_mass,nt,t% ntrack

                        call calculate_phases_and_times(t)

                        t% times_new = t% times
                        t% times = times_old
                        t% nuc_time = nuc_old

                        t% tr(i_age,:) = age_list
                        t% tr(i_he_core,:) = hecorelist
                        t% tr(i_co_core,:) = ccorelist
                        t% tr(i_logL,:) = Lum_list
                        !TEST: R remains unchanged if Mconv is significant
                        if (consvR) then
                            t% tr(i_logR,:) = rlist
!                        print*, 'conserving r', kw
                        endif
                        deallocate(hecorelist,ccorelist,Lum_list,age_list,rlist)
                        
                    endif
                    t% pars% delta = delta
                endif
            endif
            call calculate_SSE_parameters(t,zpars,tscls,lums,GB,tm,tn)
        case(TPAGB)
            call calculate_SSE_parameters(t,zpars,tscls,lums,GB,tm,tn)
        case(He_MS:He_GB)
            call calculate_He_timescales(t)
            call calculate_SSE_He_star(t,tscls,lums,GB,tm,tn)
        case(HeWD:Massless_Rem)
            tm = 1.0d+10
            tscls(1) = t% MS_time
            tn = 1.0d+10
    end select
    t% MS_time = tm
    t% nuc_time = tn

    if (debug) print*, "in star end", mt,delta,kw,tm,tscls(1),tn
    nullify(t)
    return
end subroutine METISSE_star

