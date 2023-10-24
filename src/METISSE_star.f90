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

    real(dp) :: times_old(11), nuc_old, delta,dtm, delta1,delta_wind

    integer :: idd, nt ,ierr
    logical :: debug, interpolate_all,consvR
    type(track), pointer :: t

    idd = 1
    if(present(id)) idd = id
    t => tarr(idd)
        
    ierr=0
    
    debug = .false.
!    if ((id == 1) .and. (kw<=7))debug = .true.

    if (debug) print*, '-----------STAR---------------'
    if (debug) print*,"in star", mass,mt,kw,id

    consvR = .false.
    delta = 0.d0
    t% zams_mass = mass

    select case(kw)
        case(low_mass_MS:TPAGB)
            if (t% pars% mass < 0) then
                ! Initial call- just do normal interpolation
                if (debug) print*, "normal interpolate mass", t% initial_mass,t% zams_mass,t% pars% mass,mt,kw
                t% initial_mass = mass
                call interpolate_mass(t% initial_mass,id)
                call calculate_timescales(t)

                !Todo: explain what age 2 and Times_new are
                t% times_new = t% times
                t% tr(i_age,:) = t% tr(i_age2,:)
                t% ms_old = t% times(MS)
                mt = t% tr(i_mass,ZAMS_EEP)
        !            call write_eep_track(t,t% initial_mass)
            else
                !first check whether star lost its envelope during binary interaction
                !to avoid unneccesssary call to interpolation routine
                if (t% pars% core_mass.ge.mt) then
                    if (debug) print*, 'star has lost envelope, exiting star',t% pars% core_mass,mt
                else
                    ! Next check if mass has changed since last time star was called.
                    ! For tracks that already have wind mass loss,
                    ! exclude mass loss due to winds
                    if (t% has_mass_loss) then
                        delta_wind = (t% pars% dms*abs(dtm)*1.0d+06)
                    else
                        delta_wind = 0.d0
                    endif
                    
                    !this is to avoid extra mass loss if star subroutine
                    ! is called multiple times in same step
                    delta = (t% pars% mass-mt) - delta_wind
                    delta = delta -t% pars% delta
                    delta1 = 1.0d-04*mt
                    if (abs(delta) .ge. delta1) then
                
                        if (debug) print*, "mass loss in interpolate mass called for", &
                                                        t% initial_mass,mt,delta,id
                        t% times_new = -1.0
                        nt = t% ntrack
                        call get_initial_mass_for_new_track(t,delta,interpolate_all,idd)
                        if (debug)print*, 'initial mass for the new track',t% initial_mass
                        t% ms_old = t% times(MS)
                        if (interpolate_all) then
                            ! kw=0,1: main-sequence star, rewrite all columns with new track
                            if (debug)print*, 'main-sequence star, rewrite with new track'
                            call interpolate_mass(t% initial_mass,id)
                            if (t% ntrack<nt) print*, '***WARNING: track length reduced***',t% initial_mass,nt,t% ntrack
                            
                            ! Calculate timescales and assign SSE phases (Hurley et al.2000)
                            call calculate_timescales(t)
                            t% times_new = t% times
                            t% tr(i_age,:) = t% tr(i_age2,:)
                        else
                            !store core properties for post-main sequence evolution
                            if (debug)print*, 'post-main-sequence star'
                            allocate(hecorelist(nt),ccorelist(nt),Lum_list(nt),age_list(nt))
                            hecorelist = t% tr(i_he_core,:)
                            ccorelist =  t% tr(i_co_core,:)
                            Lum_list = t% tr(i_logL,:)

                            times_old = t% times
                            nuc_old = t% nuc_time
                            age_list = t% tr(i_age,:)
                            rlist = t% tr(i_logR,:)
                            if ((t% pars% mcenv/t% pars% mass).ge.0.2) consvR= .true.
    !                            .and.(t% pars% env_frac.ge.0.2)
                            call interpolate_mass(t% initial_mass,id)
                           !write the mass interpolated track if write_eep_file is true
                            if (kw>=1 .and. kw<=4 .and. .false.) then
                                call write_eep_track(t,mt)
                            end if
                            
                            if (t% ntrack<nt) print*, '**WARNING: track length reduced**',t% initial_mass,nt,t% ntrack

                            call calculate_timescales(t)
                            t% times_new = t% times
                            t% times = times_old
                            
                            t% nuc_time = nuc_old

                            t% tr(i_age,:) = age_list(1:t% ntrack)
                            t% tr(i_he_core,:) = hecorelist(1:t% ntrack)
                            t% tr(i_co_core,:) = ccorelist(1:t% ntrack)
                            t% tr(i_logL,:) = Lum_list(1:t% ntrack)
                            !TEST: R remains unchanged if Mconv is significant
                            if (consvR) then
                                t% tr(i_logR,:) = rlist(1:t% ntrack)
                            endif
                            deallocate(hecorelist,ccorelist,Lum_list,age_list,rlist)
                        endif
                        t% pars% delta = delta
                    endif
                endif
            endif
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

    if (debug) print*, "in star end", mt,delta,kw,tm,tscls(1),tn,t%initial_mass
    nullify(t)
    return
end subroutine METISSE_star

