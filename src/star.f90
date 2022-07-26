subroutine star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars,dtm,id)

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

    integer :: str, idd, nt
    character(len=strlen) :: eep_filename

    logical :: debug, mass_loss
    type(track), pointer :: t

    t => NULL()
    if(present(id))then
        t => tarr(id)
        idd = id
    else
        t => tarr(1)
        idd = 1
    endif
    mass_loss = .false.
    debug = .false.

    if (kw<10 .and. debug) print*, "in star", mt,kw
    
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
            t% times_new = t% times
            t% tr(i_age,:) = t% tr(i_age2,:)
            mt = t% tr(i_mass,ZAMS_EEP)
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

            if (abs(delta) .gt. 1.0d-8) then
                mass_loss = .true.
!            if (kw<10 .and. debug) print*, "mass loss in interpolate mass called for", t% initial_mass,delta,t% tr(i_mass,1),id

                nt = t% ntrack

                call get_initial_mass_for_new_track(idd, delta)
                if (kw>1) then
                    !store core properties for post-main sequence evolution
                    allocate(hecorelist(nt), ccorelist(nt),Lum_list(nt),age_list (nt))
                    hecorelist = t% tr(i_he_core,:)
                    ccorelist =  t% tr(i_co_core,:)
                    Lum_list = t% tr(i_logL,:)

                    times_old = t% times
                    nuc_old = t% nuc_time
                    age_list = t% tr(i_age,:)
                    rlist = t% tr(i_logR,:)

                    call interpolate_mass(t% initial_mass,idd)
!                    write the mass interpolated track if write_eep_file is true
                    if (kw>=1 .and. kw<=4 .and. .false.) then
                        str = int(mt*100)
                        print*,'writing',str,'M.track.eep'
                        write(eep_filename,"(a,a,i5.5,a)") trim(METISSE_DIR),"/output_eep/",str,"M.track.eep"
                        call write_eep_track(eep_filename,t)
                    end if
                    
                    if (t% ntrack<nt) print*, '***WARNING: track length reduced***',nt,t% ntrack

                    call calculate_phases_and_times(t)
                    !above should use i_age2 and not i_age1
                
                    t% times_new = t% times
                    t% times = times_old
                    t% nuc_time = nuc_old

                    t% tr(i_age,:) = age_list
                    t% tr(i_he_core,:) = hecorelist
                    t% tr(i_co_core,:) = ccorelist
                    t% tr(i_logL,:) = Lum_list
                    !TEST: forcing r to remain unchanged
!                    if (kw>=4 .and. mass<=20) t% tr(i_logR,:) = rlist
                    deallocate(hecorelist,ccorelist,Lum_list,age_list,rlist)
                else
                    ! kw=0,1: main-sequence star, rewrite with new track
                    call interpolate_mass(t% initial_mass,idd)
                    if (t% ntrack<nt) print*, '***WARNING: track length reduced***',nt,t% ntrack
                    !   Calculate timescales and assign SSE phases (Hurley et al.2000)
                    call calculate_phases_and_times(t)
                    t% times_new = t% times
                    t% tr(i_age,:) = t% tr(i_age2,:)
                endif
                !TEST
!                tnext = t% pars% age+ dtm
!                call interpolate_age(t, tnext, i_mass, mnext)
!                print*, "next in star",mnext, mt , mnext-mt
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


    if (kw<10 .and. debug) print*, "in star end", mt,tm,delta,kw
!    print*, "in star", tm, t% MS_time, tn, t% nuc_time, kw
nullify(t)
    return
end subroutine star

