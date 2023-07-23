 subroutine METISSE_hrdiag(mass,aj,mt,tm,tn,tscls,lums,GB,zpars,&
                                r,lum,kw,mc,rc,menv,renv,k2,mcx,id)
    use track_support
    use interp_support, only: interpolate_age
    use sse_support
    use remnant_support

    implicit none
    integer, intent(in), optional :: id
    real(dp) :: mass,aj,mt,tm,tn,tscls(20),lums(10),GB(10),zpars(20)
    real(dp) :: r,lum,mc,rc,menv,renv,k2,mcx
    
    integer :: kw,i,idd,j_bagb
    real(dp) :: rg,rzams,rtms
    real(dp) :: Mcbagb
    real(dp) :: bhspin ! only for cosmic
    type(star_parameters) :: old_pars

    logical :: debug, post_agb
    logical :: has_become_remnant
    type(track), pointer :: t

    INTEGER irecord
    COMMON /REC/ irecord

    debug = .false.
    idd = 1
    if(present(id)) idd = id
    t => tarr(idd)
    
!    if ((id == 2))debug = .true.

    if (debug) print*, '-----------HRDIAG-------------'
    if (debug) print*,"started hrdiag",mt,mc,aj,kw,tn

    end_of_file = .false. !this is just the end of eep track
    has_become_remnant = .false.
    
    !save input state
    post_agb = t% post_agb
    old_pars = t% pars
  
    t% pars% mass = mt
    t% pars% phase = kw
!    if (kw <=MS) then
!        t% pars% age = aj*(t% MS_time/t% ms_old)
!
!!            them_new = t% times_new(kw)
!!            frac = (input_age)/t% times(ms)*t% times_new(kw)
!!            age2 = (frac*them_new)
!    else
!        t% pars% age = aj
!    endif
     t% pars% age = aj
    t% irecord = irecord

    !if(irecord>0) print*,"In Hrdiag aj,tn ",t% pars% age,mt,t% tr(i_age,t% ntrack),t% tr(i_age2,t% ntrack)
    !print*, 'age, final time',t% pars% age,t% tr(i_age,t% ntrack),abs(t% pars% age-t% tr(i_age,t% ntrack))

    select case(t% pars% phase)
        case(low_mass_MS:TPAGB)
            if (t% post_agb) then
                ! contruct an artficial track until WD cooling phase is reached
                call evolve_after_agb(t)
            elseif (check_ge(t% pars% age,t% tr(i_age,t% ntrack))) then
                !have reached the end of the eep track; self explanatory
                
                if (kw<5 .and.verbose) &
                        print*,'WARNING: possible early end of file due to incomplete track beyond phase',kw,t% initial_mass
                if (debug) print*,"end of file:aj,tn ",t% pars% age,t% tr(i_age,t% ntrack),t% tr(i_age2,t% ntrack)

                end_of_file = .true.
                j_bagb = min(t% ntrack, TA_cHeB_EEP)
                Mcbagb = t% tr(i_he_core, j_bagb)
                if (check_remnant_phase(t% pars, mcbagb)) has_become_remnant = .true.
            else
                !interpolation in age and other checks for nuclear burning phases
                t% pars% core_radius = -1.0
                call interpolate_age(t,t% pars% age)
                if (debug)print*, "mt difference",t% pars% mass, mt, mt - t% pars% mass,kw

                t% pars% mass = mt
                !check if phase/type/kw of the star has changed
                do i = t% pars% phase,6
                    if (i== 0 .or. (.not. defined(t% times(i+1)))) exit
                    if (check_ge(t% pars% age,t% times(i))) then
                        t% pars% phase = i+1
                        if (debug) print*,"phase change",t% pars% age,t% times(i),i+1
                    endif
                end do

                !check if envelope has been lost
                if ((t% pars% core_mass.ge.t% pars% mass) .or. &
                    (t% initial_mass>=10.0 .and. abs(t% pars% core_mass-t% pars% mass)<0.01)) then
                        
                    if (debug) print*, "envelope lost at", t% pars% mass, t% pars% age,t% pars% phase
                    if (t% pars% phase == TPAGB) then
                        ! TPAGB star becomes a WD upon losing envelope
                        ! calling assign_remnant_phase to determine CO-WD/ONe WD
                        j_bagb = min(t% ntrack, TA_cHeB_EEP)
                        Mcbagb = t% tr(i_he_core, j_bagb)
                        has_become_remnant = .true.
                        !TODO: add a check if it's not a white dwarf
                    else
                        call assign_stripped_star_phase(t)
                        call evolve_after_envelope_loss(t)
                        call calculate_SSE_He_star(t,tscls,lums,GB,tm,tn)
                    endif
                endif
            endif
        case(He_MS)
            call evolve_after_envelope_loss(t)
            rzams = t% He_pars% Rzams
        case(He_HG:He_GB)
            call evolve_after_envelope_loss(t)
            rzams = t% He_pars% Rzams
            Mcbagb = t% zams_mass
            if(check_remnant_phase(t% pars, mcbagb)) has_become_remnant = .true.
            
        case(HeWD:Massless_Rem)
            if (front_end == main .or. front_end == BSE) then
                call evolve_remnants_METISSE(t% pars)
            elseif (front_end == COSMIC) then
                call hrdiag_remnant(zpars,t% pars% mass,t% pars% core_mass,t% pars% luminosity,&
                                t% pars% radius,t% pars% age,t% pars% phase)
            endif
            rc = r
            menv = 1.0d-10
            renv = 1.0d-10
            k2 = 0.21d0
        
    end select
      
      
    if(has_become_remnant) then
        if (front_end == main .or. front_end == BSE) then
            call assign_remnant_METISSE(t% pars, mcbagb)
            call post_agb_parameters(t,kw)
        elseif (front_end == COSMIC) then
            call assign_remnant(zpars,t% pars% core_mass,mcbagb,t% zams_mass,&
                                t% pars% mass,t% pars% phase,bhspin,id)
            call post_agb_parameters(t,kw)
            ! maybe instead set mt here
            call hrdiag_remnant(zpars,t% pars% mass,t% pars% core_mass,t% pars% luminosity,&
                                t% pars% radius,t% pars% age,t% pars% phase)
        endif
        has_become_remnant = .false.
        !has_become_remnant is only for assigning remnants, setting it to false now
    endif



    lum = t% pars% luminosity
    r =  t% pars% radius
    mc = t% pars% core_mass
    rc = t% pars% core_radius
    mcx = t% pars% McCO

    mt = t% pars% mass
    kw = t% pars% phase
    aj = t% pars% age
    mass= t% zams_mass


    ! Calculate mass and radius of convective envelope, and envelope gyration radius.
    !this should be separate to above loop as phases may change during the evolution step

    select case(t% pars% phase)
    case(low_mass_MS:TPAGB)
        !rc is calculated during age interpolation if neccessary columns are present,
        !revert to SSE method if those columns are not present
        if (rc<0) CALL calculate_rc(t,tscls,zpars,rc)
        t% pars% core_radius = rc
        if (i_mcenv>0) then
            !calculate mass of convective envelope
            call interpolate_age(t,t% pars% age,i_mcenv,menv)
            menv = min(menv,mt-mc)  ! limit it to the total envelope mass

            if (i_rcenv>0) then
                call interpolate_age(t,t% pars% age,i_rcenv,renv)
            else
                if((mt-mc)>0) then
                    renv = (r - rc)*menv/(mt - mc)
                else
                    renv = 0.d0
                endif
            endif

            !following two lines have been copied from SSE's mrenv.f
            !1.0d-10 can be replaced by 0.d0
            menv = MAX(menv,1.0d-10)
            renv = MAX(renv,1.0d-10)
            k2 = 0.21
            !TODO: add similar lines K2(M.I.) if column is present)
        else
            !revert to SSE method
            rzams = 10.d0**t% tr(i_logR, ZAMS_EEP)
            rtms = 10.d0**t% tr(i_logR, TAMS_EEP)
            CALL calculate_rg(t,rg)
            CALL mrenv(kw,mass,mt,mc,lum,r,rc,aj,tm,lums(2),lums(3),&
            lums(4),rzams,rtms,rg,menv,renv,k2)
        endif

    case(He_MS:He_GB)
        CALL calculate_rc(t,tscls,zpars,rc)
        CALL calculate_rg(t,rg)
        CALL mrenv(kw,mass,mt,mc,lum,r,rc,aj,tm,lums(2),lums(3),&
                    lums(4),rzams,rtms,rg,menv,renv,k2)
    case(HeWD:Massless_Rem)
            rc = r
            menv = 1.0d-10
            renv = 1.0d-10
            k2 = 0.21d0
    end select


    t% pars% cenv_frac = menv/mt
    t% pars% env_frac = (mt-t% pars% McHe)/mt

    if (irecord>0) then
!        if (kw<10 .and. debug) print*, "delta in hrdiag",t% pars% delta
        t% pars% delta = 0.d0
        if (kw<10 .and. debug) print*,"saving pars",mt,mc,aj,kw,r
    else
        t% pars = old_pars
        t% post_agb = post_agb
        !tm and tn get calculated in star.f90
    endif
    if (irecord>0 .and. debug) print*,"finished hrdiag",mt,mc,aj,kw,id
!if (kw>1 .and. mc<=0.0) stop
!    print*,"finished hrdiag",t% pars% mass, t% pars% core_mass,t% pars% age
!    if ((id == 1) .and. (kw>1 .and. mc<=0.0)) stop

    
    nullify(t)
    end subroutine METISSE_hrdiag
