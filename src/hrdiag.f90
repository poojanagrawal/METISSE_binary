 subroutine hrdiag(mass,aj,mt,tm,tn,tscls,lums,GB,zpars,&
                                r,lum,kw,mc,rc,menv,renv,k2,mcx,id,irecord)
    use track_support
    use interp_support, only: interpolate_age
    use sse_support
    use remnant_support

    implicit none
    integer, intent(in), optional :: id
    real(dp) :: mass,aj,mt,tm,tn,tscls(20),lums(10),GB(10),zpars(20)
    real(dp) :: r,lum,mc,rc,menv,renv,k2,mcx
    integer :: kw,i, irecord, idd
    real(dp) :: rg,rzams,rtms
    type(star_parameters) :: old_pars

    logical :: debug,lost_envelope, post_agb
    type(track), pointer :: t

    debug = .false.
!    if (kw>1) print*,"started hrdiag",mt,mc,aj,kw
    t => NULL()
    if(present(id))then
        t => tarr(id)
        idd = id
    else
        t => tarr(1)
        idd = 1
    endif
    
    
    end_of_file = .false. !this is just the end of eep track, different to t_end defined in evolv1.f90

    !save input state
    lost_envelope = t% lost_envelope
    post_agb = t% post_agb
    old_pars = t% pars
  
    t% pars% mass = mt
    t% pars% phase = kw
    t% pars% age = aj
    t% irecord = irecord

!if(irecord>0) print*,"In Hrdiag aj,tn ",t% pars% age,mt,t% tr(i_age,t% ntrack),t% tr(i_age2,t% ntrack)


    select case(t% pars% phase)
        case(low_mass_MS:TPAGB)
            if (t% post_agb) then
                call evolve_after_agb(t)
            elseif (check_ge(t% pars% age,t% tr(i_age,t% ntrack))) then
                if (kw<5) print*,'WARNING: possible early end of file due to incomplete track beyond phase',kw
                end_of_file = .true.
!                    if (debug)
                print*,"end of file:aj,tn ",t% pars% age,t% tr(i_age,t% ntrack),t% tr(i_age2,t% ntrack)
                call check_remnant_phase(t)
            else
                t% pars% core_radius = -1.0
                call interpolate_age(t,t% pars% age)
                if (irecord>0 .and. debug)print*, "mt difference",t% pars% mass, mt, mt - t% pars% mass,kw
                t% pars% mass = mt
!                if (kw>=5 .and. irecord>0) print*, "mass",mt,mc,t% pars% core_mass,t% pars% McCO,kw
                !check if phase has changed
                do i = t% pars% phase,6
!                if (irecord>0) print*,"hr", t% times(i),t% times(i+1),aj,i
                    if (i== 0 .or. (.not. defined(t% times(i+1)))) exit
                    if (check_ge(t% pars% age,t% times(i))) then
                        t% pars% phase = i+1
                        if (debug) print*,"phase change",t% pars% age,t%times(i),i+1
                    endif
                end do

                !check envelope loss
                call check_env_loss(t)
                
                if (t% lost_envelope .and. t% pars% phase>6) then
                if (debug) print*, "envelope lost at", t% pars% mass, t% pars% age,t% pars% phase
                    call calculate_SSE_He_star(t,tscls,lums,GB,tm,tn)
                endif
            endif
            
            
        case(He_MS)
            call evolve_after_envelope_loss(t)
            rzams = t% He_pars% Rzams
        case(He_HG:He_GB)
            call evolve_after_envelope_loss(t)
            rzams = t% He_pars% Rzams
            call check_remnant_phase(t)
        case(HeWD:ONeWD)
            ! add post agb phase for he wd?
            call evolve_white_dwarf(t% pars)
        case(NS)
            call evolve_neutron_star(t% pars)
        case(BH)
            call evolve_black_hole(t% pars)
        case(Massless_Rem)
            ! return with whatever is the last value
    end select
      
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
        !TODO: should t% pars% core_radius be updated as well?
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
        t% lost_envelope = lost_envelope
        t% post_agb = post_agb
        !tm and tn get calculated in star.f90
    endif
    if (irecord>0 .and. debug) print*,"finished hrdiag",mt,mc,aj,kw,id
!    print*,"finished hrdiag",t% pars% mass, t% pars% core_mass,t% pars% age

    
    nullify(t)
    end subroutine hrdiag
