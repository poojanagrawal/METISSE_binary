 module remnant_support
    use track_support
    use sse_support
!    use z_support, only: Mup_core, Mec_core
    implicit none
    
    logical :: end_of_file, debug_rem

    !flags to be used while making decision for which method to use
    !for calculating properties of remnnants: neutron stars and black holes
    integer, parameter :: original_SSE = 0
    integer, parameter :: Belczynski2002 = 1
    integer, parameter :: Belczynski2008 = 2
    integer, parameter :: Eldridge_Tout2004 = 3
    integer, parameter :: Fryer2012 = 4

    !for white dwarfs
    integer, parameter :: Mestel = 0
    integer, parameter :: Modified_mestel = 1

    real(dp), parameter :: A_He = 4.d0
    real(dp), parameter :: A_CO = 16.d0
    real(dp), parameter :: tmin= 0.1d0

    !flags from SSE
    integer :: ns_flag = 0
    integer :: wd_flag = 0
    integer :: ec_flag = 0
    integer :: if_flag = 0

    real(dp) :: mc1 = 5.d0, mc2 = 7.6d0 !mass cutoffs for Belczynski methods

    contains

    subroutine check_remnant_phase(t)
    use track_support
    implicit none
    type(track),pointer, intent(inout) :: t
    real(dp) :: mc_max,mc_threshold,Mc11,Mcbagb
    integer :: j_bagb

    debug_rem = .false.
    mc_max= 0.d0
        !mcmax1 = 0.d0
    j_bagb=min(t% ntrack, TA_cHeB_EEP)
    Mcbagb = t% tr(i_he_core, j_bagb)

        if (t% pars% phase <=6) then       !without envelope loss
            Mc11 = 0.773* Mcbagb-0.35
            mc_max = MAX(M_ch,Mc11)
            !mc = MAX(mc_max,mc_threshold)
            !mc_max = MIN(pars% mass,mc_max)
            mc_threshold = t% pars% McCO
        else
            mc_max = max_core_mass_he(t% pars% mass, t% zams_mass)
            mc_threshold = t% pars% core_mass
            Mcbagb= t% zams_mass
        endif

        if (mc_max<=0.0) then
            print*,"fatal error: incorrect mc_max ", mc_max
            STOP
        endif

        if(mc_threshold>mc_max .or. abs(mc_max- mc_threshold)<tiny .or. end_of_file)then
            !mc = MIN(mc_max,mc_threshold)
            t% pars% core_mass = mc_threshold
            t% pars% age_old = t% pars% age

            if (debug_rem) then
                print*, "In remnant section"
                print*, "mass, core_mass, McCO, mc_max, 0.773*mcbagb-0.35"
                print*, t% pars% mass, t% pars% core_mass, t% pars% McCO, mc_max, Mc11
            end if

            !Mup_core, Mec_core are calculated in set_zparmeters routine of zfuncs
            if(t% pars% core_mass < M_ch)then
                    if(Mcbagb< Mup_core)then
                        t% pars% phase = CO_WD        !Zero-age Carbon/Oxygen White Dwarf
                        call post_agb_parameters(t)
                    else
                            if (ec_flag>0 .and. t% pars% McCO >= 1.372)then   !1.372<= mc< 1.44
                                t% pars% phase = NS
                                call initialize_ECSNe(t% pars)
                                if (debug_rem) print*,"ECSNe I: Mc< Mch, Mup> Mcbagb"
                            else
                                t% pars% phase = ONeWD    !Zero-age Oxygen/Neon White Dwarf
                                call post_agb_parameters(t)
                            endif
                    endif
            else !(mc>mch)
                !supernova 
                if(Mcbagb < Mup_core)then
                    ! Star is not massive enough to ignite C burning.
                    ! so no remnant is left after the SN
                    t% pars% phase = Massless_REM
                    call initialize_massless_rem(t% pars)

                else if(Mcbagb>= Mup_core .and. Mcbagb<= Mec_core)then
                    !Check for an electron-capture collapse of an ONe core.
                    !pars% McCO >= 1.372  !this used to be there
                    if(ec_flag>0) then
                        t% pars% phase = NS
                        call initialize_ECSNe(t% pars)
                        if (debug_rem) print*,"ECSNe II: Mc> Mch, Mup< Mbagb< Mec"
                    else
                        call check_ns_bh(t% pars)
                    endif
                else
                    call check_ns_bh(t% pars)
                endif
            endif
        if (debug_rem .and. t% pars% phase>9) print*,"In remnant phase", phase_label(t% pars% phase+1)," , mass", t% pars% mass
        endif
    end subroutine check_remnant_phase

    subroutine post_agb_parameters(t)
        type(track), pointer, intent(inout) :: t
        
        if(construct_wd_track) then
            t% agb% phase_wd = t% pars% phase
            t% pars% phase = TPAGB
            t% post_agb = .true.
            t% pars% age_old = 0.0

            t% agb% age = t% pars% age
            t% agb% lum = t% pars% luminosity
            t% agb% radius = t% pars% radius
            call evolve_after_agb(t)
            if (debug_rem) print*, "In post-agb phase, mass = ", t% pars% mass
!            print*,t% pars% luminosity, t% pars% radius
        else
            t% zams_mass = t% pars% mass
            call initialize_white_dwarf(t% pars)
        endif
    end subroutine

    subroutine evolve_after_agb(t)
        type(track),pointer, intent(inout) :: t

        real(dp) :: alfa, beta,dt,r3
        real(dp) :: radius_wd,lum_wd,mass_wd
        real(dp) :: t1,t2,t_post_agb
        !this phase doesn't account for accretion, HeWD -TODO check

        t1 = fit_Z_t1(initial_Z)*t% MS_time
        t2 = 10*t1
        t_post_agb = t1+t2
        t% times(TPAGB) = t% nuc_time + t_post_agb
        mass_wd = t% pars% core_mass
        radius_wd = calculate_wd_radius(mass_wd)
        lum_wd = calculate_wd_lum(mass_wd, 0.d0, A_CO)  !xx = a_co
            dt= t% pars% age- t% agb% age
!            print*,"aj",t% pars% age,age_agb, t% pars% age- age_agb
!            print*, t% nuc_time,t% times(TPAGB),dt
            alfa = 0d0; beta = 0d0
            r3 = 0.3*(t% agb% radius+radius_wd)

            if (dt<t1) then
                alfa = dt/t1
                beta = 1d0-alfa
                t% pars% mass = alfa* mass_wd + beta*t% pars% mass
                t% pars% radius = alfa* r3 + beta* t% agb% radius
                t% pars% luminosity = alfa* 0.9*t% agb% lum + beta*t% agb% lum
                t% pars% extra =1
            else
                alfa = (dt-t1)/t2
                beta = 1d0-alfa
                t% pars% radius = (radius_wd*r3)/(alfa*(r3-radius_wd)+radius_wd)
                t% pars% luminosity = alfa* lum_wd + beta* 0.9*t% agb% lum
                t% pars% extra = 2
            endif
            if (check_ge(t% pars% age,t% times(TPAGB))) then
                t% post_agb = .false.
                t% pars% extra = 0
                t% pars% age_old = t% pars% age
                t% pars% phase = t% agb% phase_wd   !TODO: Should phase_wd be recalculated?
                t% zams_mass = t% pars% mass
                call initialize_white_dwarf(t% pars)
            endif
    end subroutine

    subroutine initialize_white_dwarf(pars)
    type(star_parameters),intent(inout) :: pars

        if(if_flag>=1) call check_IFMR(pars% mass, pars% core_mass)
        pars% mass = pars% core_mass
        !pars% McHe = 0.0
        !pars% McCO = 0.0
        pars% age = 0.0

        call evolve_white_dwarf(pars)
        if (debug_rem) print*, "Remnant phase = ", phase_label(pars% phase+1), ", mass =", pars% mass
    end subroutine

    subroutine check_IFMR(mass, mc)
        real(dp), intent(in) :: mass
        real(dp) :: mc

        ! Invoke WD IFMR from HPE, 1995, MNRAS, 272, 800.
          if(Z04>= 1E-08)then
             mc = MIN(0.36+0.104*mass,0.58+0.061*mass)
             mc = MAX(0.54+0.042*mass,mc)
             if(mass<1.0) mc= 0.46
          else
             mc= MIN(0.29+0.178*mass,0.65+0.062*mass)
             mc= MAX(0.54+0.073*mass,mc)
          endif
          mc= MIN(M_ch,mc)
    end subroutine

    subroutine check_ns_bh(pars)
    type(star_parameters) :: pars
    real(dp):: Mrem
    logical :: debug

    debug = .false.
            if (debug) print*,"In CCSNe section"
            pars% age = 0.0
            Mrem = calculate_NSBH_mass(pars% core_mass,pars% mass)
            if(debug) print*, "Mrem= ", Mrem

            if(Mrem <= Max_NS_mass)then
                pars% phase = NS       !Zero-age Neutron star
                pars% mass = calculate_gravitational_mass(Mrem, pars% phase)
                call evolve_neutron_star(pars)
            else
                pars% phase = BH       !Zero-age Black hole
                pars% mass = calculate_gravitational_mass(Mrem, pars% phase)
                call evolve_black_hole(pars)

            endif
            if (debug) print*, "NS/BH mass from", BHNS_mass_scheme,"scheme = ",pars% mass
    end subroutine

    real(dp) function calculate_NSBH_mass(Mc,Mt) result(Mrem)
        real(dp), intent(in):: Mc,Mt
        real(dp) :: Mc_FeNi

        Mrem = Mc
        select case(ns_flag)
            case(Belczynski2002)
                !Use FeNi core mass given by Belczynski et al. 2002, ApJ, 572, 407.
                if(Mc< 2.5d0)then
                    Mc_FeNi = 0.161767d0*Mc+ 1.067055d0
                else
                    Mc_FeNi = 0.314154d0*Mc+ 0.686088d0
                endif
                Mrem = calculate_remnant_mass(Mc,Mc_FeNi,Mt)

            case(Belczynski2008)
                !Use FeNi core mass given by Belczynski et al. 2008, ApJSS, 174, 223.
                if(Mc< 4.82d0)then
                    Mc_FeNi = 1.5d0
                elseif(Mc< 6.31d0)then
                    Mc_FeNi = 2.11d0
                elseif(Mc< 6.75d0)then
                    Mc_FeNi = 0.69d0*Mc- 2.26d0
                else
                    Mc_FeNi = 0.37d0*Mc- 0.0828d0
                endif
                Mrem = calculate_remnant_mass(Mc,Mc_FeNi,Mt)

            case(Eldridge_Tout2004)
                !Use remnant masses based on Eldridge & Tout 2004, MNRAS, 353, 87.
                !NOT checked
                if(Mc<6.d0)then
                    Mc_FeNi = 1.44d0
                else
                    Mc_FeNi = (1.4512017d0*Mc)-(6.5913737d-03*Mc*Mc)-6.1073371d0
                endif
                Mrem = Mc_FeNi

            case(original_SSE)
                !Use the original SSE NS/BH mass.
                Mrem = 1.17d0 + 0.09d0*Mc
        end select
    end function

    real(dp) function calculate_remnant_mass(mc,mcfeni,Mt) result (Mrem)
        real(dp), intent(in) :: mc,mcfeni,mt

        !For Belczynski methods calculate the remnant mass from the FeNi core.
        if(mc<= mc1)then
            Mrem = mcfeni
        elseif(mc< mc2)then
            Mrem = mcfeni + (mc- mc1)*(mt - mcfeni)/(mc2-mc1)
        else
            Mrem = Mt
        endif
    end function calculate_remnant_mass

    real(dp) function calculate_gravitational_mass(mass,phase)
        real(dp) :: mass, mt
        integer,intent(in) :: phase
        mt = 0.0
        if(ns_flag>=2)then
            if (phase == NS) mt = quadratic(0.075d0, 1.d0, -mass)
            if (phase == BH) mt = 0.9d0* mass
        endif
        calculate_gravitational_mass = mt
    end function
    
    subroutine evolve_white_dwarf(pars)
    type(star_parameters):: pars
    real(dp) :: xx
    logical :: debug

    debug = .false.
    
    if (debug) print*,"In evolve_white_dwarf"
        pars% core_mass = pars% mass             !to account for accretion
        if(pars% phase == HeWD)then
            xx = a_he
        else
            xx = a_co
        endif

        if(pars% core_mass >= M_ch)then    !Mch or 1.372 !check for ec_flag
            if(pars% phase == ONeWD)then
                if (debug) print*, "ECSNe:NS from accn onto WD", pars% mass
                pars% phase = NS
                call initialize_ECSNe(pars)
            else
                !Accretion induced supernova with no remnant
                if (debug) print*, "ECSNe: massless remnant from accn onto WD", pars% mass
                pars% phase = Massless_REM
                call initialize_massless_rem(pars)
            endif
        else
            pars% luminosity = calculate_wd_lum(pars% mass, pars% age, xx)
            pars% radius = calculate_wd_radius(pars% mass)
            if(pars% mass < 0.0005) pars% radius= MIN(pars% radius,0.01d0)
            if (debug) print*, "Evolving WD", pars% mass, pars% luminosity, pars% radius
        endif
    end subroutine

    real(dp) function calculate_wd_lum(mass, age, xx) result(lum)
    real(dp), intent(in) :: mass, age, xx
    real(dp) :: fac
        lum = 0.d0
        select case(wd_flag)
        case(Mestel)            ! Mestel cooling
            Lum = 635.d0* mass* Z04/(xx*(age+tmin))**1.4
        case(Modified_mestel)       ! modified-Mestel cooling
            if(age < 9000.0)then
                Lum = 300* mass* Z04/(xx*(age+tmin))**1.18
            else
                fac = (9000.1*xx)**5.3
                Lum = 300*fac* mass* Z04/(xx*(age+tmin))**6.48
            endif
        end select

        !tmin= ((635.d0*pars% mass*Z04/lum)**(1.0/1.4))/xx
        !lum = (635.d0*mass*(Z04))/((xx*tmin)**1.4)
        !if (t% post_agb) lum = min(lum,t% agb% lum)
    end function

    real(dp) function calculate_wd_radius(mass) result(radius)
    real(dp), intent(in) :: mass
        radius = sqrt((M_ch/mass)**pow-(mass/M_ch)**pow)
        radius = max(1.4d-5,0.0115*radius)
        !radius  = 0.0115*SQRT(MAX(1.48204d-06,(M_ch/mass)**pow-(mass/M_ch)**pow))
        radius = MIN(0.1d0,radius )
    end function

    subroutine evolve_neutron_star(pars)
    type(star_parameters) :: pars
        !Neutron Star
        pars% core_mass= pars% mass
        if(pars% core_mass > Max_NS_mass)then
            pars% phase = BH  !Accretion induced Black Hole
            !pars% age_old = pars% age
            pars% age = 0.0
            call evolve_black_hole(pars)
        else
            pars% luminosity = 0.02*(pars% mass**0.67)/(MAX(pars% age,0.1d0))**2
            pars% radius= 1.4E-05
            !pars% McCO = 0.0
            !pars% McHe = 0.0
        endif

        !print*,"I am in evolve_neutron_star"
    end subroutine

    subroutine evolve_black_hole(pars)
    type(star_parameters) :: pars
        !Black hole
        !pars% mass has been calculated during remnant check
        pars% core_mass = pars% mass 
        pars% luminosity = 1.0E-10
        pars% radius = 4.24E-06*pars% mass
        !pars% McCO = 0.0
        !pars% McHe = 0.0

    end subroutine

    subroutine initialize_massless_rem(pars)
        type(star_parameters) :: pars
        pars% age = 0.0
        pars% mass = 0.0
        pars% luminosity = 1E-10
        pars% radius = 1E-10
        pars% Teff = 1E-10
        pars% core_mass = 0.0
        pars% McCO = 0.0
        pars% McHe = 0.0
    end subroutine

    subroutine initialize_ECSNe(pars)
    type(star_parameters) :: pars
        !pars% age_old = pars% age
        pars% age = 0.0
        pars% mass = 1.26
        pars% core_mass= pars% mass
        pars% luminosity = 0.02*(pars% mass**0.67)/(MAX(pars% age,0.1d0))**2
        pars% radius= 1.4E-05
    end subroutine

    subroutine check_env_loss(t)

    type(track), pointer, intent(inout) :: t
    real(dp) :: HeI_time, HeB_time
    logical :: debug

    debug = .false.
    t% lost_envelope = (t% initial_mass <10.0 .and. t% pars% core_mass>= t% pars% mass) &
                        .or. (t% initial_mass>=10.0 .and. abs(t% pars% core_mass-t% pars% mass)<0.01)
    
    !if (t% irecord>0)print*,"age, core mass, mass", t%  pars% age, t% pars% core_mass ,t% pars% mass
    if(t% lost_envelope)then
        if (debug) print*,"Lost envelope at phase", t% pars% phase
        if (debug) print*,"age, core mass, mass", t% pars% age, t% pars% core_mass, t% pars% mass

        t% pars% age_old = t% pars% age

        select case(t% pars% phase)
            case(HG:RGB)   !HG or RGB
                if(t% pars% mass< Mhef)then
                    t% pars% phase = HeWD      !Zero-age helium white dwarf
                    t% pars% core_mass = t% pars% mass
                    call initialize_white_dwarf(t% pars)
                else
                    t% pars% phase = He_MS       !Zero-age helium star
                    t% pars% age = 0.d0
                    t% zams_mass = t% pars% mass
                    t% pars% core_mass = 0.d0
                    t% pars% McCO = 0.d0
                    t% pars% McHe = 0.d0

                    call calculate_he_timescales(t)
                    t% lost_envelope = .true.
                    call evolve_after_envelope_loss(t)
                endif

            case(HeBurn)  !core he Burning
                t% pars% phase = He_MS
                t% zams_mass = t% pars% mass
                t% pars% core_mass = 0.d0
                t% pars% McCO = 0.d0
                t% pars% McHe = 0.d0
                call calculate_he_timescales(t)
                HeI_time = t% times(3)
                HeB_time = t% times(4)-t% times(3)
                t% pars% age = t% MS_time*((t% pars% age- HeI_time)/HeB_time)

                t% lost_envelope = .true.
                call evolve_after_envelope_loss(t)
                
            case(EAGB) !eAGB
                t% pars% phase = He_GB       !Evolved naked He star
                t% pars% mass = t% pars% core_mass
                t% zams_mass = t% pars% mass
                t% pars% McHe = t% pars% mass
                t% pars% core_mass = t% pars% McCO

                call calculate_he_timescales(t)
                t% pars% age = He_GB_age(t% pars% core_mass,t% times(8), &
                                t% times(9),t% He_pars% D, t% He_pars% Mx)
                t% pars% age = MAX(t% pars% age,t% MS_time)
                t% lost_envelope = .true.
                call evolve_after_envelope_loss(t)
            case(TPAGB)
                call check_remnant_phase(t) !TODO -  check if this works
        end select
    endif
    end subroutine

    subroutine evolve_after_envelope_loss(t)
    type(track),pointer, intent(inout) :: t

    real(dp) :: rg,tau  !,McHeI
    logical :: debug

    debug = .false.
    if (debug) print*,"In evolve_after_envelope_loss: phase age",t% pars% phase, t% pars% age,t% ms_time,t% zams_mass
            
    t% He_pars% Rzams = radius_He_ZAMS(t% pars% mass)
    !if (check_le(t% pars% age ,t% MS_time)) then
    if(t% pars% age <t% MS_time .and. abs(t% pars% age-t% MS_time)>tiny)then
        !Helium Main Sequence
        !From SSE: "Star has no core mass and hence no memory of its past
        !which is why we subject mass and mt to mass loss for this phase."

        t% pars% phase = He_MS
        tau = t% pars% age/t% MS_time
        t% He_pars% Lzams = lum_He_ZAMS(t% zams_mass)
        t% pars% luminosity = lum_He_MS(t% zams_mass, t% He_pars% Lzams,tau)
        t% pars% radius = radius_He_MS(t% pars% mass,t% He_pars% Rzams,tau)
        rg = t% He_pars% Rzams
        t% pars% core_mass = 0.0
        t% pars% age_old = t% pars% age
        !McHeI= He_McDu = core_mass_He_GB(lums(4), D,lx)   !core mass at He Ignition
        !print*, t% pars% mass, t% McHeI

        !TODO -- this step needs to checked specially for binaries
        ! TODO -- also evolve this HeWD
!                 j_HeI = min(t% ntrack, cHeIgnition_EEP)
!                McHeI = t% tr(i_he_core, j_HeI)
!                if(t% pars% mass < McHeI ) t% pars% phase = HeWD
     else
        !Helium Shell Burning
        t% pars% phase = He_HG
        t% pars% luminosity = lum_He_GB(t% pars% age,t% times(8),t% times(9),&
                                t% times(10),t% He_pars% D)
        t% pars% radius = radius_He_HG(t% pars% mass,t% pars% luminosity,&
                                t% He_pars% Rzams,t% He_pars% LtMS)
        rg = radius_He_GB(t% pars% luminosity)  !radius on he giant branch
        if(t% pars% radius >= rg)then
           t% pars% phase = He_GB
           t% pars% radius = rg
        endif
        t% pars% core_mass = core_mass_He_GB(t% pars% luminosity,t% He_pars% Lx,t% He_pars% D)
        t% pars% McHe = t% pars% mass
        t% pars% McCO = t% pars% core_mass
    endif
    if (debug) print*,"End: Phase",t% pars% phase ," core mass",t% pars% core_mass
!        print*,"lum", t% pars% luminosity, "rad", t% pars% radius

    end subroutine

    subroutine set_remnant_scheme()
        if (WD_mass_scheme == 'Mestel') then
            wd_flag = Mestel
        elseif (WD_mass_scheme == 'Modified_mestel') then
            wd_flag = Modified_mestel
        else
            print*,"Error: Invalid Option for WD_mass_scheme."
            print*,"Choose from 'Mestel' and 'Modified_mestel'. "
            STOP
        endif

        if (BHNS_mass_scheme == 'original_SSE') then
            ns_flag = original_SSE
        elseif (BHNS_mass_scheme == 'Belczynski2002') then
            ns_flag = Belczynski2002
        elseif (BHNS_mass_scheme == 'Belczynski2008') then
            ns_flag = Belczynski2008
        elseif (BHNS_mass_scheme == 'Eldridge_Tout2004') then
            ns_flag = Eldridge_Tout2004
        else
            print*,"Error: Invalid Option for BHNS_mass_scheme."
            print*,"Choose from 'original_SSE','Belczynski2002','Belczynski2008','Eldridge_Tout2004'. "
            STOP
        endif
        !call cutoffs_for_Belzynski_methods(ns_flag,mc1,mc2)

        !if (verbose) print*,"For Belzynski methods", mc1,mc2
        if (allow_electron_capture) ec_flag = 1
        if (Use_Initial_final_mass_relation) if_flag = 1
    end subroutine set_remnant_scheme

    !for constructing post-main sequence phase if construct_wd_track is true
    real(dp) function fit_Z_t1(Z)
    real(dp):: Z,x
        x = log10(Z)
        !y = 10**(-0.1*x*x-0.46*x-5.1)
        fit_Z_t1 = 10**(-0.136*x*x-1.117*x-6.256)
        return
    end function


    subroutine calculate_rc(t, tscls,zpars,rc)
     ! Calculate the core radius 
     implicit none
     type(track), pointer :: t
     real(dp) :: tscls(20), zpars(20)
     real(dp) :: tau, lx,rx, rc,am,mt,mc,aj
     real(dp) :: tbagb,mass,lums1,lums2,tn
     
        mass = t% zams_mass
        mt = t% pars% mass
        mc = t% pars% core_mass
        aj = t% pars% age
        tn = t% nuc_time
        tau = 0.d0
       select case(t% pars% phase)
         case(low_mass_MS:MS)
             rc = 0.d0
         case(HG:RGB)
             if(mass.gt.zpars(2))then
                 rc = radius_He_ZAMS(mc)
             else
                 rx = calculate_wd_radius(mc)
                 rc = 5.d0*rx
             endif
          case(HeBurn)
              tau = (aj - tscls(2))/tscls(3)     !tau = (aj - t_HeI)/t_He
              rx = radius_He_ZAMS(mc)
              !Following is akin to rzhef of SSE
              !TODO: following two lines may cause dt error
              ! might need interpolation
              am = MAX(0.d0,0.4d0-0.22d0*LOG10(mc))
              rc = rx*(1.d0+am*(tau-tau**6))
         case(EAGB)
    !            kwp = 9
             mc = t% pars% McHe
             tbagb = t% times(HeBurn)
             if(tn.gt.tbagb) tau = 3.d0*(aj-tbagb)/(tn-tbagb)
             lx = lmcgbf(t% pars% McCO,t% He_pars% D, t% He_pars% Mx)
             lums1 = lum_He_ZAMS(mc)
             lums2 = lum_He_MS(mc,lums1,1.d0)
             if(tau.lt.1.d0) lx = lums2*(lx/lums2)**tau
             rx = radius_He_HG(mc,lx,radius_He_ZAMS(mc),lums2)
             rc = MIN(rx,radius_He_GB(lx))
         case(TPAGB: He_GB)
             if (t% pars% phase == He_MS) then
                 rc = 0.d0
             else
                 rx = calculate_wd_radius(mc)
                 rc = 5.d0*rx
             endif
         end select
        rc = MIN(rc,t% pars% radius)
    end subroutine

!    subroutine cutoffs_for_Belzynski_methods(ns_flag,mc1,mc2)
!        use track_support
!        use z_support, only : Mcrit
!        implicit none
!
!        integer, intent(in) :: ns_flag
!        real(dp), intent(out)  :: mc1, mc2
!        type(track) :: x
!        integer :: j_ntrack
!
!        if (ns_flag == Belczynski2002 .or. ns_flag == Belczynski2008) then
!            if(Mcrit(9)% mass>=20d0) then
!                call star(20d0,x)
!                j_ntrack = min(final_eep,x% ntrack)
!                mc1 = x% tr(i_co_core,j_ntrack)
!                call dealloc_track (x)
!            endif
!
!            if(Mcrit(9)% mass>=42d0) then
!                call star(42d0,x)
!                j_ntrack = min(final_eep,x% ntrack)
!                mc2 = x% tr(i_co_core,j_ntrack)
!                call dealloc_track (x)
!            endif
!        endif
!    end subroutine
end module remnant_support



