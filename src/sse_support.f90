module sse_support
!module to help interface metisse with sse
    use track_support
    use z_support, only: Mcrit
    implicit none
    
 !SSE parameters required if envelope is lost
    real(dp) :: AHe = 8.0d-05   !7.66d-5 !Msun Lsun**-1 MYR**-1
    real(dp) :: B = 4.1d+04
    real(dp) :: p = 5.0
    real(dp) :: q = 3.0
    real(dp), parameter :: pow = 2.d0/3.d0
    real(dp), parameter :: M_ch = 1.44d0

!
!    some of these needs recalculation at every step and others need to be a part of t
!    initialize them when star loses envelope or becomes WD
!    real(dp) :: Tinf1,Tinf2,Tx

contains

    subroutine calculate_timescales(t)
        !subroutine to assign sse phases to the interpolated track
        !also calculate timescales associated with different phases
        implicit none
        
        type(track), pointer,intent(inout) :: t
        integer :: i,j_bgb
        real(dp) :: mass
        real(dp), pointer :: age(:)=> NULL()
        logical :: debug
        debug = .false.

        age => t% tr(i_age2,:)
        mass = t% initial_mass

        t% MS_time = age(TAMS_EEP) - age(ZAMS_EEP)
        do i = 1, t% neep
            if (t% eep(i) == TAMS_EEP) then    !MS
                t% times(MS) = age(TAMS_EEP)

            elseif (t% eep(i) == cHeIgnition_EEP) then
                !Herztsprung gap
                t% times(HG) = age(cHeIgnition_EEP)
                !t% times(HG) gets modified to t(BGB)
                ! if RGB phase is present
                t% times(RGB) = t% times(HG)

            elseif (t% eep(i) == TA_cHeB_EEP) then
                !red_HB_clump /core He Burning
                t% times(HeBurn) = age(TA_cHeB_EEP)

            elseif (t% eep(i) == cCBurn_EEP) then
                !EAGB/ core C burning
                t% times(EAGB) = age(cCBurn_EEP)
                
            elseif (t% eep(i) == TPAGB_EEP) then
                !AGB
                t% times(EAGB) = age(TPAGB_EEP)
                
            elseif (t% eep(i) == post_AGB_EEP) then
                !TP-AGB :only for low_inter mass stars
                t% times(TPAGB) = age(post_AGB_EEP)
            endif
        enddo
        !print*,"bgb",BGB_EEP,identified(BGB_EEP)

        t% times(11) = age(min(Final_EEP,t% ntrack))
        !Todo: nuc_time should be for WR phase
        t% nuc_time = t% times(11)
        !determine the base of the giant branch times, if present
        if (mass > Mcrit(2)% mass .and. mass< Mcrit(5)% mass) then
            if (identified(BGB_EEP)) then
                !Red giant Branch
                t% times(HG) = age(BGB_EEP)
            elseif (t% ntrack >TAMS_EEP) then
                j_bgb =  base_GB(t)
                if (j_bgb>0) then
                    j_bgb = j_bgb+TAMS_EEP-1
                    !Red giant Branch
                    t% times(HG) = age(j_bgb)
                elseif (debug) then
                    print*, "Unable to locate BGB ", j_bgb
!                    STOP
                end if
            endif
        endif

        nullify(age)
    end subroutine calculate_timescales


    subroutine calculate_SSE_parameters(t,zpars,tscls,lums,GB,tm,tn)
        type(track), pointer, intent(in) :: t
        real(dp) :: tscls(20),lums(10),GB(10),zpars(20),tm,tn
    
        ! Other SSE parameters
        call calculate_SSE_tscls(t, tscls,tm,tn)
        call calculate_SSE_GB(t, zpars, GB)
        call calculate_SSE_lums(t, zpars, lums)
        lums(6) = GB(4)*GB(7)**GB(5)
        
    end subroutine calculate_SSE_parameters
    
    !TSCLS (all timescales are in Myr units)
    !1; BGB      Base of giant  branch      2; He ignition   3; He burning
    !4; Giant t(inf1)    5; Giant t(inf2) 6; Giant t(Mx)
    !7; FAGB t(inf1)     8; FAGB t(inf2)  9; FAGB  t(Mx)   First Giant Branch
    !10; SAGB t(inf1)    11; SAGB t(inf2) 12; SAGB  t(Mx) Second Giant Branch
    !13; TP            Thermal Pulsations  14; t(Mcmax)  maximum age
    
    subroutine calculate_SSE_tscls(t, tscls,tm,tn)
    !only for phase upto 6 (TPAGB)
    implicit none
    type(track), pointer, intent(in) :: t
    real(dp) , intent(out):: tscls(20),tm,tn

    tscls = -1.d0

    tm = t% times(MS)
    tscls(1) = t% times(HG)
    tscls(2) = t% times(RGB)        !for stars that don't have RGB phase, this corresponds to t% times(HG)
    tscls(3) = t% times(HeBurn)-t% times(RGB)

   
    tscls(13) = t% times(EAGB)
    
    tscls(14) = t% times(11)
    tn = tscls(14)
    if (tscls(13)<0) tscls(13) = tn

    !filling in the gaps
    tscls(4) = t% times(RGB)
    tscls(5) = tscls(4)
    tscls(6) = tscls(4)
    
    tscls(7) = t% times(5)
    tscls(8) = tscls(7)
    tscls(9) = tscls(7)
    
    tscls(10) = min(tn,tscls(13))
    tscls(11) = tscls(10)
    tscls(12) = tscls(10)

    end subroutine

   ! LUMS (all luminosities are in solar units)
    !1; ZAMS             2; End MS        3; BGB
    !4; He ignition      5; He burning    6; L(Mx)
    !7; BAGB             8; TP
    subroutine calculate_SSE_lums(t, zpars,lums)
    !only for phase upto 6 (TPAGB)
    implicit none
     type(track), pointer, intent(in) :: t
    real(dp),  intent(in) :: zpars(20)
    
    integer :: i,j_bgb
    real(dp) :: lums(10), mass
    real(dp), allocatable :: Lumt(:)
    
    logical :: debug
    debug = .false.
    Lumt = 10**(t% tr(i_logL,:))
    mass = t% initial_mass
            
    do i = 1, t% neep
        if (t% eep(i) == ZAMS_EEP) lums(1) = Lumt(ZAMS_EEP)
        if (t% eep(i) == TAMS_EEP) lums(2) = Lumt(TAMS_EEP)
        if (t% eep(i) == cHeIgnition_EEP) lums(4) = Lumt(cHeIgnition_EEP)
        if (t% eep(i) == cHeBurn_EEP) lums(5) = Lumt(cHeBurn_EEP)
        if (t% eep(i) == TA_cHeB_EEP) lums(7)= Lumt(TA_cHeB_EEP)
        if (t% eep(i) == TPAGB_EEP) lums(8) = Lumt(TPAGB_EEP)
    enddo
    lums(6) = 0.d0
    
    !TODO check this
     if (mass .lt. zpars(3) .and. identified(BGB_EEP)) then
            j_bgb = min(BGB_EEP, t% ntrack)
            lums(3) = Lumt(j_bgb)
    else
            lums(3) = lums(4)
    endif
    
    deallocate(Lumt)
    end subroutine

    ! GB = giant branch parameters
    !1; effective A(H)   2; A(H,He)       3; B
    !4; D                5; p             6; q
    !7; Mx               8; A(He)         9; Mc,BGB
    subroutine calculate_SSE_GB(t,zpars, GB)
     !only for phase upto 6 (TPAGB)
    implicit none

    type(track), pointer, intent(in) :: t
    real(dp),  intent(in) :: zpars(20)
    real(dp) :: GB(10)
    real(dp) :: mass,dx
    integer :: j
    
    mass = t% pars% mass
    
    GB(1) = MAX(-4.8d0,MIN(-5.7d0+0.8d0*mass,-4.1d0+0.14d0*mass))
    GB(1) = 10.d0**GB(1)
    GB(2) = 1.27d-05
    GB(8) = 8.0d-05
    GB(3) = MAX(3.0d+04,500.d0 + 1.75d+04*mass**0.6d0)
    if(mass.le.2.0)then
         GB(4) = zpars(6)
         GB(5) = 6.d0
         GB(6) = 3.d0
    elseif(mass.lt.2.5)then
         dx = zpars(6) - (0.975d0*zpars(6) - 0.18d0*2.5d0)
         GB(4) = zpars(6) - dx*(mass - 2.d0)/(0.5d0)
         GB(5) = 6.d0 - (mass - 2.d0)/(0.5d0)
         GB(6) = 3.d0 - (mass - 2.d0)/(0.5d0)
    else
         GB(4) = MAX(-1.d0,0.5d0*zpars(6) - 0.06d0*mass)
         GB(4) = MAX(GB(4),0.975d0*zpars(6) - 0.18d0*mass)
         GB(5) = 5.d0
         GB(6) = 2.d0
    endif
    GB(4) = 10.d0**GB(4)
    GB(7) = (GB(3)/GB(4))**(1.d0/(GB(5)-GB(6)))
      
      
     if(mass.le.zpars(2).and. identified(BGB_EEP))then                   !MHeF
            j = min(t% ntrack, BGB_EEP)
            GB(9) = t% tr(i_he_core, j )
      elseif(mass.le.zpars(3))then         !MFGB
            j = min(t% ntrack, cHeIgnition_EEP)
            GB(9) = t% tr(i_he_core, j )
      else
            j = min(t% ntrack, cHeIgnition_EEP)
            GB(9) = t% tr(i_he_core, j )
      endif
    end subroutine

!     subroutine calculate_he_timescales(t,t% He_pars% LtMS, Mx, Tinf1, Tx, Tinf2)
     subroutine calculate_He_timescales(t)
      !only for phases 7 to 9
            type(track), pointer, intent(inout) :: t
            real(dp) :: Tinf1, Tx, Tinf2
            real(dp) :: D,Mx,Lx,LtMS
            real(dp) :: mc1, Tmax, mc_max,tm
            logical :: debug

            debug = .false.
            
    !        Calculate Helium star Main Sequence lifetime.
            t% MS_time = time_He_MS(t% zams_mass)
            tm = t% MS_time
            
    !       Zero- and terminal age Helium star main sequence luminosity
            t% He_pars% Lzams = lum_He_ZAMS(t% zams_mass)
            LtMS = lum_He_MS(t% zams_mass, t% He_pars% Lzams,1.d0)
            t% He_pars% LtMS = LtMS
            
    !       Set the Helium star GB parameters
            
            D = 5.5d+04/(1.d0+0.4d0* t% zams_mass**4)
            t% He_pars% D = D
            Mx = (B/D)**(1.d0/(p-q))
            t% He_pars% Mx = Mx
            ! Change in slope of giant L-Mc relation.
            lx= D*Mx**p
            t% He_pars% lx = lx

    !       Set Helium star GB timescales
            !core_mass_THeMS; mc1
            if(LtMS <= lx)then
                mc1 = (LtMS/D)**(1.d0/p)
            else
                mc1 = (LtMS/B)**(1.d0/q)
            endif

            Tinf1 = tm + (1.d0/((p-1.d0)*AHe*D))*mc1**(1.d0-p)
            Tx = Tinf1 - (Tinf1 - tm)*((Mx/mc1)**(1.d0-p))
            Tinf2 = Tx + (1.d0/((q-1.d0)*AHe*B))*Mx**(1.d0-q)

    !        Get an idea of when Mc = MIN(Mt,Mc,C,max) on the GB
        
            mc_max = max_core_mass_he(t% pars% mass, t% zams_mass)
!TODO: Next line requires a check
!Note that this is done to avoid negative timesteps that result from more massive cores than what sse formulae predict
            mc_max = max(mc_max,t% pars% core_mass+1.d-10)
            if(debug)print*, 'core',mc_max,t% pars% core_mass,t% pars% mass,t% zams_mass

            if(mc_max.le.Mx)then
                Tmax = Tinf1 - (1.d0/((p-1.d0)*AHe*D))*(mc_max**(1.d0-p))
            else
                Tmax = Tinf2 - (1.d0/((q-1.d0)*AHe*B))*(mc_max**(1.d0-q))
            endif
            if (debug) print*, 'tmax', tmax, mc_max,mx,mc1
            Tmax = MAX(Tmax,t% MS_time)
            t% nuc_time = Tmax
            
            t% times(He_MS)= t% MS_time
            t% times(8) = Tinf1     !these 8-10 here do not represent the phase numbers
            t% times(9) = Tinf2
            t% times(10) = Tx
            if(debug) print*,"He timescales", tinf1, tx, tinf2, tmax,t% MS_time,t% nuc_time
        return
        end subroutine calculate_He_timescales

    subroutine calculate_SSE_He_star(t,tscls,lums,GB,tm,tn)
            type(track), pointer, intent(in) :: t
            real(dp) :: tscls(20),lums(10),GB(10),tm,tn

                GB(3) = 4.1d+04
                GB(4) = t% He_pars% D
                GB(5) = 5.d0
                GB(6) = 3.d0
                GB(7) = t% He_pars% Mx
                GB(8) = 8.0d-05
                
                lums(1) = t% He_pars% Lzams
                lums(2) = t% He_pars% LtMS
                lums(6) = t% He_pars% lx
                
                tscls(1) = t% MS_time
                tscls(4) = t% times(8)
                tscls(5) = t% times(9)
                tscls(6) = t% times(10)
                tscls(14) = t% nuc_time
                tm = t% MS_time
                tn = t% nuc_time
    end subroutine calculate_SSE_He_star
    
    
    real(dp) FUNCTION lum_He_ZAMS(m)
    !old name = lzhef(mass)
    implicit none
    real(dp):: m, top, bottom

        !A function to evaluate Naked Helium star 'ZAMS' luminosity

        top = 1.5262d+04*m**(41.d0/4.d0)
        bottom = 0.0469d0 + 31.18d0*m**6 + 29.54d0*m**7.5 + m**9
        lum_He_ZAMS= top/bottom
    return
    end

    real(dp) FUNCTION radius_He_ZAMS(m)
        ! old name = rzhef(mt)
        implicit none
        real(dp) :: m,top,bottom
        !Function to evaluate Helium star 'ZAMS' radius
        top = 0.2391*m**4.6
        bottom = 0.0065 + 0.162*m**3 + m**4
        radius_He_ZAMS = top/bottom
    return
    end

    real(dp) FUNCTION time_He_MS(m)
    !old name = themsf
    implicit none
    real(dp) :: m,top,bottom
        !A function to evaluate Helium star main sequence lifetime

        top = 0.4129d0 + 18.81d0*m**4 + 1.853d0*m**6
        bottom = m**6.5d0
        time_He_MS = top/bottom
        return
    end
    !themsf = (0.4129d0 + 18.81d0*m**4 + 1.853d0*m**6)/m**(13.d0/2.d0)


    real(dp) FUNCTION lum_He_MS(mass,lum,tau)
    implicit none
    real(dp), intent(in) :: mass,lum,tau
    real(dp) :: alpha
    !derived from lzhef and lums(2) calculations for He star from star.f

    !A function to evaluate Helium star luminosity on the MS

        alpha = MAX(0.d0,0.85d0-0.08d0*mass)
        lum_He_MS = lum *(1.d0+ 0.45d0*tau+ alpha*tau**2)
    return
    end


    real(dp) FUNCTION radius_He_MS(mass,radius,tau)
    implicit none
    real(dp), intent(in) :: mass,radius,tau
    real(dp) :: beta
    !rzhef: A function to evaluate Helium star radius on the MS

        beta = MAX(0.d0,0.4d0-0.22d0*LOG10(mass))
        radius_He_MS = radius *(1.d0+beta*(tau-tau**6))
    return
    end

    real(dp) FUNCTION radius_He_HG(m,lum,rx,lum2)
    
    implicit none
    real(dp):: m,lum,rx,cm,lum2

    !rhehgf: A function to evaluate Helium star radius on the Hertzsprung gap
    !from its mass and luminosity: see eqn 86 & 87 of SSE

        cm = 2.0d-03*m**(5.d0/2.d0)/(2.d0 + m**5)   !1/lambda
        radius_He_HG = rx*(lum/lum2)**0.2+0.02*(EXP(cm*lum)-EXP(cm*lum2))

    return
    end

    real(dp) FUNCTION radius_He_GB(lum)
    implicit none
    real(dp) lum
    !rhegbf: A function to evaluate Helium star radius on the giant branch.
        radius_He_GB = 0.08d0*lum**(3.d0/4.d0)
    return
    end


    !lum_GB(pars% age,t% zams_mass,t% He_pars% LtMS, t% MS_time)
    real(dp) FUNCTION lum_He_GB(age,Tinf1,Tinf2,Tx,D)
    implicit none
    real(dp), intent(in) :: age,Tinf1,Tinf2,Tx,D
       !lgbtf: A function to evaluate L given t for GB, AGB NHe star
        if(age <= Tx)then
         lum_He_GB = D*((p-1.0)*AHe*D*(Tinf1-age))**(p/(1.0-p))
        else
         lum_He_GB = B*((q-1.0)*AHe*B*(Tinf2-age))**(q/(1.0-q))
        endif
    return
    end

    !core_mass_He_GB(pars% luminosity )
    real(dp) FUNCTION core_mass_He_GB(lum,lx,D)
        implicit none
        real(dp),intent(in) :: lum,lx,D
        
       !mcgbf(lum,GB,lx): A function to evaluate Mc
    !given L for GB, AGB and NHe stars
          if(lum <= Lx)then
             core_mass_He_GB = (lum/D)**(1.d0/p)
          else
             core_mass_He_GB = (lum/B)**(1.d0/q)
          endif
    end

    real(dp) FUNCTION max_core_mass_he(mt,mass)
    implicit none
    real(dp),intent(in) :: mt, mass
    real(dp) :: mtc
    !A function to evaluate max allowed core mass for NHe stars
        mtc = MIN(mt,1.45d0*mt-0.31d0)
        if(mtc <= 0.d0) mtc = mt
        max_core_mass_he = MIN(mtc,MAX(M_ch,0.773d0*mass-0.35d0))
    end

    real(dp) FUNCTION He_GB_age(mc,Tinf1,Tinf2,D,Mx)
    implicit none
    real(dp),intent(in) :: mc,Tinf1,Tinf2,D,Mx

        if(mc.le.Mx)then
            He_GB_age = Tinf1 - (1.d0/((p-1.d0)*AHe*D))*(mc**(1.d0-p))
        else
            He_GB_age = Tinf2 - (1.d0/((q-1.d0)*AHe*B))*(mc**(1.d0-q))
        endif
    end

    real(dp) FUNCTION lmcgbf(mc,D,Mx)
    implicit none
    real(dp), intent(in) :: mc,D,Mx
        !lmcgbf: A function to evaluate L given Mc for GB, AGB and NHe stars

        if(mc.le.Mx)then
         lmcgbf = D*(mc**p)
        else
         lmcgbf = B*(mc**q)
        endif
      return
    end

!    subroutine calculate_mchei(zpars)
!        real(dp) :: zpars(20),GB(10)
!
!      lums(6) = GB(4)*GB(7)**GB(5)
!
!              zpars(10) = mcgbf(lums(4),GB,lums(6))
!
!
!          lums(4) = lHeIf(mass,zpars(2))
!
!
!    end


    subroutine calculate_rg(t,rg)
    !  rg = giant branch or Hayashi track radius, appropiate for the type.
    !       For kw=1 or 2 this is radius at BGB, and for kw=4 either GB or
    !       AGB radius at present luminosity.
    implicit none
    type(track), pointer, intent(in) :: t
    real(dp), intent(out):: rg
    real(dp) :: Rbgb, Rbagb, Lbgb, Lbagb, L, alfa
    integer :: j
    logical :: debug
    
        debug = .false.
        select case(t% pars% phase)
            case(low_mass_MS: HG)
                if (identified(BGB_EEP)) then
                    j = min(BGB_EEP,t% ntrack)
                    rg = t% tr(i_logR,j)   !TODO: check BGB  (or lums(3))
                else
                    j = min(cHeIgnition_EEP,t% ntrack)
                    rg = t% tr(i_logR,j)
                endif
                rg = 10.d0**rg
            case(RGB)
                rg = t% pars% radius
            case(HeBurn)
                !Linear interpolation between r(bgb) and r(bagb)
                !wrt luminosity l(bgb) and l(bagb) and L
!                if (identified(BGB_EEP)) then
!                    j = min(BGB_EEP,t% ntrack)
!                else
!                    j = min(cHeIgnition_EEP,t% ntrack)
!                endif
                j = min(cHeIgnition_EEP,t% ntrack)
                Rbgb = t% tr(i_logR,j)
                Lbgb = t% tr(i_logL,j)
                j = min(TA_cHeB_EEP,t% ntrack)
                Rbagb = t% tr(i_logR,j)
                Lbagb = t% tr(i_logL,j)
                L = log10(t% pars% luminosity)
                alfa = (L - Lbgb)/(Lbagb-Lbgb)
                rg = (alfa*Rbagb)+((1d0 - alfa)*Rbgb)
                rg = 10**rg
                
                if (rg .lt. t% pars% radius) then
                    write(UNIT=err_unit,fmt=*)'Error in calculate_rg: Rg, R',rg,t% pars% radius,Rbgb, Rbagb
                    write(UNIT=err_unit,fmt=*) t% pars% age, alfa, L, Lbgb, Lbagb
                endif
                    
            case(EAGB:TPAGB)
                rg = t% pars% radius
            case(He_MS)
                rg = radius_He_ZAMS(t% pars% mass)
            case(He_HG: He_GB)
                rg = radius_He_GB(t% pars% luminosity)
        end select
    end subroutine

end module sse_support
