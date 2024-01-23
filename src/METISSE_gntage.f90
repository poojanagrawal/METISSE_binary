subroutine METISSE_gntage(mc,mt,kw,zpars,m0,aj,id)

    use track_support
    use interp_support
    use z_support, only: Mcrit
    
    implicit none
    integer, intent(in), optional :: id
    real(dp) :: m0,aj,mt,mc
    real(dp):: tscls(20),lums(10),GB(10),zpars(20)
    real(dp):: dtm,mcy,tm,tn,mt0
    integer :: kw,idd,j
    type(track), pointer :: t
    logical :: debug
    
    debug = .false.
        if (debug) print*, 'in gntage',kw,m0,mt,mc,aj

    
    idd = 1
    if(present(id)) idd = id
    t => tarr(idd)
    
    
    
    dtm = 0.d0
    mcy = 0.d0
    mt0 = mt
    ! for very low mass stars MS extends beyond the age of the universe
    ! so higher phases are often missing
    ! below simplification avoids error in length
    if(mt<Mmin_array(TA_cHeB_EEP) .and. kw>1) then != very_low_mass_limit
        mt = Mmin_array(TA_cHeB_EEP)
!        print*,mt,mt0
    endif
    
    !this is just to signal star that gnatge is calling it
    !pars% phase will get updated to its correct value in star
!    if (t% pars% phase<=kw) t% pars% phase = kw+1
    t% pars% phase = -15

    !TODO: provide a backup in case one of the mcrits are not defined

    
    if(kw.eq.4)then
    ! Set the minimum CHeB core mass (at BGB or He ignition)
    ! using M = Mflash/Mhef
         if(Mcrit(4)% loc >0) then
             j = min(s(Mcrit(4)% loc)% ntrack,cHeIgnition_EEP)
             mcy = s(Mcrit(4)% loc)% tr(i_he_core,j)
         endif
         if(mc.le.mcy) then
            kw = 3
            if (debug) WRITE(*,*)' GNTAGE4: changed to 3'
         endif
         
      endif

    !Next we check that we don't have a GB star for M => Mfgb
      if(kw.eq.3)then
        ! Set the maximum GB core mass using M = Mfgb
        if(Mcrit(5)% loc >0) then
            j = min(s(Mcrit(5)% loc)% ntrack,cHeIgnition_EEP)
            mcy = s(Mcrit(5)% loc)% tr(i_he_core,j)
        endif
        if(mc.ge.mcy)then
            kw = 4
            aj = 0.d0
           if (debug) WRITE(*,*)' GNTAGE3: changed to 4'
         endif
      endif

    select case(kw)
    
    case(6)
        ! We try to start the star from the start of the SAGB
        ! by setting Mc = Mc,TP.
        t% pars% age = t% times(EAGB)
        CALL star(kw,m0,mt,tm,tn,tscls,lums,GB,zpars,dtm,id)
        aj = tscls(13)
            
    case(5)
        ! We fit a Helium core mass at the base of the AGB.
        t% pars% age = t% times(HeBurn)
        CALL star(kw,m0,mt,tm,tn,tscls,lums,GB,zpars,dtm,id)
        aj = tscls(2) + tscls(3)
    case(4)
        ! The supplied age is actually the fractional age, fage, of CHeB lifetime
        ! that has been completed, ie. 0 <= aj <= 1.
        if(aj.lt.0.d0.or.aj.gt.1.d0) aj = 0.d0
        t% pars% age = t% times(RGB) + aj*(t% times(HeBurn)-t% times(RGB))
        
        !for stars that don't have RGB phase, times(RGB) corresponds to times(HG)
        CALL star(kw,m0,mt,tm,tn,tscls,lums,GB,zpars,dtm,id)
        aj = tscls(2) + aj*tscls(3)
    case(3)
        !Place the star at the BGB
        t% pars% age = t% times(RGB)
        CALL star(kw,m0,mt,tm,tn,tscls,lums,GB,zpars,dtm,id)
        aj = tscls(1) + 1.0d-06*(tscls(2) - tscls(1))
!    case(9)

    end select
    
    t% pars% age = aj
    t% pars% phase = kw
    mt = mt0
    nullify(t)
    
    if (debug) print*, 'exit gntage',kw,m0,mt,mc,aj
end subroutine METISSE_gntage
