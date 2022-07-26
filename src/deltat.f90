    subroutine deltat(kw,age,tm,tn,tscls,dt,dtr,id)
    
    !calculates timestep for evolution
    use track_support
    integer, intent(in), optional :: id

    INTEGER :: kw
    REAL(dp) :: age,tm,tn,tscls(20)
    REAL(dp) :: pts1,pts2,pts3
    COMMON /POINTS/ pts1,pts2,pts3
    
    real(dp) :: dtr,dt, tMS_hook
    type(track), pointer :: t => NULL()
    
    if(present(id))then
        t=> tarr(id)
    else
        t=> tarr(1)
    endif
    
        dtr = -1.d0
        dt = -1.d0
        tMS_hook = (0.95* t% times(1))
        !Base new time scale for changes in radius & mass on stellar type.
        select case (t% pars% phase)
            case(low_mass_MS:MS)
                if (age< tMS_hook) then
                    dt = pts1* t% times(1)
                else
                    dt = pts2*(t% times(1)-tMS_hook)
                endif
                    dtr = t% times(1)-age
            case(HG)
                dt = pts1*(t% times(2)-t% times(1))
!                dt = pts1*(t% times(2)-tm)
                dtr = t% times(2)-age
!                dt = pts1*(tscls(1) - tm)
!         dtr = tscls(1) - age
            case(RGB)
                dt = pts2*(t% times(3)- t% times(2))
                dtr = MIN(t% nuc_time,t% times(3))-age
!                print*,"dt", t% times(2),t% times(3),age

            case(HeBurn)
                dt = pts2*(t% times(4)- t% times(3))
                dtr = MIN(t% nuc_time,t% times(4))-age
!                print*, t% times(3),t% times(4),t% nuc_time
!                print*, dt,dtr,age
            case(EAGB)
                dt = pts3*(t% times(5)- t% times(4))
                dtr = MIN(t% nuc_time,t% times(5))-age
!                print*, t% nuc_time, t% times(5), age
            case(TPAGB)
                if (t% post_agb) then
                    dt = pts3*(t% times(6) - t% agb% age)
                    dtr = t% times(6) - age
!                    print*,"aj,tn",age,t% nuc_time,t% times(TPAGB),t% agb% age
                !t% times(6) gets recalculated for post_agb phase- see remnant_support.f
                else
                    dt = pts3*(t% times(6)-t% times(5))
                    dt = MIN(dt,0.001d0)
                    dtr = MIN(t% nuc_time,t% times(6))- age
                endif
            case(He_MS)
                dt = pts1* t% times(7)
                dtr = t% times(7) - age
                !this gets modified if star is losing too much mass
            case(He_HG:He_GB)
                if(age < t% times(10))then
                    dt = pts2*(t% times(8) - age)
                else
                    dt = pts2*(t% times(9) - age)
                endif
                dtr = t% nuc_time -age
            case(HeWD:NS)
                dt = MAX(0.1d0,age*10.0)
                dt = MAX(0.1d0,dt*10.0)
                !dt = MIN(dt,5.0E+02)
                dt = MIN(dt,1.d0)
                dtr = dt
            case(BH:Massless_REM)
                !dt = MIN(dt,1.d0)
                dt = 10.d0
                dtr = dt
        end select

            t% pars% dt = min(dt,dtr)
            
!           if (kw==5) print*,"dt dtr phase",dt,dtr,t% pars%phase
!            if (kw==5) print*, t% nuc_time, tn
            if (dtr<=0.0) then
                print*,"fatal error: invalid timestep", dtr ,"for phase", t% pars% phase
                print*,"t% zams_mass, t% nuc_time, age, t% pars% mass"
                print*,t% zams_mass, t% nuc_time, age, t% pars% mass
                stop
            endif
            
            nullify(t)
      end subroutine deltat
