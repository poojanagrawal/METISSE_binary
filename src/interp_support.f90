module interp_support

    use track_support
    use z_support, only: Mcrit, m_cutoff
    use sse_support, only: calculate_phases_and_times
    implicit none

    integer, parameter :: no_interpolation = 0
    integer, parameter :: linear = 1
    integer, parameter :: Steffen1990 = 2
    !integer, parameter :: extrapolation = 3
    logical :: debug_mass

    contains

    !   interpolates a new track given initial mass
    subroutine interpolate_mass(mass,t)
        implicit none
        real(dp), intent(in) :: mass
        type(track), pointer :: t
!        integer, intent(in), optional :: id
        
        integer :: iseg, m,keyword,min_index
        type(eep_track), pointer :: a(:)
        real(dp) :: f(3), dx, x(4), y(4), alfa, beta
        integer :: i, j, k, mlo, mhi,ierr
        integer, allocatable :: min_eeps(:)
        
    
        a => NULL()
!        t => tarr(1)
!        if(present(id)) t => tarr(id)
    
        debug_mass = .false.
        ierr=0

        !this line is to avoid array length problem while multiple calls to fix-track
        if (allocated(t% tr)) call deallocate_arrays(t)
        
        !   takes a set of EEP-defined tracks and find tracks for interpolation (a)
        call findtracks_for_interpolation(mass,a,iseg,m,min_index,keyword)

        !   interpolate the new track for given initial mass
        dx=0d0; alfa=0d0; beta=0d0; x=0d0; y=0d0
        mlo = 1
        mhi = size(a)

        k = minloc(a(mlo:mhi)% ntrack,dim=1)
        if (a(m)% initial_mass .ge. Mcrit(2)% mass) then
            t% ntrack = min(Final_eep, a(k)% ntrack)
        else
            t% ntrack = min(TAMS_EEP, a(k)% ntrack)
        endif
        
        call write_header(t,a(m))

        !	doing interpolation based on keyword
        select case(keyword)
        case(no_interpolation)
            t% tr(1:t% ncol,1:t% ntrack) = a(m)% tr(1:t% ncol,1:t% ntrack)
            ! t has an extra age column (i_age), so cannot use assumed shape arrays in the above
        case(linear)
            !print*, "case1: mlo mhi", mlo, mhi
            alfa = (t% initial_mass - a(mlo)% initial_mass)/(a(mhi)% initial_mass - a(mlo)% initial_mass)
            beta = 1d0 - alfa
            do i=1,t% ntrack
                do j=1,t% ncol
                    t% tr(j,i) = alfa*a(mhi)% tr(j,i) + beta*a(mlo)% tr(j,i)
                enddo
            enddo

        case(Steffen1990)
            x = a(mlo:mhi)% initial_mass
            dx = t% initial_mass - x(2)
            do i=1,t% ntrack
                do j=1,t% ncol
!                    if (j == age_col) cycle
                    do k=1,4
                        y(k) = a(k)% tr(j,i)
                    enddo

                    call interp_4pt_pm(x, y, f)
                    t% tr(j,i) = y(2) + dx*(f(1) + dx*(f(2) + dx*f(3)))
                enddo
!                alfa = dx/(x(3) - x(2))
!                beta = 1d0 - alfa
!                t% tr(age_col,i) = alfa*a(3)% tr(age_col,i) + beta*a(2)% tr(age_col,i)
            enddo
            
        end select
        

        if (fix_track) call fix_incomplete_tracks(iseg,t,min_index)

        ! get eeps
        
!        if (mass > Mcrit(2)% mass .and. mass< Mcrit(5)% mass) then
!        if (a(m)% initial_mass .ge. Mcrit(2)% mass) then
!            t% ntrack = min(Final_eep, a(k)% ntrack)
!            t% neep = min(size(key_eeps),a(k)% neep)
!        else
!            t% ntrack = min(TAMS_EEP, a(k)% ntrack)
!            t% neep = min(4,a(k)% neep)
!        endif
!        allocate(t% eep(t% neep))
!        t% eep = key_eeps(1:t% neep)  !
        
        allocate(min_eeps(size(key_eeps)))
        min_eeps= -1
        do i = 1, size(key_eeps)
            if (key_eeps(i)<= t% ntrack) min_eeps(i) = key_eeps(i)
!            if (identified(BGB_EEP) .and. key_eeps(i)= BGB_EEP) then
!                if (t% initial_mass >= Mcrit(5)% mass) min_eeps(i)-1
!            endif
        end do

        t% neep = count(min_eeps>0,1)

        allocate(t% eep(t% neep))
        t% eep = pack(min_eeps, min_eeps>0)
        if (debug_mass) print*, 'eeps',t% neep,t% eep(t% neep), key_eeps(size(key_eeps))
        
        !   check if mass and age are monotonic
        call smooth_track(t)

        !   recalibrate age from ZAMS and calculate timescales
        !   also assign SSE phases (Hurley et al.2000)
!        call calculate_phases_and_times(t)
        allocate(t% phase(t% ntrack))
        allocate(t% times(11))
        
        t% phase = undefined
        t% times = undefined
!   recalibrate age from ZAMS

        t% tr(i_age2,:) = t% tr(i_age2,:)- t% tr(i_age2,ZAMS_EEP)
        t% tr(i_age2,:) = t% tr(i_age2,:)*1E-6          !Myrs
        if(ierr/=0) then
        write(0,*) 'interpolate_track: interpolation failed for ', mass
        endif
 	nullify(a)
    end subroutine interpolate_mass

    subroutine findtracks_for_interpolation(mass,a,iseg,m,min_index,keyword)
        !takes a set of EEP-defined tracks and find tracks for interpolation
        real(dp), intent(in) :: mass
        integer, intent(out) :: iseg,m,min_index,keyword
        type(eep_track), pointer, intent(out) :: a(:)
        type(eep_track), pointer :: sa(:)
        real(dp) :: min_mass

        integer :: m_low,m_high,num_list
        real(dp), pointer :: mass_list(:)

        sa => NULL()
        a => NULL()
        m = 1
        m_low = 0
        num_list = 0
        
        !mcutoff array is defined in zpars.f90
        do iseg = 1,size(m_cutoff)
            m_low = m_cutoff(iseg)
            m_high = m_cutoff(iseg+1)-1
            ! to avoid using m_low twice- can cause error!!
            ! if it does, change the upper limit to size(m_cutoff)-1

            if (m_high == m_low) m_high = m_low+1  !TODO: - check this
            if (mass > s(m_high)% initial_mass*1.01) cycle
                num_list = m_high-m_low+1      !to include count for m_low point
                sa => s(m_low:m_high)
                allocate(mass_list(num_list))
                mass_list = sa% initial_mass
            exit
        end do
        
        !search database for nearest in mass track
        call index_search(num_list,mass_list,mass,min_index)
        
        if(min_index>num_list) min_index = num_list
        min_mass=s(m_low+min_index-1)% initial_mass

        if (debug_mass) print*,"min_index for mass", m_low+min_index-1,min_index
        if (debug_mass) print*,"mass(min_index)", min_mass
        
        if(abs(mass-sa(min_index)% initial_mass)< mass_accuracy_limit) then       !if track is already in the database
            keyword = no_interpolation
            !allocate(interpolate(1))
            !interpolate = min_index
            a => sa(min_index:min_index)
            if (debug_mass) print*,"Interpolation NOT required during interpolate mass",mass

        elseif(mass < sa(2)% initial_mass) then     !index adjustments for mass at boundaries
            keyword = linear
            !allocate(interpolate(2))
            !interpolate = [1,2]
            a => sa(1:2)
            if (debug_mass) print*, "Close to lower limit"

        elseif(mass > sa(num_list-1)% initial_mass) then
            keyword = linear
            !allocate(interpolate(2))
            !interpolate = [num_list-1,num_list]
            a => sa(num_list-1:num_list)
            if (debug_mass) print*, "Close to uppper limit"
        else
            keyword = Steffen1990
!            allocate(interpolate(4))
            if (sa(min_index)% initial_mass> mass) then
                !interpolate = (/(j, j=min_index-2,min_index+1)/)
                a => sa(min_index-2:min_index+1)
                m = 3
            else
                !interpolate = (/(j, j=min_index-1,min_index+2)/)
                a => sa(min_index-1:min_index+2)
                m = 2
            end if
        endif
        
        if (debug_mass) print*,"mass, keyword", mass,keyword
        if (debug_mass) print*,"interpolate mass" , a% initial_mass
        if (debug_mass) print*,"interpolate ntrack" , a% ntrack
        if (debug_mass .and. mass < sa(1)% initial_mass) print*, "doing extrapolation"
!        if (mass < sa(1)% initial_mass) print*, "doing extrapolation", mass, sa(1)% initial_mass, a% initial_mass
!        print*, mass, sa(1)% initial_mass

        nullify(sa)
		deallocate(mass_list)
    end subroutine findtracks_for_interpolation

    subroutine write_header(b,a)
        implicit none
        type(eep_track):: a
        type(track), pointer :: b
        
        b% star_type = a% star_type
        !b% version_string = a% version_string
        b% initial_Z = a% initial_Z
        b% initial_Y = a% initial_Y
        b% Fe_div_H = a% Fe_div_H
        b% alpha_div_Fe = a% alpha_div_Fe
        b% v_div_vcrit = a% v_div_vcrit
        b% ncol = a% ncol
        allocate(b% cols(b% ncol+1))
        b% cols(1: b% ncol) = a% cols
        b% cols(i_age2)% name = a% cols(i_age2)% name
         
        b% cols(b% ncol+1)% name = 'age_old'
!        if (.not. allocated(b% tr))
        allocate(b% tr(b% ncol+1, b% ntrack))
        !allocate(b% eep(b% neep))
        !b% eep = a% eep(1:b% neep)
        !if(a% has_phase) b% phase = a% phase
        b% tr = 0d0
        b% complete = .true.
        b% has_mass_loss = a% has_mass_loss

    end subroutine

    subroutine fix_incomplete_tracks(iseg,t,min_index)
        type(track), intent(inout) :: t
        integer, intent(in) :: min_index,iseg

        type(eep_track), pointer :: sa(:)
        real(dp), pointer :: mass_list(:)
        integer :: m_low,m_high,num_list
        integer :: i, m1, up_count,low_count,temp(4)
        integer :: min_ntrack,temp_ntrack, low_lim, upp_lim
        real(dp) :: upper_tol, lower_tol
        type(eep_track) :: a(2)
        real(dp), allocatable :: c(:,:)
        
        min_ntrack = get_min_ntrack(t% initial_mass, t% star_type)
        !check length
        if (debug_mass) print*,"checking length now"
        if ((t% ntrack >= min_ntrack) ) then
          if (debug_mass) print*,"length ok", t% initial_mass, t% ntrack
          return
        else
            if (debug_mass) print*,"not complete", t% initial_mass, t% ntrack, min_ntrack
            
            m_low = m_cutoff(iseg)
            m_high = m_cutoff(iseg+1)-1
            num_list = m_high-m_low+1      !to include count for m_low point
            sa => s(m_low:m_high)
            allocate(mass_list(num_list))
            mass_list = s(m_low:m_high)% initial_mass


            temp =0
            up_count =0; low_count =0

            upper_tol = t% initial_mass + t% initial_mass* lookup_index
            call index_search(num_list,mass_list,upper_tol,upp_lim)
            upp_lim= min(num_list,upp_lim)

            lower_tol = t% initial_mass - t% initial_mass* lookup_index
            call index_search(num_list,mass_list,lower_tol,low_lim)
            low_lim = max(1,low_lim)

            if (sa(min_index)% initial_mass> t% initial_mass) then
            m1 = min_index
            else
            m1 = min_index+1
            end if

            do i= m1-1, low_lim,-1
                if (sa(i)% ntrack >= min_ntrack) then
                    low_count = low_count+1
                    temp(low_count)=i        !new_a(count) = sa(i)
                if (low_count == 2) exit
                endif
            end do
            do i= m1, upp_lim
                if (sa(i)% ntrack >= min_ntrack) then
                    up_count = up_count+1
                    temp(low_count+up_count)=i        !new_a(count) = sa(i)
                    if (up_count == 2) exit
                endif
            end do
            
            select case(low_count)
            case(0)
                if (up_count==2) then  !extrapolate
                    !m =1
                    a = sa(pack(temp,mask = temp .ne. 0))
                    if (debug_mass) print*,"0,2"
                else
                    t% complete = .false.
                endif

            case(1)
                if(up_count>0) then !interpolate
                    temp(3)=0       !if already not so
                    a = sa(pack(temp,mask = temp .ne. 0))
                    if (debug_mass) print*,"1,1"
                else
                    t% complete = .false.
                endif

            case(2)
                if(up_count>0) then !interpolate
                    temp(2) = 0  !low_count =2 and we don't want to use more distant track
                    temp(4) = 0 !if already not so
                    a = sa(pack(temp,mask = temp .ne. 0))
                    if (debug_mass) print*,"2,1"
                else  !extrapolate
                    !m=1
                    temp(3) = temp(1) !linearly increasing order
                    temp(1) = 0
                    a = sa(pack(temp,mask = temp .ne. 0))
                    if (debug_mass) print*,"2,0"
                endif
            end select
            
            if (t% complete) then
                if (debug_mass) print*, "new masses for interpolate", a% initial_mass
                
                !store orginal track tr in c
                temp_ntrack = t% ntrack
                

!                print*, t% ntrack, size(t% tr, dim=2),t% initial_mass
                allocate(c(t% ncol, t% ntrack))
                c(:,:) = t% tr(1: t% ncol,:)

                !reallocate tr for rewriting with new length
                deallocate(t% tr)
                t% ntrack = min_ntrack

                allocate(t% tr(t% ncol+1, t% ntrack))

                t% tr = 0d0
                t% tr(1: t% ncol,1: temp_ntrack) = c(:,1:temp_ntrack)
                
                !complete the track between temp_ntrack and t% ntrack
                call linear_interp(a,t,temp_ntrack)
                deallocate(c,mass_list)
                nullify(mass_list,sa)

                if (debug_mass) print*, "new length", t% ntrack
            end if
        endif
        !t% tr = reshape(t%tr, (/t% ncol, t% ntrack/))
        
        
    end subroutine fix_incomplete_tracks
    
    
    integer function get_min_ntrack(initial_mass, star_type)
    real(dp) :: initial_mass
    integer:: star_type
        get_min_ntrack = 0
        !calculating min required length to the new track
        if (initial_mass >= Mcrit(2)% mass) then
            if (star_type == star_low_mass) then     !use Mec here?
                get_min_ntrack = min(low_mass_final_eep,final_eep)
            else if (star_type == star_high_mass)then
                get_min_ntrack = min(high_mass_final_eep,final_eep)
            endif
        else
            get_min_ntrack = TAMS_EEP    !a(m)% ntrack
        endif
    end function
    
    
    subroutine linear_interp(a,t,n)
    !this has been modified for use with fix_icomplete_tracks only
        implicit none
        type(eep_track), intent(in) :: a(:)
        type(track), intent(inout) :: t
        integer, intent(in) :: n
        real(dp) :: alfa,beta,bprime(t% ncol)
        integer :: i,j

        alfa=0d0; beta=0d0
        bprime = 0.d0
        alfa = (t% initial_mass - a(1)% initial_mass)/(a(2)% initial_mass - a(1)% initial_mass)
        beta = 1d0 - alfa

        !determining the offest from previously calculated value
        do j=1,t% ncol
            bprime(j) = t% tr(j,n)-(alfa*a(2)% tr(j,n) + beta*a(1)% tr(j,n))
        enddo

        do i= n, t% ntrack
            do j=1,t% ncol
                t% tr(j,i) = alfa*a(2)% tr(j,i) + beta*a(1)% tr(j,i) +bprime(j)
            enddo
        enddo
    end subroutine linear_interp

    subroutine smooth_track(t)
    implicit none
    type(track), intent(inout) :: t
    integer :: i

        call mod_PAV(t% tr(i_age2,:))
        
        do i=2,t% ntrack
            t% tr(i_mass,i) = min(t% tr(i_mass,i), t% tr(i_mass,i-1))
        end do
    end subroutine smooth_track
    
    subroutine mod_PAV(y)
        !PAV from ISO, modified for steps in time instead of smoothing over with average values
        real(dp), intent(inout) :: y(:)
        integer :: i,j, n, start, last, m,old_start
        real(dp), allocatable :: d(:)
        real(dp) :: diff, h
        logical :: debug

        debug = .false.
        n = size(y)
        allocate(d(n-1))
        old_start = 0

        do while(.true.)
           d = y(2:n)-y(1:n-1)
           if(all(d>=0)) exit !test for monotonicity
           i = locate(d) !finds the first point in d that is < 0
           start = i
            last = n
            if (debug) print*, 'in mod_PAV', i, d(i)

            if (start == old_start) start=old_start-1 !if start is the greatest value not in the end, take one step up and redo
            do while(.true.)
            if (y(n)- y(start)>0) exit
            start=start-1
            if (start <1)  then
                print*,"error in mod_PAV"
                stop
            endif
            end do
            if (debug) print*,'start, old_start', start, old_start
            do j= start+1,n
                diff= y(j)-y(start)
                if (debug) print*,"i and diff",i,diff
                if (diff>0) then
                    last = j
                    if (debug) print*,j,y(j),y(i)
                    if (y(last)>y(n)) last = n
                    exit
                endif
            enddo

           m = last - start
           h = (y(last)-y(start))/real(m)
            do j= 1,m
                y(j+start) = y(start)+ j*h
            end do
            old_start = start
            if (debug) print*,last,start,y(last),y(start),h,m
        end do
        
        deallocate(d)
    end subroutine mod_PAV
  
    integer function locate(y)
    real(dp), intent(in) :: y(:)
    integer :: i, n
    n=size(y)
    do i=1,n
       if(y(i)<0.0)then
          locate=i
          return
       endif
    end do
    locate=0
    end function locate
    
    subroutine interpolate_age(t, input_age, icolumn, val)
        implicit none
        real(dp), intent(in) :: input_age


        integer, intent(in), optional :: icolumn
        real(dp), intent(out), optional :: val
        integer :: jstart, jend, age_col, kw
        real(dp) ::  them, them_new, frac, age, age2


        real(dp) :: f(3), dx, x(4), y(4), alfa, beta
        integer :: j, k, mlo, mhi, pass, n_pass
        
        real(dp), allocatable :: new_line(:,:)
        integer, allocatable :: min_eeps(:)

        logical :: debug, interpolate
        type(track), pointer :: t

        debug = .false.

!        print*,"in  interpolate age"
        jstart = 1
        jend = t% ncol
        if (present(icolumn)) then
            jstart = icolumn
            jend = icolumn
        endif
        allocate (new_line(t% ncol,1))
        new_line = -1.d0
        dx=0d0; alfa=0d0; beta=0d0; x=0d0; y=0d0

        kw = t% pars% phase
        !TODO: this is temporary, to avoid NaN during interpolation
        ! it happens due to gntage
        if (kw >1) then
            if(t% times(kw)-t% times(kw-1)<1d-12) kw = kw+1
        endif
        if (kw>1 .and. kw<6) then
            !scale the input age for the new track
            !TODO: check tscls(kw) and (kw-1) are >0 (defined)
            them = t% times(kw)-t% times(kw-1)
            them_new = t% times_new(kw)-t% times_new(kw-1)
            frac = (input_age-t% times(kw-1))/them
            age2 = t% times_new(kw-1)+(frac*them_new)
!            age2 = new_age(t% times(kw),t% times(kw-1),t% times_new(kw),t% times_new(kw-1),input_age)
            
            n_pass = 2
            if (t% irecord>0 .and. debug) print*, "in interp2", age2, t% pars% age,kw
    !print*,t% times(kw),t% times(kw-1),t% times_new(kw),t% times_new(kw-1)
        !t% pars% mass,t% pars% core_mass,t% times_new(kw)
        else
            age2 = input_age
            n_pass = 1
!            if (kw ==1) age2= input_age*(t% MS_time/t% ms_old)
!            print*, 'age MS', age2
        endif

        do pass = 1, n_pass
            if (pass == 1) then
                age = age2
                age_col = i_age2     !i_age2 = new age, age col in the main array
            else
                age = input_age
                age_col = i_age     !i_age = old age, stored at t% ncol+1
            endif

    !        print*, age,age_col
            call find_nearest_eeps(t,min_eeps, age, age_col)
            mlo = minval(min_eeps)
            mhi = maxval(min_eeps)

            if (debug) print*, 'pass', n_pass,pass, age, kw
            if (debug) print*,"neighbouring_eeps", min_eeps

            if (mhi == mlo) then
                if (t% irecord>0 .and. debug) print*, "no interp in age needed"
                do j=jstart,jend
                    interpolate = check_core_quant(j,n_pass, pass)
                    if (interpolate) new_line(j,1) = t% tr(j,mlo)
                end do

            elseif ((mhi-mlo)<4) then
                !linear interpolation
                if (debug) print*, "doing linear interp in age"
                alfa = (age - t% tr(age_col,mlo))/(t% tr(age_col,mhi) - t% tr(age_col,mlo))
                beta = 1d0 - alfa
                do j=jstart,jend
                    interpolate = check_core_quant(j,n_pass, pass)
                    if (interpolate) then
                        new_line(j,1) = alfa*t% tr(j,mhi) + beta*t% tr(j,mlo)
                        if (new_line(j,1)/= new_line(j,1)) then
                        print*, '**Warning: NaN encountered during interpolation age** ',t% initial_mass,input_age,j,mhi,mlo
    !                    stop
                        endif
                    endif
                end do
                if (debug) print*, "ending linear interp in age"

            else
                if (debug) print*, "doing cubic interp in age"

    !            mlo = 1!min(max(1,m-1),n-3)  !1
    !            mhi = 4!max(min(m+2,n),4)     !4

                x = t% tr(age_col,mlo:mhi)
                dx = new_line(age_col,1) - x(2)

                if (age< x(2) .or. age> x(3)) then
                   print*,"Error in age interpolation in interp_support"
                   stop
                endif

                do j=jstart,jend
                    interpolate = check_core_quant(j,n_pass, pass)
                    if (interpolate) then
                        do k=1,4
                            y(k) = t% tr(j,mlo-1+k)
                        enddo
                        call interp_4pt_pm(x, y, f)
                        new_line(j,1) = y(2) + dx*(f(1) + dx*(f(2) + dx*f(3)))
                    end if
                enddo
            endif
            
        end do
        if (present(icolumn)) then
            val = new_line(icolumn,1)
        else
            call save_values(new_line,t% pars)
        endif

        if (t% pars% mass <0.0) then
            print*,"fatal error: mass <0 "
            stop
        endif
        deallocate(new_line,min_eeps)
         
        if (debug) print*, 'exiting interpolate_age'
    end subroutine interpolate_age
    
    real(dp) function new_age(tc,tprev,tnew,tnew_prev,age)
        real(dp), intent(in) :: tc,tprev,tnew,tnew_prev,age
        real(dp) ::  frac
    
!    tc = t% times(kw)
!    tprev = t% times(kw-1)
!    tnew = t% times_new(kw)
!    tnew_prev = t% times_new(kw-1)
        
        frac = (age-tprev)/(tc-tprev)
        new_age = tnew_prev+(frac*(tnew-tnew_prev))
    end function
            
    
    logical function check_core_quant(j,n_pass, pass)
        integer :: j, pass, n_pass
        check_core_quant = .true.

        if (any(j .eq. core_cols,1) .or. j == i_age2) then
            if (n_pass>1 .and. pass==1) check_core_quant = .false.    !do not interpolate in core quantities
            check_core_quant = .true.
        else
            check_core_quant = .true.
        endif
        !print*, j, n_pass, pass,check_core_quant
        return
    end function check_core_quant


    subroutine find_nearest_eeps(t,min_eeps,age,age_col)
        implicit none
        type(track), pointer :: t
        integer, allocatable, intent(out) :: min_eeps(:)
        real(dp), intent(in) :: age
        integer :: i, j, len_eep, min_index, age_col
        real(dp), allocatable :: age_list(:)
        real(dp) :: last_age
        logical :: debug
        
        debug =  .false.
        !Todo: min_eeps-> nbr_eeps

        last_age = 0.d0
        do i = 1,t% neep-1
            if (debug) print*,"loc_low", t% eep(i), t% eep(i+1),t%neep

            last_age = t% tr(age_col,t% eep(i+1))
            if (debug) print*,"ages", age,last_age, t% eep(i)< initial_eep

            if ((t% eep(i)< initial_eep) .or. (age .gt. last_age)) cycle
            len_eep = t% eep(i+1)-t% eep(i)+1
            allocate(age_list(len_eep))
            age_list = t% tr(age_col,t% eep(i):t% eep(i+1))

            if (debug) print*,"len_eep:",len_eep,"bounds:",age_list(1),age_list(len_eep)

            call index_search(len_eep,age_list,age,min_index)
            !min_index = binary_search(len_eep,age_list,age)
!                    if (debug) print*,"in interp_support, age, min_index, age at min index"
!                    if (debug) print*,t% pars% age, min_index,age_list(min_index)

            if(abs(age_list(min_index)-age)< tiny) then        !less than a year
                    if (debug) print*,"no interpolation, min_index", age_list(min_index)
                    allocate(min_eeps(1))
                    min_eeps = min_index
!                            a => b(:,min_index: min_index)
            elseif(age< age_list(2)) then
                    if (debug) print*,"age< age_list(2)", age_list(1:2)
                    allocate(min_eeps(2))
                    min_eeps = [1,2]
!                            a => b(:,1:2)
            elseif(age> age_list(len_eep-1)) then
                    if (debug) print*,"age> age_list(len_eep-1)", age_list(len_eep-1:len_eep)
                    allocate(min_eeps(2))
                    min_eeps = [len_eep-1,len_eep]
!                            a => b(:,len_eep-1:len_eep)
            elseif(age < age_list(min_index)) then
                    if (debug) print*,"age< min_index", age_list(min_index-2:min_index+1)

!                            a => b(:,min_index-2:min_index+1)
                    allocate(min_eeps(2))
                    min_eeps = (/(j, j=min_index-1,min_index)/)
!                            a => b(:,min_index-1:min_index) !forcing linear

            else
                    if (debug) print*,"age> min_index",age_list(min_index-1:min_index+2)
!                            a => b(:,min_index-1:min_index+2)
                    allocate(min_eeps(2))
                    min_eeps = (/(j, j=min_index,min_index+1)/)
!                            a => b(:,min_index:min_index+1)
            endif
            min_eeps = min_eeps +t% eep(i)-1    !original min_eeps were only a given primary eep, this is to scale them back to full track
            deallocate(age_list)
            exit
        end do
        if (.not. allocated(min_eeps)) print*, 'error finding nearest eeps'
    end subroutine find_nearest_eeps


    !from MESA-r7503/1d_interp/
    subroutine interp_4pt_pm(x, y, a)
        ! returns coefficients for monotonic cubic interpolation from x(2) to x(3)
        real(dp), intent(in)    :: x(4)    ! junction points, strictly monotonic
        real(dp), intent(in)    :: y(4)    ! data values at x's
        real(dp), intent(inout)   :: a(3)    ! coefficients
        real(dp) :: h1, h2, h3, s1, s2, s3, p2, p3, as2, ss2, yp2, yp3

        !integer, parameter :: pm_work_size = 3  !from mesa/interp_1d_def.f90

        ! for x(2) <= x <= x(3) and dx = x-x(2),
        ! y(x) = y(2) + dx*(a(1) + dx*(a(2) + dx*a(3)))
        h1 = x(2)-x(1)
        h2 = x(3)-x(2)
        h3 = x(4)-x(3)
        s1 = (y(2)-y(1))/h1
        s2 = (y(3)-y(2))/h2
        s3 = (y(4)-y(3))/h3
        p2 = (s1*h2+s2*h1)/(h1+h2)
        p3 = (s2*h3+s3*h2)/(h2+h3)
        as2 = abs(s2)
        ss2 = sign(1d0, s2)
        yp2 = (sign(1d0, s1)+ss2)*min(abs(s1), as2, 0.5d0*abs(p2))
        yp3 = (ss2+sign(1d0, s3))*min(as2, abs(s3), 0.5d0*abs(p3))
        a(1) = yp2
        a(2) = (3*s2-2*yp2-yp3)/h2
        a(3) = (yp2+yp3-2*s2)/(h2*h2)
    end subroutine interp_4pt_pm
    

    subroutine get_initial_mass_for_new_track(t, delta,interpolate_all)

    real(dp), intent(in) :: delta
    type(track), pointer :: t

    integer :: min_index,num_list,Mupp,Mlow,i,j,nt,kw
    real(dp), allocatable :: mlist(:),mlist1(:), age_list(:)
    real(dp) :: Mnew,alfa,beta,them, themold, age
    integer :: eep_m, eep_core

    logical :: debug, interpolate_all

    debug = .false.
    interpolate_all = .false.

    !nt is the length of the track before new interpolation
    eep_m = -1
    nt = t% ntrack
    allocate(age_list(nt))
    age_list = t% tr(i_age,:)

    kw = t% pars% phase
    if (t% pars% phase>1 .and. .false.) then
        !scale the age for the new track
        !TODO: check tscls(kw) and (kw-1) are >0 (defined)
        themold = t% times(kw)-t% times(kw-1)
        them = t% times_new(kw)-t% times_new(kw-1)
        age = t% times_new(kw-1)+(t% pars% age-t% times(kw-1))*(them/themold)
    else
    age = t% pars% age
    endif

    if (kw>4 .and.debug) print*,"clear0",age,age_list(1),age_list(t% ntrack)

    call index_search(nt,age_list, age,eep_m)

    if (kw>4 .and.debug) print*,"clear1: eep_m, nt :",eep_m,nt,kw

    !if (age_list(eep_m)<age) eep_m = eep_m+1
    if (eep_m > nt) eep_m = nt

    if (kw>4 .and.debug) print*,"clear2",age,age_list(eep_m),age_list(1), age_list(nt)


!    Mnew = t% tr(i_mass,eep_m)-delta  !TEST
    Mnew = t% pars% mass-delta  !TEST

    !print*,"mnew, eepm",Mnew,t% tr(i_mass,eep_m),eep_m

    
    num_list = count(s% ntrack >= eep_m,1)
    allocate(mlist(num_list),mlist1(num_list))
    ! NOTE: cannot use s(:)% tr(i_mass,eep_m), gives error
    j = 1
    do i =1, size(s)
        if (s(i)% ntrack >= eep_m) then
        mlist1(j) = s(i)% initial_mass
        mlist(j) = s(i)% tr(i_mass,eep_m)   !this eep_m can vary with fractional age
        ! only do until mlist(j)< mnew
        j= j+1
        endif
    end do

    ! get nearest mass for eep_m
    call index_search (num_list,mlist,Mnew,min_index)
    if (min_index > num_list) min_index = num_list

    if (kw<10 .and. debug)  print*,"clear2, min index mass", mlist(min_index)
    if (mlist(min_index)<Mnew) then
        Mlow = min_index
        Mupp = min_index+1
    else
        Mlow = min_index-1
        Mupp = min_index
    endif
    if(Mlow<1 .or. Mupp > num_list) then
        if (debug) print*,"Error: beyond the bounds for interpolation"
        if (debug) print*, "Mlow,Mupp,num_list,mnew,eep_m,kw", &
                Mlow,Mupp,num_list,mnew,eep_m,kw
        if (Mlow<1) Mlow = 1
        if (Mupp> num_list) Mupp = num_list
        return
    endif
    if (debug)  print*, "Mup =", mlist(Mupp), "mlow", mlist(Mlow),"mnew",mnew
    alfa = (Mnew - mlist(Mlow))/(mlist(Mupp) - mlist(Mlow))
    beta = 1d0 - alfa
    t% initial_mass = alfa*mlist1(Mupp) + beta*mlist1(Mlow)
    if (debug) print*, "new ini mass",t% initial_mass, mlist1(Mupp), mlist1(Mlow), eep_m


    if (kw<=1)then
        interpolate_all = .true.
!        eep_core = 0.75* TAMS_EEP
!        if(identified(IAMS_EEP)) eep_core = IAMS_EEP
!        if (eep_m>=eep_core) interpolate_all = .false.
!        print*, 'eep_core', eep_core, IAMS_EEP,interpolate_all
!TODO: need to scale age as well in the age interpolation
    endif
    deallocate(mlist,mlist1)
    deallocate(age_list)

    end subroutine
    end module interp_support
