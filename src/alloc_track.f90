subroutine allocate_track(mass)
    use track_support
    implicit none

!for 1 star, make separate file for 2 stars alloc_track2.f ?
    real (dp) , intent(in) :: mass(2)
    integer:: n,i
        

    n = size(mass)
            print*,"I am in alloc_track with ", mass,n

!    n= 1
    allocate(tarr(n))
    do i =1,n

    tarr(i)% initial_mass = undefined
    tarr(i)% pars% extra = 0
    tarr(i)% pars% age_old = 0.d0
    tarr(i)% pars% mass = undefined
!    tarr(i)% old_pars% mass = undefined
    tarr(i)% pars% delta = 0.d0

    end do
    
    
end subroutine allocate_track

