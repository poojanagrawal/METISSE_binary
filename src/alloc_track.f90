subroutine allocate_track(n,mass)
    use track_support
    implicit none

    integer, intent(in):: n
    real(dp), intent(in) :: mass(n)
    integer :: i
        
!    n= 1
!    n = size(mass)
!    print*,"I am in alloc_track with ", n,mass

    allocate(tarr(n))
    do i = 1,n
        tarr(i)% initial_mass = undefined
        tarr(i)% pars% extra = 0
        tarr(i)% pars% age_old = 0.d0
        tarr(i)% pars% mass = undefined
    !    tarr(i)% old_pars% mass = undefined
        tarr(i)% pars% delta = 0.d0
        tarr(i)% pars% bhspin = 0.d0
        tarr(i)% ierr = 0 
    end do
    
end subroutine allocate_track

