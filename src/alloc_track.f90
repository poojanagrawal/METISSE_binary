subroutine allocate_track(n,mass)
    use track_support
    implicit none

    integer, intent(in):: n
    real(dp), intent(in), optional :: mass(n)
        
!    n= 1
!    n = size(mass)
!    print*,"I am in alloc_track with ", n,mass

    allocate(tarr(n))
    tarr% star_type = unknown
    tarr% pars% age = 0.d0
end subroutine allocate_track

