subroutine allocate_track(n,mass)
    use track_support
    implicit none

    integer, intent(in):: n
    real(dp), intent(in), optional :: mass(:)
        
!    n= 1
!    n = size(mass)
!    print*,"I am in alloc_track with ", n,mass

    allocate(tarr(n))
    tarr% star_type = unknown
    tarr% pars% age = 0.d0
    tarr% pars% extra = 0
    tarr% pars% bhspin = 0.d0
    tarr% ierr = 0
    tarr% pars% dms = 0.d0
    tarr% pars% delta = 0.d0
    

end subroutine allocate_track

