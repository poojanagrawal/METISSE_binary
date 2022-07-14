subroutine dealloc_track()
    use track_support
    implicit none
    
    !for 1 star, make separate file for 2 stars dealloc_track2.f ?
!    real (dp) , intent(in) :: mass
    integer:: n,i
        
!    n = size(mass)
    n= 2

do i =1,n
        deallocate(tarr(i)% eep)
        deallocate(tarr(i)% cols)
        deallocate(tarr(i)% tr)
        deallocate(tarr(i)% phase)
        deallocate(tarr(i)% times)
end do
deallocate(tarr)
    end subroutine dealloc_track
