subroutine dealloc_track()
    use track_support
    implicit none
    
    integer:: n,i

    n = size(tarr)
!    print*," deallocating", n
!    n = 1

    do i =1,n
        deallocate(tarr(i)% eep)
        deallocate(tarr(i)% cols)
        deallocate(tarr(i)% tr)
        deallocate(tarr(i)% phase)
        deallocate(tarr(i)% times)
    end do
    deallocate(tarr)
end subroutine dealloc_track
