subroutine dealloc_track()
    use track_support
    implicit none
    
    integer:: n,i

    n = size(tarr)
!    print*," deallocating", n
!    n = 1
    do i = 1,n
!        deallocate(tarr(i)% times_new)
!        deallocate(tarr(i)% times)
        deallocate(tarr(i)% eep)
        deallocate(tarr(i)% tr)
        deallocate(tarr(i)% cols)
        deallocate(tarr(i)% bounds)
        if((tarr(i)% ierr/=0).and.verbose)print*,'Error in evolving the system; check error file',i

    end do
    deallocate(tarr)
end subroutine dealloc_track
