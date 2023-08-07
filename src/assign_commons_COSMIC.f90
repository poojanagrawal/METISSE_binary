subroutine assign_commons()
    use track_support
    implicit none
    
    !to assign common variables when METISSE is used with COSMIC
      
    REAL(dp) :: ecsn,ecsn_mlow
    COMMON /SNVARS1/ ecsn,ecsn_mlow
     
    real(dp) :: d
    
    
    if(front_end == COSMIC) then
    ! use inputs from COSMIC
    
        if (Mec_core > 0.d0) ecsn = Mec_core
        d = (Mec_core-Mup_core)
        if (Mup_core > 0.d0 .and. d>tiny ) ecsn_mlow = Mup_core
        print*, 'assigning commons cosmic', Mup_core,Mec_core,ecsn, ecsn_mlow
        
    else
        print*,'Error: Front end mismtach in assign commons'
        print*,'expected 2 (COSMIC); got ', front_end
    endif

    end subroutine

