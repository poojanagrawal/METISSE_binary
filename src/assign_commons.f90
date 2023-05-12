subroutine assign_commons()
    use track_support,only: dp,pts_1,pts_2,pts_3
    use remnant_support, only:ns_flag,wd_flag,if_flag,Max_NS_mass,ec_flag
    implicit none
    
!   used to assign common variables when METISSE is used in standalone mode

REAL(dp) :: pts1,pts2,pts3
    COMMON /POINTS/ pts1,pts2,pts3
  
INTEGER :: ceflag,tflag,ifflag,nsflag,wdflag
COMMON /FLAGS/ ceflag,tflag,ifflag,nsflag,wdflag

REAL(dp) :: neta,bwind,hewind,mxns
      COMMON /VALUE1/ neta,bwind,hewind,mxns
    
ns_flag = nsflag
    wd_flag = wdflag
    if_flag = ifflag
    Max_NS_mass = mxns
    pts_1 = pts1
    pts_2 = pts2
    pts_3 = pts3
    ec_flag = 1
    !neta,bwind,hewind,sigma

    end subroutine
