subroutine comenv_lambda(KW,M0,L,R,MENVD,LAMBDA,id,LAMBF)
!real(dp) FUNCTION
use track_support
implicit none
real(dp), intent(in):: M0,L,R,MENVD,LAMBDA
integer, intent(in) :: KW
integer, intent(in), optional :: id
real(dp), intent(out) :: LAMBF
real(dp):: RZAMS
integer :: idd
REAL(dp) :: CELAMF
EXTERNAL CELAMF


!if (SSE_FLAG.eqv..TRUE.) THEN
!  RZAMS = RZAMSF(M01)
!  LAMB1 = CELAMF(KW,M01,L1,R1,RZAMS,MENVD,LAMBDA,J1)
!ELSE
if(present(id))then
    idd = id
else
    idd = 1
endif
RZAMS = 10.d0**tarr(idd)% tr(i_logR,ZAMS_EEP)
LAMBF = CELAMF(KW,M0,L,R,RZAMS,MENVD,LAMBDA)
!comenv_lambda = LAMBF
end
