      SUBROUTINE zcnsts(z,zpars)
      IMPLICIT NONE
      INCLUDE 'const_bse.h'
      
      real*8 z,zpars(20)
      
      if (SSE_FLAG.eqv..TRUE.) then
          !WRITE(*,*) 'Calling SSE_zcnsts'
          CALL SSE_zcnsts(z,zpars)
      else!if (METISSE_FLAG.eqv..TRUE.) then
          !WRITE(*,*) 'Calling METISSE_zcnsts',using_METISSE
          CALL METISSE_zcnsts(z,zpars)
      endif
      END
