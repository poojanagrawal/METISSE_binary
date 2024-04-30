      SUBROUTINE zcnsts(z,zpars,path_to_tracks,path_to_he_tracks)
      IMPLICIT NONE
      INCLUDE 'const_bse.h'
      
      real*8 z,zpars(20)
      CHARACTER*256 path_to_tracks,path_to_he_tracks

      if (SSE_FLAG.eqv..TRUE.) then
          !WRITE(*,*) 'Calling SSE_zcnsts'
          CALL SSE_zcnsts(z,zpars)
      else!if (METISSE_FLAG.eqv..TRUE.) then
          !WRITE(*,*) 'Calling METISSE_zcnsts',using_METISSE
          
          !SSE_zcnsts also sets some coefficients used in gntage and in the evolution of He stars
          !updating gntage and he stars for METISSE will remove the need to call this
          CALL SSE_zcnsts(z,zpars)
          CALL METISSE_zcnsts(z,zpars,path_to_tracks,path_to_he_tracks)
          
      endif
      END
