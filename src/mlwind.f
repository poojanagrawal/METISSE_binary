      real*8 FUNCTION mlwind(kw,lum,r,mt,mc,rl,z,id)
      IMPLICIT NONE
      INCLUDE 'const_bse.h'
      
      integer kw,id
      real*8 lum,r,mt,mc,rl,z,dtm
      
      real*8 SSE_mlwind, METISSE_mlwind
      external SSE_mlwind, METISSE_mlwind
    
      if (SSE_FLAG.eqv..TRUE.) then
          !WRITE(*,*) 'Calling SSE_mlwind'
          mlwind = SSE_mlwind(kw,lum,r,mt,mc,rl,z)
      else!if (METISSE_FLAG.eqv..TRUE.) then
          !WRITE(*,*) 'Calling METISSE_mlwind'
          mlwind = METISSE_mlwind(kw,lum,r,mt,mc,rl,z,id)
      endif
      END
