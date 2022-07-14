***
      SUBROUTINE kick(kw,m1,m1n,m2,ecc,sep,jorb,vs)
      implicit none
*
* There are various choices that can be made for kicks. 
*
* For regular NSs the kick choice is controlled by the value of 
* sigma set in the input file. Choices are: 
*    sigma < 0.0 - kick drawn randomly between 0 - ABS(sigma) km/s
*    sigma = 0.0 - no kick
*    sigma > 0.0 - kick drawn from a Maxwellian with dispersion sigma. 
*
* Then for an electron capture supernova or an accretion-induced 
* collapse the choice is determined by the value of ecsig set 
* internally here. Choices are: 
*    ecsig = 0.0 - no kick
*    ecsig > 0.0 - kick drawn from a Maxwellian, dispersion ecsig
*    ecsig < 0.0 - same as for the regular NSs but scaled by ABS(ECSIG). 
* These supernova are identified by their mass of 1.26 Msun. 
*
* For BHs the kick choice is controlled by the value of bhflag 
* set in the input file. Choices are: 
*    bhflag = 0 - no kick
*    bhflag = 1 - same as for the regular NSs
*    bhflag = 2 - same as for the regular NSs but scaled by fallback. 
*
* Small kicks for WDs can also be set. These are determined by the 
* values set internally here for: 
*   wdsig1 - He and COWDs
*   wdsig2 - ONeWDs 
* where the kick will be zero if this value is zero, otherwise it will 
* be the dispersion in a Maxwellian. 
* Currently a limit of wdkmax = 6 km/s is set. 
*
      integer kw,k
      INTEGER idum
      COMMON /VALUE3/ idum
      INTEGER idum2,iy,ir(32)
      COMMON /RAND3/ idum2,iy,ir
      integer bhflag
      real*8 m1,m2,m1n,ecc,sep,jorb,ecc2
      real*8 pi,twopi,gmrkm,yearsc,rsunkm
      parameter(yearsc=3.1557d+07,rsunkm=6.96d+05)
      real*8 mm,em,dif,der,del,r
      real*8 u1,u2,vk,v(4),s,theta,phi
      real*8 sphi,cphi,stheta,ctheta,salpha,calpha
      real*8 vr,vr2,vk2,vk0,vn2,hn2
      real*8 mu,cmu,vs(3),v1,v2,mx1,mx2
      real*8 sigma,sigma0,ecsig,wdsig1,wdsig2,wdkmax
      COMMON /VALUE4/ sigma,bhflag
      real ran3,xx
      external ran3
      logical iflat
*
      ecsig = 20.d0
      wdsig1 = 0.d0
      wdsig2 = 0.d0
      wdkmax = 6.d0
*
      do k = 1,3
         vs(k) = 0.d0
      enddo
*
      pi = ACOS(-1.d0)
      twopi = 2.d0*pi
* Conversion factor to ensure velocities are in km/s using mass and
* radius in solar units.
      gmrkm = 1.906125d+05
*
* Find the initial separation by randomly choosing a mean anomaly.
      if(sep.gt.0.d0.and.ecc.ge.0.d0)then
         xx = RAN3(idum)
         mm = xx*twopi
         em = mm
 2       dif = em - ecc*SIN(em) - mm
         if(ABS(dif/mm).le.1.0d-04) goto 3
         der = 1.d0 - ecc*COS(em)
         del = dif/der
         em = em - del
         goto 2
 3       continue
         r = sep*(1.d0 - ecc*COS(em))
*
* Find the initial relative velocity vector.
         salpha = SQRT((sep*sep*(1.d0-ecc*ecc))/(r*(2.d0*sep-r)))
         calpha = (-1.d0*ecc*SIN(em))/SQRT(1.d0-ecc*ecc*COS(em)*COS(em))
         vr2 = gmrkm*(m1+m2)*(2.d0/r - 1.d0/sep)
         vr = SQRT(vr2)
      else
         vr = 0.d0
         vr2 = 0.d0
         salpha = 0.d0
         calpha = 0.d0
      endif
*
      iflat = .false.
      if(sigma.lt.-0.01)then
         if(kw.eq.13.and.m1n.ge.1.28) iflat = .true.
         if(kw.eq.13.and.m1n.lt.1.28.and.ecsig.lt.-0.01) iflat = .true.
         if(kw.eq.14.and.bhflag.gt.0) iflat = .true.
      endif
*
      if(iflat)then
*
* Generate the kick velocity from a flat distribution between 0 
* and ABS(sigma) km/s. 
*
         sigma0 = ABS(sigma)
         if(kw.eq.13.and.m1n.lt.1.28) sigma0 = sigma0*ABS(ecsig)
         vk = RAN3(idum)*sigma0
         theta = RAN3(idum)*twopi
         sphi = RAN3(idum)
         phi = ASIN(sphi)
         cphi = COS(phi)
         stheta = SIN(theta)
         ctheta = COS(theta)
         v(1) = ctheta*cphi*vk
         v(2) = stheta*cphi*vk
         v(3) = sphi*vk
         vk2 = vk*vk
      else
*
* Generate Kick Velocity using Maxwellian Distribution (Phinney 1992).
* Use Henon's method for pairwise components (Douglas Heggie 22/5/97).
*
         sigma0 = MAX(sigma,0.d0)
         if(kw.eq.10.or.kw.eq.11) sigma0 = MAX(wdsig1,0.d0)
         if(kw.eq.12) sigma0 = MAX(wdsig2,0.d0)
         if(kw.eq.13.and.m1n.lt.1.28)then
            if(ecsig.lt.-0.01)then
               sigma0 = sigma0*ABS(ecsig)
            else
               sigma0 = MAX(ecsig,0.d0)
            endif
         endif
         if(kw.eq.14.and.bhflag.eq.0) sigma0 = 0.d0
*
         do 20 k = 1,2
            u1 = RAN3(idum)
            u2 = RAN3(idum)
* Generate two velocities from polar coordinates S & THETA.
            s = sigma0*SQRT(-2.d0*LOG(1.d0 - u1))
            theta = twopi*u2
            v(2*k-1) = s*COS(theta)
            v(2*k) = s*SIN(theta)
 20      continue
*
         if(sigma0.gt.0.001)then
            vk2 = v(1)*v(1) + v(2)*v(2) + v(3)*v(3)
            vk = SQRT(vk2)
         else
            vk2 = 0.d0
            vk = 0.d0
         endif
*
         u1 = RAN3(idum)
         sphi = -1.d0 + 2.d0*u1
         phi = ASIN(sphi)
         cphi = COS(phi)
         stheta = SIN(theta)
         ctheta = COS(theta)
*
      endif
*
      vk0 = vk
*
* Impose the maximum WD kick velocity. 
*
      if(kw.le.12.and.vk.gt.wdkmax)then
         vk = wdkmax
         do k = 1,3
            vs(k) = vs(k)*vk/vk0
         enddo
         vk2 = vk*vk
      endif
*
* Restrict the BH kick velocity by fallback. 
* This could be done better but in the N-body code we only have 
* limited information. 
*
      if(kw.eq.14.and.bhflag.gt.1)then
         vk = vk0*(m1-m1n)/m1
         do k = 1,3
            vs(k) = vs(k)*vk/vk0
         enddo
         vk2 = vk*vk
      endif
*      write(66,*)' VK ',vk0,vk,m1,m1n
*
*     WRITE(66,*)' KICK VK PHI THETA ',vk,phi,theta
*
      if(sep.le.0.d0.or.ecc.lt.0.d0) goto 90
*
* Determine the magnitude of the new relative velocity.
      vn2 = vk2+vr2-2.d0*vk*vr*(ctheta*cphi*salpha-stheta*cphi*calpha)
* Calculate the new semi-major axis.
      sep = 2.d0/r - vn2/(gmrkm*(m1n+m2))
      sep = 1.d0/sep
*     if(sep.le.0.d0)then
*        ecc = 1.1d0
*        goto 90
*     endif
* Determine the magnitude of the cross product of the separation vector
* and the new relative velocity.
      v1 = vk2*sphi*sphi
      v2 = (vk*ctheta*cphi-vr*salpha)**2
      hn2 = r*r*(v1 + v2)
* Calculate the new eccentricity.
      ecc2 = 1.d0 - hn2/(gmrkm*sep*(m1n+m2))
      ecc2 = MAX(ecc2,0.d0)
      ecc = SQRT(ecc2)
* Calculate the new orbital angular momentum taking care to convert
* hn to units of Rsun^2/yr.
      jorb = (m1n*m2/(m1n+m2))*SQRT(hn2)*(yearsc/rsunkm)
* Determine the angle between the new and old orbital angular
* momentum vectors.
      cmu = (vr*salpha-vk*ctheta*cphi)/SQRT(v1 + v2)
      mu = ACOS(cmu)
*
 90   continue
*
      if(ecc.le.1.0)then
* Calculate the components of the velocity of the new centre-of-mass.
         mx1 = vk*m1n/(m1n+m2)
         mx2 = vr*(m1-m1n)*m2/((m1n+m2)*(m1+m2))
         vs(1) = mx1*ctheta*cphi + mx2*salpha
         vs(2) = mx1*stheta*cphi + mx2*calpha
         vs(3) = mx1*sphi
*        write(*,*)' vr vk vs ',vr,vk,mx2
      else
* Calculate the relative hyperbolic velocity at infinity.
         sep = r/(ecc-1.d0)
*        cmu = SQRT(ecc-1.d0)
*        mu = ATAN(cmu)
         mu = ACOS(1.d0/ecc)
         vr2 = gmrkm*(m1n+m2)/sep
         vr = SQRT(vr2)
         vs(1) = vr*SIN(mu)
         vs(2) = vr*COS(mu)
         vs(3) = 0.d0
         ecc = MIN(ecc,99.99d0)
      endif
*
 95   continue
*
      RETURN
      END
***
