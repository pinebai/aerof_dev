      SUBROUTINE BOUNDFLUX5(type, gam, normal, sig, ua, ub, flux)
c     
c     FLUX = A+(Wi)*Wi + A-(Winf)*Winf
c
      IMPLICIT NONE
c     num are the vertex numbers for the facet
c     ua is the state in primitive,
c     ub is the facet external boundary value passed in conservative variables
      REAL*8 ua(*), ub(*), vnx, vny, vnz
      REAL*8 normal(3), flux(*)
c     Sig is the dot product of the normal with the velocity of the facet
c     times the area of the facet
      REAL*8  sig
c     Plus and minus fluxes
      REAL*8 fgp(10), fgm(10)
      REAL*8 cap
      REAL*8 u, v, w, ro, usro, vit, pm, c2, c, usc
      REAL*8 vp4, vp5, vpp, vpm
      REAL*8 gam, gam1, rhoE, updir
      REAL*8 capd1, capd2, capd3, cc, unx, vp1, xn1, xn2, ff, gg, phiro
      INTEGER type
c
      gam1 = gam - 1.d0
c
      vnx = normal(1)
      vny = normal(2)
      vnz = normal(3)
c     
      cap		     = SQRT(vnx*vnx+ vny*vny + vnz*vnz)
      capd1		     = vnx/cap
      capd2		     = vny/cap
      capd3		     = vnz/cap
c
      ro		     = ua(1)
      usro		     = 1.0/ro     
c     u		     = ua(2,is)*usro
c     v		     = ua(3,is)*usro
c     w		     = ua(4,is)*usro     
      u                   = ua(2)
      v                   = ua(3)
      w                   = ua(4)
      vit		     = u*u + v*v + w*w
c     pm		     = gam1*(ua(5,is) - 0.5*ro*vit)
      pm                  = ua(5)         
      rhoE                = pm/gam1 + 0.5*ro*vit
c     
      c2		     = gam*pm*usro
      c		     = SQRT(c2)
      usc		     = 1.0/c
      cc		     = 0.5*vit + c2/gam1
c     
      unx		     = capd1*u + capd2*v + capd3*w
c     
      vp1		     = MAX(cap*unx - sig, 0.0d0)
      vp4		     = MAX(cap*(unx + c) - sig, 0.0d0)
      vp5		     = MAX(cap*(unx - c) - sig, 0.0d0)
c     
      vpp		     = vp4 + vp5
      vpm		     = vp4 - vp5
c     
      fgp(1) 	     = (gam1*ro*vp1 + 0.5*ro*vpp)/gam
      fgp(2) 	     = fgp(1)*u + 0.5*pm*usc*capd1*vpm
      fgp(3) 	     = fgp(1)*v + 0.5*pm*usc*capd2*vpm
      fgp(4) 	     = fgp(1)*w + 0.5*pm*usc*capd3*vpm
      fgp(5) 	     = 0.5*((gam1/gam)*ro*vit*
     &  			    (vp1 - 0.5*vpp) + 
     &  			    rhoE*vpp + 
     &  			    ro*unx*c*vpm/gam)

c
      vp1		     = MIN(cap*unx - sig, 0.0d0)
      vp4		     = MIN(cap*(unx + c) - sig, 0.0d0)
      vp5		     = MIN(cap*(unx - c) - sig, 0.0d0)
c
      vpp		     = vp1 - 0.5*(vp4 + vp5)
      vpm		     = 0.5*(vp5 - vp4)
c
      xn1		     = gam1*usc*(-0.5*vit*ub(1) + 
     &  		       u*ub(2)+v*ub(3)+w*ub(4)-ub(5))
      xn2		     = unx*ub(1) - 
     &  		       capd1*ub(2)-capd2*ub(3)-capd3*ub(4)
c
      ff		     = usc*(vpp*xn1 + vpm*xn2)
      gg		     = vpm*xn1 + vpp*xn2
c
      fgm(1) 	     = vp1*ub(1) + ff
      fgm(2) 	     = vp1*ub(2) + u*ff  + capd1*gg
      fgm(3) 	     = vp1*ub(3) + v*ff  + capd2*gg
      fgm(4) 	     = vp1*ub(4) + w*ff  + capd3*gg
      fgm(5) 	     = vp1*ub(5) + cc*ff + unx*gg
c
      phiro = fgp(1) + fgm(1)
c
      flux(1) = phiro 
      flux(2) = fgp(2) + fgm(2)
      flux(3) = fgp(3) + fgm(3)
      flux(4) = fgp(4) + fgm(4)
      flux(5) = fgp(5) + fgm(5)
c
      if (type.eq.1) then
         updir = 0.5d0 + dsign(0.5d0, phiro)
         usro = 1.0d0 / ub(1)
         flux(6) = phiro * (updir * ua(6) + (1.0d0-updir) * ub(6)*usro)
      else if (type.eq.2) then
         updir = 0.5d0 + dsign(0.5d0, phiro)
         usro = 1.0d0 / ub(1)
         flux(6) = phiro * (updir * ua(6) + (1.0d0-updir) * ub(6)*usro)
         flux(7) = phiro * (updir * ua(7) + (1.0d0-updir) * ub(7)*usro)
      endif
c
      END
