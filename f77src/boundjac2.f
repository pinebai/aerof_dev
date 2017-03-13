      SUBROUTINE BOUNDJAC2(type, gam, normal, sig, ua, ub, jac)
c     
c     FLUX = A+(Wi)*Wi + A-(Wi)*Winf
c
      IMPLICIT NONE
c     num are the vertex numbers for the facet
c     ua is the state in primitive,
c     ub is the facet external boundary value passed in conservative variables
      REAL*8 ua(*), ub(*), vnx, vny, vnz
      REAL*8 normal(3), jac(*), jac1
c     Sig is the dot product of the normal with the velocity of the facet
c     times the area of the facet
      REAL*8  sig
c     Plus and minus fluxes
c$$$      REAL*8 fgp(10), fgm(10)
      REAL*8 fgp1, fgm1
      REAL*8 cap
      REAL*8 u, v, w, ro, usro, vit, pm, c2, c, usc
      REAL*8 vp4, vp5, vpp, vpm
      REAL*8 gam, gam1, updir
      REAL*8 capd1, capd2, capd3, unx, vp1, xn1, xn2, ff, phiro
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
c$$$      rhoE                = pm/gam1 + 0.5*ro*vit
c     
      c2		     = gam*pm*usro
      c		     = SQRT(c2)
      usc		     = 1.0/c
c$$$      cc		     = 0.5*vit + c2/gam1
c     
      unx		     = capd1*u + capd2*v + capd3*w
c     
      vp1		     = MAX(cap*unx - sig, 0.0d0)
      vp4		     = MAX(cap*(unx + c) - sig, 0.0d0)
      vp5		     = MAX(cap*(unx - c) - sig, 0.0d0)
c     
      vpp		     = vp4 + vp5
c$$$      vpm		     = vp4 - vp5
c     
      fgp1 	     = (gam1*ro*vp1 + 0.5*ro*vpp)/gam
c$$$      fgp(2) 	     = fgp(1)*u + 0.5*pm*usc*capd1*vpm
c$$$      fgp(3) 	     = fgp(1)*v + 0.5*pm*usc*capd2*vpm
c$$$      fgp(4) 	     = fgp(1)*w + 0.5*pm*usc*capd3*vpm
c$$$      fgp(5) 	     = 0.5*((gam1/gam)*ro*vit*
c$$$     &  			    (vp1 - 0.5*vpp) + 
c$$$     &  			    rhoE*vpp + 
c$$$     &  			    ro*unx*c*vpm/gam)

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
c$$$      gg		     = vpm*xn1 + vpp*xn2
c
      fgm1 	     = vp1*ub(1) + ff
c$$$      fgm(2) 	     = vp1*ub(2) + u*ff  + capd1*gg
c$$$      fgm(3) 	     = vp1*ub(3) + v*ff  + capd2*gg
c$$$      fgm(4) 	     = vp1*ub(4) + w*ff  + capd3*gg
c$$$      fgm(5) 	     = vp1*ub(5) + cc*ff + unx*gg
c
      phiro = fgp1 + fgm1
c
c$$$      flux(1) = phiro 
c$$$      flux(2) = fgp(2) + fgm(2)
c$$$      flux(3) = fgp(3) + fgm(3)
c$$$      flux(4) = fgp(4) + fgm(4)
c$$$      flux(5) = fgp(5) + fgm(5)
c
      updir = 0.5d0 + dsign(0.5d0, phiro)

      jac1 = phiro * updir * usro

      if (type.eq.1) then
         jac(1) = jac1
      else if (type.eq.2) then
         jac(1) = jac1
         jac(2) = 0.0
         jac(3) = 0.0
         jac(4) = jac(1)
      endif
c
      END
