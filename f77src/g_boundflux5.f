      SUBROUTINE GXBOUNDFLUX5(type, gam, normal, dnormal, sig, dsig,
     &                        ua, ub, dub, flux, dflux)
c     
c     FLUX = A+(Wi)*Wi + A-(Winf)*Winf
c
      IMPLICIT NONE
c     num are the vertex numbers for the facet
c     ua is the state in primitive,
c     ub is the facet external boundary value passed in conservative variables
      REAL*8 ua(*), ub(*), vnx, vny, vnz
      REAL*8 dua(7), dub(*), dvnx, dvny, dvnz
      REAL*8 normal(3), flux(*)
      REAL*8 dnormal(3), dflux(*)
c     Sig is the dot product of the normal with the velocity of the facet
c     times the area of the facet
      REAL*8  sig
      REAL*8  dsig
c     Plus and minus fluxes
      REAL*8 fgp(10), fgm(10)
      REAL*8 dfgp(10), dfgm(10)
      REAL*8 cap, capd1, capd2, capd3 
      REAL*8 dcap, dcapd1, dcapd2, dcapd3
      REAL*8 u, v, w, ro, usro, vit, pm, c2, c, usc
      REAL*8 du, dv, dw, dro, dusro, dvit, dpm, dc2, dc, dusc
      REAL*8 vp1, vp4, vp5, vpp, vpm
      REAL*8 dvp1, dvp4, dvp5, dvpp, dvpm
      REAL*8 cc, unx, xn1, xn2, ff, gg, phiro, rhoE 
      REAL*8 dcc, dunx, dxn1, dxn2, dff, dgg, dphiro, drhoE
      REAL*8 gam, gam1, updir
      INTEGER type
c
      gam1 = gam - 1.0d0
c
      vnx = normal(1)
      dvnx = dnormal(1)
      vny = normal(2)
      dvny = dnormal(2)
      vnz = normal(3)
      dvnz = dnormal(3)
c     
      cap = SQRT(vnx*vnx+ vny*vny + vnz*vnz)
      dcap = 1.0d0 / SQRT(vnx*vnx+ vny*vny + vnz*vnz) * 
     &       (vnx*dvnx+ vny*dvny + vnz*dvnz)
      capd1 = vnx/cap
      dcapd1 = (dvnx*cap-vnx*dcap)/(cap*cap)
      capd2 = vny/cap
      dcapd2 = (dvny*cap-vny*dcap)/(cap*cap)
      capd3 = vnz/cap
      dcapd3 = (dvnz*cap-vnz*dcap)/(cap*cap)
c
      dua(1) = 0.0d0
      dua(2) = 0.0d0
      dua(3) = 0.0d0
      dua(4) = 0.0d0
      dua(5) = 0.0d0
      dua(6) = 0.0d0
      dua(7) = 0.0d0
c
      ro = ua(1)
      dro = dua(1)
      usro = 1.0d0/ro     
      dusro = -1.0d0/(ro*ro) * dro
      u = ua(2)
      du = dua(2)
      v = ua(3)
      dv = dua(3)
      w = ua(4)
      dw = dua(4)
      vit = u*u + v*v + w*w
      dvit = 2.0d0*(u*du + v*dv + w*dw)
      pm = ua(5)         
      dpm = dua(5)         
      rhoE = pm/gam1 + 0.5d0*ro*vit
      drhoE = dpm/gam1 + 0.5d0*dro*vit + 0.5d0*ro*dvit
c     
      c2 = gam*pm*usro
      dc2 = gam*dpm*usro + gam*pm*dusro
      c = SQRT(c2)
      dc = 1.0d0/(2.0d0*SQRT(c2)) * dc2
      usc = 1.0d0/c
      dusc = -1.0d0/(c*c) * dc
      cc = 0.5*vit + c2/gam1
      dcc = 0.5*dvit + dc2/gam1
c     
      unx = capd1*u + capd2*v + capd3*w
      dunx = dcapd1*u + capd1*du + dcapd2*v + 
     &       capd2*dv + dcapd3*w + capd3*dw
c     
      vp1 = MAX(cap*unx - sig, 0.0d0)
      if (MAX(cap*unx - sig, 0.0d0) .eq. 0.0d0) then
        dvp1 = 0.0d0
      else
        dvp1 = dcap*unx + cap*dunx - dsig
      end if
      vp4 = MAX(cap*(unx + c) - sig, 0.0d0)
      if (MAX(cap*(unx + c) - sig, 0.0d0) .eq. 0.0d0) then
        dvp4 = 0.0d0
      else
        dvp4 = dcap*(unx + c) + cap*(dunx + dc) - dsig
      end if
      vp5 = MAX(cap*(unx - c) - sig, 0.0d0)
      if (MAX(cap*(unx - c) - sig, 0.0d0) .eq. 0.0d0) then
        dvp5 = 0.0d0
      else
        dvp5 = dcap*(unx - c) + cap*(dunx - dc) - dsig
      end if
c     
      vpp = vp4 + vp5
      dvpp = dvp4 + dvp5
      vpm = vp4 - vp5
      dvpm = dvp4 - dvp5
c     
      fgp(1) = (gam1*ro*vp1 + 0.5d0*ro*vpp)/gam
      dfgp(1) = (gam1*dro*vp1 + gam1*ro*dvp1 + 
     &          0.5d0*dro*vpp + 0.5d0*ro*dvpp)/gam
      fgp(2) = fgp(1)*u + 0.5d0*pm*usc*capd1*vpm
      dfgp(2) = dfgp(1)*u + fgp(1)*du + 0.5d0*dpm*usc*capd1*vpm +
     &          0.5d0*pm*dusc*capd1*vpm + 0.5d0*pm*usc*dcapd1*vpm +
     &          0.5d0*pm*usc*capd1*dvpm
      fgp(3) = fgp(1)*v + 0.5d0*pm*usc*capd2*vpm
      dfgp(3) = dfgp(1)*v + fgp(1)*dv + 0.5d0*dpm*usc*capd2*vpm + 
     &          0.5d0*pm*dusc*capd2*vpm + 0.5d0*pm*usc*dcapd2*vpm + 
     &          0.5d0*pm*usc*capd2*dvpm
      fgp(4) = fgp(1)*w + 0.5d0*pm*usc*capd3*vpm
      dfgp(4) = dfgp(1)*w + fgp(1)*dw + 0.5d0*dpm*usc*capd3*vpm + 
     &          0.5d0*pm*dusc*capd3*vpm + 0.5d0*pm*usc*dcapd3*vpm + 
     &          0.5d0*pm*usc*capd3*dvpm
      fgp(5) = 0.5d0*((gam1/gam)*ro*vit*
     &         (vp1 - 0.5d0*vpp) + 
     &         rhoE*vpp + 
     &         ro*unx*c*vpm/gam)
      dfgp(5) = 0.5d0*((gam1/gam)*dro*vit*(vp1 - 0.5d0*vpp) + 
     &          (gam1/gam)*ro*dvit*(vp1 - 0.5d0*vpp) +
     &          (gam1/gam)*ro*vit*(dvp1 - 0.5d0*dvpp) +
     &          drhoE*vpp + rhoE*dvpp + 
     &          dro*unx*c*vpm/gam + ro*dunx*c*vpm/gam + 
     &          ro*unx*dc*vpm/gam + ro*unx*c*dvpm/gam)

c
      vp1 = MIN(cap*unx - sig, 0.0d0)
      if (MIN(cap*unx - sig, 0.0d0) .eq. 0.0d0) then
        dvp1 = 0.0d0
      else
        dvp1 = dcap*unx + cap*dunx - dsig
      end if
      vp4 = MIN(cap*(unx + c) - sig, 0.0d0)
      if (MIN(cap*(unx + c) - sig, 0.0d0) .eq. 0.0d0) then
        dvp4 = 0.0d0
      else
        dvp4 = dcap*(unx + c) + cap*(dunx + dc) - dsig
      end if
      vp5 = MIN(cap*(unx - c) - sig, 0.0d0)
      if (MIN(cap*(unx - c) - sig, 0.0d0) .eq. 0.0d0) then
        dvp5 = 0.0d0
      else
        dvp5 = dcap*(unx - c) + cap*(dunx - dc) - dsig
      end if
c
      vpp = vp1 - 0.5*(vp4 + vp5)
      dvpp = dvp1 - 0.5*(dvp4 + dvp5)
      vpm = 0.5*(vp5 - vp4)
      dvpm = 0.5*(dvp5 - dvp4)
c
      xn1 = gam1*usc*(-0.5*vit*ub(1) + 
     &      u*ub(2) + v*ub(3) + w*ub(4) - ub(5))
      dxn1 = gam1*dusc*(-0.5*vit*ub(1) + 
     &       u*ub(2)+v*ub(3)+w*ub(4)-ub(5)) + gam1*usc*
     &       (-0.5*dvit*ub(1) - 0.5*vit*dub(1) + du*ub(2) + u*dub(2) + 
     &       dv*ub(3) + v*dub(3) + dw*ub(4) + w*dub(4) - dub(5))
      xn2 = unx*ub(1) - capd1*ub(2) - capd2*ub(3) - capd3*ub(4)
      dxn2 = dunx*ub(1) + unx*dub(1) - dcapd1*ub(2) - capd1*dub(2) - 
     &       dcapd2*ub(3) - capd2*dub(3) - dcapd3*ub(4) - capd3*dub(4)
c
      ff = usc*(vpp*xn1 + vpm*xn2)
      dff = dusc*(vpp*xn1 + vpm*xn2) + usc* 
     &      (dvpp*xn1 + vpp*dxn1 + dvpm*xn2 + vpm*dxn2)
      gg = vpm*xn1 + vpp*xn2
      dgg = dvpm*xn1 + vpm*dxn1 + dvpp*xn2 + vpp*dxn2
c
      fgm(1) = vp1*ub(1) + ff
      dfgm(1) = dvp1*ub(1) + vp1*dub(1) + dff
      fgm(2) = vp1*ub(2) + u*ff  + capd1*gg
      dfgm(2) = dvp1*ub(2) + vp1*dub(2) + du*ff + u*dff + 
     &          dcapd1*gg + capd1*dgg
      fgm(3) = vp1*ub(3) + v*ff  + capd2*gg
      dfgm(3) = dvp1*ub(3) + vp1*dub(3) + dv*ff + v*dff + 
     &          dcapd2*gg + capd2*dgg
      fgm(4) = vp1*ub(4) + w*ff  + capd3*gg
      dfgm(4) = dvp1*ub(4) + vp1*dub(4) + dw*ff + w*dff + 
     &          dcapd3*gg + capd3*dgg
      fgm(5) = vp1*ub(5) + cc*ff + unx*gg
      dfgm(5) = dvp1*ub(5) + vp1*dub(5) + dcc*ff + cc*dff + 
     &          dunx*gg + unx*dgg
c
      phiro = fgp(1) + fgm(1)
      dphiro = dfgp(1) + dfgm(1)
c
      flux(1) = phiro 
      dflux(1) = dphiro 
      flux(2) = fgp(2) + fgm(2)
      dflux(2) = dfgp(2) + dfgm(2)
      flux(3) = fgp(3) + fgm(3)
      dflux(3) = dfgp(3) + dfgm(3)
      flux(4) = fgp(4) + fgm(4)
      dflux(4) = dfgp(4) + dfgm(4)
      flux(5) = fgp(5) + fgm(5)
      dflux(5) = dfgp(5) + dfgm(5)
c
      if (type.eq.1) then
        updir = 0.5d0 + dsign(0.5d0, phiro)
        usro = 1.0d0 / ub(1)
        dusro = -1.0d0 / (ub(1)*ub(1)) * dub(1)
        flux(6) = phiro * (updir * ua(6) + (1.0d0-updir)*ub(6)*usro)
        dflux(6) = dphiro * (updir * ua(6) + (1.0d0-updir)*ub(6)*usro) +
     &  	   phiro * (updir * dua(6) + (1.0d0-updir)*dub(6)*usro +
     &  	   (1.0d0-updir)*ub(6)*dusro)
      else if (type.eq.2) then
        updir = 0.5d0 + dsign(0.5d0, phiro)
        usro = 1.0d0 / ub(1)
        dusro = -1.0d0 / (ub(1)*ub(1)) * dub(1)
        flux(6) = phiro * (updir * ua(6) + (1.0d0-updir)*ub(6)*usro)
        dflux(6) = dphiro * (updir * ua(6) + (1.0d0-updir)*ub(6)*usro) +
     &             phiro * (updir * dua(6) + (1.0d0-updir)*dub(6)*usro +
     &             (1.0d0-updir) * ub(6)*dusro)
        flux(7) = phiro * (updir * ua(7) + (1.0d0-updir)*ub(7)*usro)
        dflux(7) = dphiro * (updir * ua(7) + (1.0d0-updir)*ub(7)*usro) +
     &             phiro * (updir * dua(7) + (1.0d0-updir)*dub(7)*usro +
     &             (1.0d0-updir)* ub(7)*dusro)
      endif
c
      END
