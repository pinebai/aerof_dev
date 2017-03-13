      SUBROUTINE COMPUTEDBOUNDFLUXOPERATOR(type, gam, normal, 
     &                                     sig, 
     &                                     ua, ub, 
     &                                     dfluxdnormal, 
     &                                     dfluxdsig,
     &                                     dfluxdub)
c     
c     FLUX = A+(Wi)*Wi + A-(Winf)*Winf
c
      IMPLICIT NONE
c     num are the vertex numbers for the facet
c     ua is the state in primitive,
c     ub is the facet external boundary value passed in conservative variables
      REAL*8 ua(*), ub(*), vnx, vny, vnz
      REAL*8 normal(3)
      REAL*8 dfluxdnormal(3,7), dfluxdsig(7), dfluxdub(7,7)
      REAL*8 dfluxdvp1(7), dfluxdvp4(7), dfluxdvp5(7)
      REAL*8 dfluxdvp10(7), dfluxdvp40(7), dfluxdvp50(7)
      REAL*8 dvp1dnormal(3), dvp4dnormal(3), dvp5dnormal(3)
      REAL*8 dvp10dnormal(3), dvp40dnormal(3), dvp50dnormal(3)
      REAL*8 dvp1dsig, dvp4dsig, dvp5dsig
      REAL*8 dvp10dsig, dvp40dsig, dvp50dsig

c     Sig is the dot product of the normal with the velocity of the facet
c     times the area of the facet
      REAL*8  sig
c     Plus and minus fluxes
      REAL*8 fgp(10), fgm(10)
      REAL*8 dfgp(10), dfgm(10)
      REAL*8 cap, capd1, capd2, capd3 
      REAL*8 u, v, w, ro, usro, vit, pm, c2, c, usc
      REAL*8 dro
      REAL*8 vp1, vp4, vp5, vpp, vpm
      REAL*8 cc, unx, xn1, xn2, ff, gg, phiro, rhoE 
      REAL*8 dphiro
      REAL*8 gam, gam1, updir
      REAL*8 cf(200)
      INTEGER type, I, J
      DO 1010, J= 1,7
        dfluxdsig(J) = 0.0d0
        dfluxdvp10(J) = 0.0d0
        dfluxdvp40(J) = 0.0d0
        dfluxdvp50(J) = 0.0d0 
        dfluxdvp1(J) = 0.0d0
        dfluxdvp4(J) = 0.0d0
        dfluxdvp5(J) = 0.0d0 
        DO 1000, I = 1,3
          dfluxdnormal(I,J) = 0.0d0
 1000 CONTINUE
        DO 1020, I = 1,7
          dfluxdub(I,J) = 0.0d0
 1020 CONTINUE
 1010 CONTINUE
      DO 1030, I=1,3
        dvp1dnormal(I) = 0.0d0
        dvp4dnormal(I) = 0.0d0
        dvp5dnormal(I) = 0.0d0
        dvp10dnormal(I) = 0.0d0
        dvp40dnormal(I) = 0.0d0
        dvp50dnormal(I) = 0.0d0
 1030 CONTINUE
      dvp10dsig = 0.0d0
      dvp40dsig = 0.0d0
      dvp50dsig = 0.0d0
      dvp1dsig = 0.0d0
      dvp4dsig = 0.0d0
      dvp5dsig = 0.0d0
c
      gam1 = gam - 1.0d0
c
      vnx = normal(1)
      vny = normal(2)
      vnz = normal(3)
c     
      cap = SQRT(vnx*vnx+ vny*vny + vnz*vnz)
      cf(1) = 1.0d0 / SQRT(vnx*vnx+ vny*vny + vnz*vnz)
      cf(2) = cf(1)*vnx
      cf(3) = cf(1)*vny
      cf(4) = cf(1)*vnz
      capd1 = vnx/cap
      cf(5) = 1.0d0/(cap*cap)
      cf(6) = cf(5)*(cap - vnx*cf(2))
      cf(7) = -cf(5)*vnx*cf(3)
      cf(8) = -cf(5)*vnx*cf(4)
      capd2 = vny/cap
      cf(9) = -cf(5)*vny*cf(2)
      cf(10) = cf(5)*(cap-vny*cf(3))
      cf(11) = - cf(5)*vny*cf(4)
      capd3 = vnz/cap
      cf(12) = -cf(5)*vnz*cf(2)
      cf(13) = -cf(5)*vnz*cf(3)
      cf(14) = cf(5)*(cap - vnz*cf(4))
c
      ro = ua(1)
      usro = 1.0d0/ro     
      u = ua(2)
      v = ua(3)
      w = ua(4)
      vit = u*u + v*v + w*w
      pm = ua(5)         
      rhoE = pm/gam1 + 0.5d0*ro*vit
c     
      c2 = gam*pm*usro
      c = SQRT(c2)
      usc = 1.0d0/c
      cc = 0.5*vit + c2/gam1
c     
      unx = capd1*u + capd2*v + capd3*w
      cf(15) = (u*cf(6)+v*cf(9)+w*cf(12))
      cf(16) = (u*cf(7)+w*cf(13)+v*cf(10))
      cf(17) = (u*cf(8)+v*cf(11)+w*cf(14))
c     
      vp1 = MAX(cap*unx - sig, 0.0d0)
      if (MAX(cap*unx - sig, 0.0d0) .eq. 0.0d0) then
      else
        cf(18) = (unx*cf(2) + cap*cf(15))
        cf(19) = (unx*cf(3) + cap*cf(16))
        cf(20) = (unx*cf(4) + cap*cf(17))
        dvp10dnormal(1) = cf(18)
        dvp10dnormal(2) = cf(19)
        dvp10dnormal(3) = cf(20)
        dvp10dsig = -1.0d0
      end if
      vp4 = MAX(cap*(unx + c) - sig, 0.0d0)
      if (MAX(cap*(unx + c) - sig, 0.0d0) .eq. 0.0d0) then
      else
        cf(21) = ((unx + c)*cf(2) + cap*cf(15))
        cf(22) = ((unx + c)*cf(3) + cap*cf(16))
        cf(23) = ((unx + c)*cf(4) + cap*cf(17))
        dvp40dnormal(1) = cf(21)
        dvp40dnormal(2) = cf(22)
        dvp40dnormal(3) = cf(23)
        dvp40dsig = -1.0d0
      end if
      vp5 = MAX(cap*(unx - c) - sig, 0.0d0)
      if (MAX(cap*(unx - c) - sig, 0.0d0) .eq. 0.0d0) then
      else
        cf(24) = ((unx - c)*cf(2) + cap*cf(15))
        cf(25) = ((unx - c)*cf(3) + cap*cf(16))
        cf(26) = ((unx - c)*cf(4) + cap*cf(17))
        dvp50dnormal(1) = cf(24)
        dvp50dnormal(2) = cf(25)
        dvp50dnormal(3) = cf(26)
        dvp50dsig = -1.0d0
      end if
c     
      vpp = vp4 + vp5
      vpm = vp4 - vp5
c     
      fgp(1) = (gam1*ro*vp1 + 0.5d0*ro*vpp)/gam
      cf(49) = gam1*ro/gam
      cf(50) = 0.5d0*ro/gam
      fgp(2) = fgp(1)*u + 0.5d0*pm*usc*capd1*vpm
      cf(27) = 0.5d0*pm*usc*vpm
      cf(51) = (u*cf(50) + 0.5d0*pm*usc*capd1)
      cf(52) = (u*cf(50) - 0.5d0*pm*usc*capd1)
      cf(53) = cf(27)*cf(6)
      cf(54) = cf(27)*cf(7)
      cf(55) = cf(27)*cf(8)
      fgp(3) = fgp(1)*v + 0.5d0*pm*usc*capd2*vpm
      cf(28) = 0.5d0*pm*usc*vpm
      cf(56) = v*cf(49)
      cf(57) = (v*cf(50)+0.5d0*pm*usc*capd2)
      cf(58) = (v*cf(50)-0.5d0*pm*usc*capd2)
      cf(59) = cf(28)*cf(9)
      cf(60) = cf(28)*cf(10)
      cf(61) = cf(28)*cf(11) 
      fgp(4) = fgp(1)*w + 0.5d0*pm*usc*capd3*vpm
      cf(29) = 0.5d0*pm*usc*vpm
      cf(62) = w*cf(49)
      cf(63) = (w*cf(50)+0.5d0*pm*usc*capd3)
      cf(64) = (w*cf(50)-0.5d0*pm*usc*capd3)
      cf(65) = cf(29)*cf(12)
      cf(66) = cf(29)*cf(13)
      cf(67) = cf(29)*cf(14)
      fgp(5) = 0.5d0*((gam1/gam)*ro*vit*
     &         (vp1 - 0.5d0*vpp) + 
     &         rhoE*vpp + 
     &         ro*unx*c*vpm/gam)
      cf(30) = ro*c*vpm/gam
      cf(31) = (gam1/gam)*ro*vit
      cf(68) = 0.5d0*cf(31)
      cf(69) = (0.5d0*(rhoE-cf(31)*0.5d0)+0.5d0*ro*unx*c/gam)
      cf(70) = (0.5d0*(rhoE-cf(31)*0.5d0)-0.5d0*ro*unx*c/gam)
      cf(71) = 0.5d0*cf(30)*cf(15)
      cf(72) = 0.5d0*cf(30)*cf(16)
      cf(73) = 0.5d0*cf(30)*cf(17)
c
      vp1 = MIN(cap*unx - sig, 0.0d0)
      if (MIN(cap*unx - sig, 0.0d0) .eq. 0.0d0) then
      else
        cf(32) = (unx*cf(2)+cap*cf(15))
        cf(33) = (unx*cf(3)+cap*cf(16))
        cf(34) = (unx*cf(4)+cap*cf(17))
        dvp1dnormal(1) = cf(32)
        dvp1dnormal(2) = cf(33)
        dvp1dnormal(3) = cf(34)
        dvp1dsig = -1.0d0
      end if
      vp4 = MIN(cap*(unx + c) - sig, 0.0d0)
      if (MIN(cap*(unx + c) - sig, 0.0d0) .eq. 0.0d0) then
      else
        cf(35) = ((unx + c)*cf(2) + cap*cf(15))
        cf(36) = ((unx + c)*cf(3) + cap*cf(16))
        cf(37) = ((unx + c)*cf(4) + cap*cf(17)) 
        dvp4dnormal(1) = cf(35)
        dvp4dnormal(2) = cf(36)
        dvp4dnormal(3) = cf(37)
        dvp4dsig = -1.0d0
      end if
      vp5 = MIN(cap*(unx - c) - sig, 0.0d0)
      if (MIN(cap*(unx - c) - sig, 0.0d0) .eq. 0.0d0) then
      else
        cf(38) = ((unx - c)*cf(2)+cap*cf(15))
        cf(39) = ((unx - c)*cf(3)+cap*cf(16))
        cf(40) = ((unx - c)*cf(4)+cap*cf(17)) 
        dvp5dnormal(1) = cf(38)
        dvp5dnormal(2) = cf(39)
        dvp5dnormal(3) = cf(40)
        dvp5dsig = -1.0d0
      end if
c
      vpp = vp1 - 0.5*(vp4 + vp5)
      vpm = 0.5*(vp5 - vp4)
c
      xn1 = gam1*usc*(-0.5*vit*ub(1) + 
     &      u*ub(2) + v*ub(3) + w*ub(4) - ub(5))
      cf(41) = -gam1*usc*0.5*vit
      cf(42) = gam1*usc*u
      cf(43) = gam1*usc*v
      cf(44) = gam1*usc*w
      cf(45) = -gam1*usc
      xn2 = unx*ub(1) - capd1*ub(2) - capd2*ub(3) - capd3*ub(4)
      cf(46)=(ub(1)*cf(15)-ub(2)*cf(6)-ub(3)*cf(9)-ub(4)*cf(12))
      cf(47)=(ub(1)*cf(16)-ub(2)*cf(7)-ub(3)*cf(10)-ub(4)*cf(13))
      cf(48)=(ub(1)*cf(17)-ub(2)*cf(8)-ub(3)*cf(11)-ub(4)*cf(14)) 
c
      ff = usc*(vpp*xn1 + vpm*xn2)
      cf(74) = usc*xn1
      cf(75) = -(usc*xn1*0.5+usc*xn2*0.5)
      cf(76) = (usc*xn2*0.5-usc*xn1*0.5)
      cf(77) = (usc*vpm*unx+usc*vpp*cf(41))
      cf(78) = (usc*vpp*cf(42)-usc*vpm*capd1)
      cf(79) = (usc*vpp*cf(43)-usc*vpm*capd2)
      cf(80) = (usc*vpp*cf(44)-usc*vpm*capd3)
      cf(81) = usc*vpp*cf(45)
      cf(82) = usc*vpm*cf(46)
      cf(83) = usc*vpm*cf(47)
      cf(84) = usc*vpm*cf(48)
      gg = vpm*xn1 + vpp*xn2
      cf(85) = - 0.5*xn1
      cf(86) = - 0.5*xn2
      cf(87) = (cf(85)+cf(86))
      cf(88) = (cf(86)-cf(85))
      cf(89) = (vpm*cf(41)+vpp*unx)
      cf(90) = (vpm*cf(42)-vpp*capd1)
      cf(91) = (vpm*cf(43)-vpp*capd2)
      cf(92) = (vpm*cf(44)-vpp*capd3)
      cf(93) = vpm*cf(45)
      cf(94) = vpp*cf(46)
      cf(95) = vpp*cf(47)
      cf(96) = vpp*cf(48)
c
      fgm(1) = vp1*ub(1) + ff
      fgm(2) = vp1*ub(2) + u*ff  + capd1*gg
      cf(97) = (ub(2)*ub(2)+u*cf(74)+capd1*xn2)
      cf(98) = (u*cf(75)+capd1*cf(87))
      cf(99) = (u*cf(76)+capd1*cf(88))
      cf(100) = (u*cf(77)+capd1*cf(89))
      cf(101) = (vp1+capd1*cf(90)+u*cf(78))
      cf(102) = (u*cf(79)+capd1*cf(91))
      cf(103) = (u*cf(80)+capd1*cf(92))
      cf(104) = (capd1*cf(93)+u*cf(81))
      cf(105) = (u*cf(82)+capd1*cf(94)+gg*cf(6))
      cf(106) = (u*cf(83)+capd1*cf(95)+gg*cf(7))
      cf(107) = (capd1*cf(96)+u*cf(84)+gg*cf(8))
c    
      fgm(3) = vp1*ub(3) + v*ff  + capd2*gg
      cf(108) = (ub(3)+v*cf(74)+capd2*xn2)
      cf(109) = (capd2*cf(87)+v*cf(75))
      cf(110) = (capd2*cf(88)+v*cf(76))
      cf(111) = (v*cf(77)+capd2*cf(89))
      cf(112) = (v*cf(78)+capd2*cf(90))
      cf(113) = (vp1+v*cf(79)+capd2*cf(91))
      cf(114) = (v*cf(80)+capd2*cf(92))
      cf(115) = (v*cf(81)+capd2*cf(93))
      cf(116) = (v*cf(82)+gg*cf(9)+capd2*cf(94))
      cf(117) = (v*cf(83)+gg*cf(10)+capd2*cf(95))
      cf(118) = (v*cf(84)+gg*cf(11)+capd2*cf(96))
c
      fgm(4) = vp1*ub(4) + w*ff  + capd3*gg
      cf(119) = (ub(4)+w*cf(74)+capd3*xn2)
      cf(120) = (capd3*cf(87)+w*cf(75))
      cf(121) = (capd3*cf(88)+w*cf(76))
      cf(122) = (w*cf(77)+capd3*cf(89))
      cf(123) = (w*cf(78)+capd3*cf(90))
      cf(124) = (w*cf(79)+capd3*cf(91))
      cf(125) = (w*cf(80)+vp1+capd3*cf(92))
      cf(126) = (w*cf(81)+capd3*cf(93))
      cf(127) = (w*cf(82)+gg*cf(12)+capd3*cf(94))
      cf(128) = (capd3*cf(95)+gg*cf(13)+w*cf(83))
      cf(129) = (capd3*cf(96)+gg*cf(14)+w*cf(84))
c
      fgm(5) = vp1*ub(5) + cc*ff + unx*gg
      cf(130) = (ub(5)+unx*xn2+cc*cf(74))
      cf(131) = (unx*cf(87)+cc*cf(75))
      cf(132) = (unx*cf(88)+cc*cf(76))
      cf(133) = (cc*cf(77)+unx*cf(89))
      cf(134) = (unx*cf(90)+cc*cf(78))
      cf(135) = (unx*cf(91)+cc*cf(79))
      cf(136) = (unx*cf(92)+cc*cf(80))
      cf(137) = (unx*cf(93)+cc*cf(81)+vp1)
      cf(138) = (cc*cf(82)+gg*cf(15)+unx*cf(94))
      cf(139) = (cc*cf(83)+gg*cf(16)+unx*cf(95))
      cf(140) = (cc*cf(84)+gg*cf(17)+unx*cf(96))
c
      phiro = fgp(1) + fgm(1)
      cf(141) = (ub(1)+cf(74))
      cf(142) = (vp1+cf(77))
c
c      dflux(1) = cf(49)*dvp10
c     &        + cf(50)*dvp40
c     &        + cf(50)*dvp50
c     &        + cf(141)*dvp1
c     &        + cf(75)*dvp4
c     &        + cf(76)*dvp5
c     &        + cf(142)*dub(1)
c     &        + cf(78)*dub(2)
c     &        + cf(79)*dub(3)
c     &        + cf(80)*dub(4)
c     &        + cf(81)*dub(5)
c     &        + cf(82)*dnormal(1)
c     &        + cf(83)*dnormal(2)
c     &        + cf(84)*dnormal(3)
      dfluxdvp10(1) = cf(49)
      dfluxdvp40(1) = cf(50)
      dfluxdvp50(1) = cf(50)
      dfluxdvp1(1) = cf(141)
      dfluxdvp4(1) = cf(75)
      dfluxdvp5(1) = cf(76)
      dfluxdub(1,1) = cf(142)
      dfluxdub(2,1) = cf(78)
      dfluxdub(3,1) = cf(79)
      dfluxdub(4,1) = cf(80)
      dfluxdub(5,1) = cf(81)
      dfluxdub(6,1) = 0.0d0
      dfluxdub(7,1) = 0.0d0
      dfluxdnormal(1,1) = cf(82)
      dfluxdnormal(2,1) = cf(83)
      dfluxdnormal(3,1) = cf(84)
  
      cf(143) = (cf(53)+cf(105))
      cf(144) = (cf(106)+cf(54))
      cf(145) = (cf(107)+cf(55))
c      dflux(2) = u*cf(49)*dvp10
c     &        + cf(51)*dvp40
c     &        + cf(52)*dvp50
c     &        + cf(97)*dvp1
c     &        + cf(98)*dvp4
c     &        + cf(99)*dvp5
c     &        + cf(143)*dnormal(1)
c     &        + cf(144)*dnormal(2)
c     &        + cf(145)*dnormal(3) 
c     &        + cf(100)*dub(1)
c     &        + cf(101)*dub(2)
c     &        + cf(102)*dub(3)
c     &        + cf(103)*dub(4)
c     &        + cf(104)*dub(5)
      dfluxdvp10(2) = cf(49)*u
      dfluxdvp40(2) = cf(51)
      dfluxdvp50(2) = cf(52)
      dfluxdvp1(2) = cf(97)
      dfluxdvp4(2) = cf(98)
      dfluxdvp5(2) = cf(99)
      dfluxdnormal(1,2) = cf(143)
      dfluxdnormal(2,2) = cf(144)
      dfluxdnormal(3,2) = cf(145)
      dfluxdub(1,2) = cf(100)
      dfluxdub(2,2) = cf(101)
      dfluxdub(3,2) = cf(102)
      dfluxdub(4,2) = cf(103) 
      dfluxdub(5,2) = cf(104)
      dfluxdub(6,2) = 0.0d0
      dfluxdub(7,2) = 0.0d0
      cf(146) = (cf(59)+cf(116))
      cf(147) = (cf(117)+cf(60))
      cf(148) = (cf(118)+cf(61))
c      dflux(3) = cf(56)*dvp10
c     &         + cf(57)*dvp40
c     &         + cf(58)*dvp50
c     &         + cf(108)*dvp1
c     &         + cf(109)*dvp4
c     &         + cf(110)*dvp5
c     &         + cf(146)*dnormal(1)
c     &         + cf(147)*dnormal(2)
c     &         + cf(148)*dnormal(3) 
c     &         + cf(111)*dub(1)
c     &         + cf(112)*dub(2)
c     &         + cf(113)*dub(3)
c     &         + cf(114)*dub(4)
c     &         + cf(115)*dub(5)
      dfluxdvp10(3) = cf(56)
      dfluxdvp40(3) = cf(57)
      dfluxdvp50(3) = cf(58)
      dfluxdvp1(3) = cf(108)
      dfluxdvp4(3) = cf(109)
      dfluxdvp5(3) = cf(110)
      dfluxdnormal(1,3) = cf(146)
      dfluxdnormal(2,3) = cf(147)
      dfluxdnormal(3,3) = cf(148)
      dfluxdub(1,3) = cf(111)
      dfluxdub(2,3) = cf(112)
      dfluxdub(3,3) = cf(113)
      dfluxdub(4,3) = cf(114)
      dfluxdub(5,3) = cf(115)
      dfluxdub(6,3) = 0.0d0
      dfluxdub(7,3) = 0.0d0
      cf(149) = (cf(65)+cf(127))
      cf(150) = (cf(128)+cf(66))
      cf(151) = (cf(129)+cf(67))
c      dflux(4) = cf(62)*dvp10
c     &        + cf(63)*dvp40
c     &        + cf(64)*dvp50
c     &        + cf(119)*dvp1
c     &        + cf(120)*dvp4
c     &        + cf(121)*dvp5
c     &        + cf(149)*dnormal(1)
c     &        + cf(150)*dnormal(2)
c     &        + cf(151)*dnormal(3) 
c     &        + cf(122)*dub(1)
c     &        + cf(123)*dub(2)
c     &        + cf(124)*dub(3)
c     &        + cf(125)*dub(4)
c     &        + cf(126)*dub(5)
      dfluxdvp10(4) = cf(62)
      dfluxdvp40(4) = cf(63)
      dfluxdvp50(4) = cf(64)
      dfluxdvp1(4) = cf(119)
      dfluxdvp4(4) = cf(120)
      dfluxdvp5(4) = cf(121)
      dfluxdnormal(1,4) = cf(149)
      dfluxdnormal(2,4) = cf(150)
      dfluxdnormal(3,4) = cf(151)
      dfluxdub(1,4) = cf(122)
      dfluxdub(2,4) = cf(123)
      dfluxdub(3,4) = cf(124)
      dfluxdub(4,4) = cf(125) 
      dfluxdub(5,4) = cf(126)
      dfluxdub(6,4) = 0.0d0
      dfluxdub(7,4) = 0.0d0
      cf(152) = (cf(71)+cf(138))
      cf(153) = (cf(139)+cf(72))
      cf(154) = (cf(140)+cf(73))
c      dflux(5) = cf(68)*dvp10
c     &        + cf(69)*dvp40
c     &        + cf(70)*dvp50
c     &        + cf(152)*dnormal(1)
c     &        + cf(153)*dnormal(2)
c     &        + cf(154)*dnormal(3) 
c     &        + cf(130)*dvp1
c     &        + cf(131)*dvp4
c     &        + cf(132)*dvp5
c     &        + cf(133)*dub(1)
c     &        + cf(134)*dub(2)
c     &        + cf(135)*dub(3)
c     &        + cf(136)*dub(4)
c     &        + cf(137)*dub(5)
      dfluxdvp10(5) = cf(68)
      dfluxdvp40(5) = cf(69)
      dfluxdvp50(5) = cf(70)
      dfluxdnormal(1,5) = cf(152)
      dfluxdnormal(2,5) = cf(153)
      dfluxdnormal(3,5) = cf(154) 
      dfluxdvp1(5) = cf(130)
      dfluxdvp4(5) = cf(131)
      dfluxdvp5(5) = cf(132)
      dfluxdub(1,5) = cf(133)
      dfluxdub(2,5) = cf(134)
      dfluxdub(3,5) = cf(135)
      dfluxdub(4,5) = cf(136)
      dfluxdub(5,5) = cf(137)
      dfluxdub(6,5) = 0.0d0;
      dfluxdub(7,5) = 0.0d0;
      DO 100, I = 1,7
        DO 200, J = 1,3
          dfluxdnormal(J,I) = dfluxdnormal(J,I)
     &                      + dfluxdvp10(I)*dvp10dnormal(J)
     &                      + dfluxdvp40(I)*dvp40dnormal(J)
     &                      + dfluxdvp50(I)*dvp50dnormal(J)
     &                      + dfluxdvp1(I)*dvp1dnormal(J)
     &                      + dfluxdvp4(I)*dvp4dnormal(J)
     &                      + dfluxdvp5(I)*dvp5dnormal(J)
  200 CONTINUE 
        dfluxdsig(I) = dfluxdsig(I)
     &               + dfluxdvp10(I)*dvp10dsig
     &               + dfluxdvp40(I)*dvp40dsig
     &               + dfluxdvp50(I)*dvp50dsig
     &               + dfluxdvp1(I)*dvp1dsig
     &               + dfluxdvp4(I)*dvp4dsig
     &               + dfluxdvp5(I)*dvp5dsig
  100 CONTINUE
c      DO 300, I = 1,5
c        dflux(I) = dflux(I) + dfluxdsig(I)*dsig
c        DO 400, J = 1,3
c          dflux(I) = dflux(I)
c     &             + dfluxdnormal(J,I)*dnormal(J)
c  400 CONTINUE
c        DO 500, J = 1,5
c          dflux(I) = dflux(I)
c     &             + dfluxdub(J,I)*dub(J)
c  500 CONTINUE
c  300 CONTINUE
c
      if (type.eq.1) then
        updir = 0.5d0 + dsign(0.5d0, phiro)
        usro = 1.0d0 / ub(1)
        cf(158) = -1.0d0 / (ub(1)*ub(1))
        cf(155) = (updir * ua(6) + (1.0d0-updir)*ub(6)*usro)
        cf(156) = phiro*(1.0d0-updir)*usro
        cf(157) = phiro*(1.0d0-updir)*ub(6)       
        cf(159) = (cf(155)*cf(142)+cf(157)*cf(158)) 
c        dflux(6) = cf(155)*cf(49)*dvp10
c     &        + cf(155)*cf(50)*dvp40
c     &        + cf(155)*cf(50)*dvp50
c     &        + cf(155)*cf(141)*dvp1
c     &        + cf(155)*cf(75)*dvp4
c     &        + cf(155)*cf(76)*dvp5
c     &        + cf(159)*dub(1)
c     &        + cf(155)*cf(78)*dub(2)
c     &        + cf(155)*cf(79)*dub(3)
c     &        + cf(155)*cf(80)*dub(4)
c     &        + cf(155)*cf(81)*dub(5)
c     &        + cf(156)*dub(6) 
c     &        + cf(155)*cf(82)*dnormal(1)
c     &        + cf(155)*cf(83)*dnormal(2)
c     &        + cf(155)*cf(84)*dnormal(3) 
        dfluxdvp10(6) = cf(155)*cf(49)
        dfluxdvp40(6) = cf(155)*cf(50)
        dfluxdvp50(6) = cf(155)*cf(50)
        dfluxdvp1(6) = cf(155)*cf(141)
        dfluxdvp4(6) = cf(155)*cf(75)
        dfluxdvp5(6) = cf(155)*cf(76)
        dfluxdub(1,6) = cf(159)
        dfluxdub(2,6) = cf(155)*cf(78)
        dfluxdub(3,6) = cf(155)*cf(79)
        dfluxdub(4,6) = cf(155)*cf(80)
        dfluxdub(5,6) = cf(155)*cf(81)
        dfluxdub(6,6) = cf(156)
        dfluxdub(7,6) = 0.0d0 
        dfluxdnormal(1,6) = cf(155)*cf(82)
        dfluxdnormal(2,6) = cf(155)*cf(83)
        dfluxdnormal(3,6) = cf(155)*cf(84)
      else if (type.eq.2) then
        updir = 0.5d0 + dsign(0.5d0, phiro)
        usro = 1.0d0 / ub(1)
        cf(164) = -1.0d0 / (ub(1)*ub(1))
        cf(160) = (updir * ua(6) + (1.0d0-updir)*ub(6)*usro)
        cf(161) = phiro*(1.0d0-updir)*usro
        cf(162) = phiro*(1.0d0-updir)*ub(6)
c        dflux(6) = cf(160)*cf(49)*dvp10
c     &        + cf(160)*cf(50)*dvp40
c     &        + cf(160)*cf(50)*dvp50
c     &        + cf(160)*cf(141)*dvp1
c     &        + cf(160)*cf(75)*dvp4
c     &        + cf(160)*cf(76)*dvp5
c     &        + cf(160)*cf(142)*dub(1)
c     &        + cf(160)*cf(78)*dub(2)
c     &        + cf(160)*cf(79)*dub(3)
c     &        + cf(160)*cf(80)*dub(4)
c     &        + cf(160)*cf(81)*dub(5)
c     &        + cf(161)*dub(6) 
c     &        + cf(160)*cf(82)*dnormal(1)
c     &        + cf(160)*cf(83)*dnormal(2)
c     &        + cf(160)*cf(84)*dnormal(3)
c     &        + cf(162)*cf(164)*dub(1)
        dfluxdvp10(6) = cf(160)*cf(49)
        dfluxdvp40(6) = cf(160)*cf(50)
        dfluxdvp50(6) = cf(160)*cf(50)
        dfluxdvp1(6) = cf(160)*cf(141)
        dfluxdvp4(6) = cf(160)*cf(75)
        dfluxdvp5(6) = cf(160)*cf(76)
        dfluxdub(1,6) = cf(160)*cf(142)+cf(162)*cf(164)
        dfluxdub(2,6) = cf(160)*cf(78)
        dfluxdub(3,6) = cf(160)*cf(79)
        dfluxdub(4,6) = cf(160)*cf(80)
        dfluxdub(5,6) = cf(160)*cf(81) 
        dfluxdub(6,6) = cf(161)
        dfluxdub(7,6) = 0.0d0
        dfluxdnormal(1,6) = cf(160)*cf(82)
        dfluxdnormal(2,6) = cf(160)*cf(83)
        dfluxdnormal(3,6) = cf(160)*cf(84)
        cf(163) = (updir * ua(7) + (1.0d0-updir)*ub(7)*usro)
c        dflux(7) = cf(163)*cf(49)*dvp10
c     &        + cf(163)*cf(50)*dvp40
c     &        + cf(163)*cf(50)*dvp50
c     &        + cf(163)*cf(141)*dvp1
c     &        + cf(163)*cf(75)*dvp4
c     &        + cf(163)*cf(76)*dvp5
c     &        + cf(163)*cf(142)*dub(1)
c     &        + cf(163)*cf(78)*dub(2)
c     &        + cf(163)*cf(79)*dub(3)
c     &        + cf(163)*cf(80)*dub(4)
c     &        + cf(163)*cf(81)*dub(5)
c     &        + phiro*(1.0d0-updir)*usro*dub(7)
c     &        + cf(163)*cf(82)*dnormal(1)
c     &        + cf(163)*cf(83)*dnormal(2)
c     &        + cf(163)*cf(84)*dnormal(3)
c     &        + phiro*(1.0d0-updir)*ub(7)*cf(164)*dub(1)
        dfluxdvp10(7) = cf(163)*cf(49)
        dfluxdvp40(7) = cf(163)*cf(50)
        dfluxdvp50(7) = cf(163)*cf(50)
        dfluxdvp1(7) = cf(163)*cf(141)
        dfluxdvp4(7) = cf(163)*cf(75)
        dfluxdvp5(7) = cf(163)*cf(76)
        dfluxdub(1,7) = cf(163)*cf(142)
     &         + phiro*(1.0d0-updir)*ub(7)*cf(164)
        dfluxdub(2,7) = cf(163)*cf(78)
        dfluxdub(3,7) = cf(163)*cf(79)
        dfluxdub(4,7) = cf(163)*cf(80)
        dfluxdub(5,7) = cf(163)*cf(81)
        dfluxdub(6,7) = 0.0d0
        dfluxdub(7,7) = phiro*(1.0d0-updir)*usro
        dfluxdnormal(1,7) = cf(163)*cf(82)
        dfluxdnormal(2,7) = cf(163)*cf(83)
        dfluxdnormal(3,7) = cf(163)*cf(84)
        dfluxdub(1,7) = phiro*(1.0d0-updir)*ub(7)*cf(164)
      endif
c
c
      END
