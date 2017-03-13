      SUBROUTINE GXROEFLUX6TEMP(type,gamma,gam,pstiff,dpstiff,
     &                      enormal,denormal,evitno,devitno,
     &                      Ugr,dUgr,Ug,dUg,Udr,dUdr,Ud,dUd,phi,dphi)
c-----------------------------------------------------------------------
c This routine computes the Roe Flux derivative taken at the vectors 
c Ug, dUg, Ud and dUd. normal is the normal of the boundary concerned 
c by the flux. phi stores the resulting flux.
c-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 Ug(*), Ud(*), normal(3), enormal(3), evitno, phi(*)
      REAL*8 dUg(*), dUd(*), denormal(3), devitno, dphi(*)
      REAL*8 Ugr(*), Udr(*), energ, enerd
      REAL*8 dUgr(*), dUdr(*)
      REAL*8 H, vitno, updir, gamma 
      REAL*8 dH
      REAL*8 VdotN, rnorm, invnorm
      REAL*8 flur1, flur2, flur3, flur4, flur5
      REAL*8 dif1, dif2, dif3, dif4, dif5
      REAL*8 cr, cr2, qir
      REAL*8 vp1, vp4, vp5
      REAL*8 dvp1, dvp4, dvp5
      REAL*8 ener1, ener2
      REAL*8 uar1, uar2, uar3, uar4, uar5
      REAL*8 duar1, duar2, duar3, duar4, duar5
      REAL*8 usro, squsr1, squsr2
      REAL*8 tet1, tet2, tet3
      REAL*8 gam, gam1, vitg2, vitd2, pstiff
      REAL*8 dpstiff
      REAL*8, DIMENSION(1000) :: cf
      INTEGER type
c
c Initialisation
c
      gam1 = gam - 1.0d0      
c
      cf(1) = 1.0d0 / SQRT(enormal(1)*enormal(1) +
     &         enormal(2)*enormal(2) + enormal(3)*enormal(3))
      rnorm = SQRT(enormal(1)*enormal(1) + enormal(2)*enormal(2) + 
     &        enormal(3)*enormal(3))
      cf(2) = cf(1)*enormal(1)
      cf(3) = cf(1)*enormal(2)
      cf(4) = cf(1)*enormal(3)
      invnorm = 1.0d0 / rnorm
      cf(5) = -1.0d0 / (rnorm*rnorm) *cf(2)
      cf(6) = -1.0d0 / (rnorm*rnorm) *cf(3)
      cf(7) = -1.0d0 / (rnorm*rnorm) *cf(4)
c
      normal(1) = enormal(1) * invnorm
      cf(8) = invnorm+enormal(1)*cf(5)
      cf(9) = enormal(1)*cf(6)
      cf(10) = enormal(1)*cf(7)
      normal(2) = enormal(2) * invnorm
      cf(11) = enormal(2)*cf(5)
      cf(12) = (invnorm+enormal(2)*cf(6))
      cf(13) = enormal(2)*cf(7) 
      normal(3) = enormal(3) * invnorm
      cf(14) = enormal(3)*cf(5)
      cf(15) = enormal(3)*cf(6)
      cf(16) = (enormal(3)*cf(7)+invnorm)
      vitno = evitno * invnorm
      cf(17) = evitno*cf(5)
      cf(18) = evitno*cf(6)
      cf(19) = evitno*cf(7)
c
c Computation of the centred terms
c
      VdotN = Ug(2)*normal(1) + Ug(3)*normal(2) + Ug(4)*normal(3)
      cf(20) = (Ug(2)*cf(8)+Ug(3)*cf(11)+Ug(4)*cf(14))
      cf(21) = (Ug(2)*cf(9)+Ug(3)*cf(12)+Ug(4)*cf(15))
      cf(22) = (Ug(2)*cf(10)+Ug(3)*cf(13)+Ug(4)*cf(16))
      cf(23) = (Ud(2)*cf(8)+Ud(3)*cf(11)+Ud(4)*cf(14))
      cf(24) = (Ud(2)*cf(9)+Ud(4)*cf(15)+Ud(3)*cf(12))
      cf(25) = (Ud(2)*cf(10)+Ud(3)*cf(13)+Ud(4)*cf(16))
      cf(26) = 2.0d0*Ug(2)
      cf(27) = 2.0d0*Ug(3)
      cf(28) = 2.0d0*Ug(4)
      vitg2 = Ug(2)*Ug(2) + Ug(3)*Ug(3) + Ug(4)*Ug(4)
      H = gam*(Ug(5)+pstiff) + 0.5*gam1*Ug(1)*vitg2
      cf(29) = 0.5*gam1*vitg2
      cf(30) = 0.5*gam1*Ug(1)*cf(26)
      cf(31) = 0.5*gam1*Ug(1)*cf(27)
      cf(32) = 0.5*gam1*Ug(1)*cf(28)
      cf(33) = gam1*Ug(1)/((gam1*Ug(1))*(gam1*Ug(1)))
      cf(34) = -H*gam1/((gam1*Ug(1))*(gam1*Ug(1)))
      cf(35) = (cf(33)*cf(29)+cf(34))
      cf(36) = cf(33)*cf(30)
      cf(37) = cf(33)*cf(31)
      cf(38) = cf(33)*cf(32)
      cf(39) = cf(33)*gam
      cf(40) = Ug(1)*(cf(20)-cf(17))
      cf(41) = Ug(1)*(cf(21)-cf(18))
      cf(42) = Ug(1)*(cf(22)-cf(19))
      cf(43) = Ug(2)*normal(1) + Ug(3)*normal(2) 
     %       + Ug(4)*normal(3) - vitno
      cf(44) = Ug(1)*normal(1)
      cf(45) = Ug(1)*normal(2)
      cf(46) = Ug(1)*normal(3)
      cf(47) = -Ug(1)*invnorm
      H = H/(gam1*Ug(1))
      phi(1) = Ug(1)*(Ug(2)*normal(1) + Ug(3)*normal(2) +
     &         Ug(4)*normal(3) - vitno)
      phi(2) = phi(1)*Ug(2) + Ug(5)*normal(1)
      cf(48) = (Ug(2)*cf(40)+Ug(5)*cf(8))
      cf(49) = (Ug(2)*cf(41)+Ug(5)*cf(9))
      cf(50) = (Ug(2)*cf(42)+Ug(5)*cf(10))
      cf(51) = Ug(2)*cf(43)
      cf(52) = Ug(2)*cf(44)
      cf(53) = Ug(2)*cf(45)
      cf(54) = Ug(2)*cf(46)
      cf(55) = Ug(2)*cf(47)
      cf(56) = phi(1)
      cf(57) = (Ug(4)*cf(40)+Ug(5)*cf(14))
      cf(58) = (Ug(4)*cf(41)+Ug(5)*cf(15))
      cf(59) = (Ug(4)*cf(42)+Ug(5)*cf(16))
      cf(60) = Ug(4)*cf(43)
      cf(61) = Ug(4)*cf(44)
      cf(62) = Ug(4)*cf(45)
      cf(63) = (Ug(4)*cf(46)+cf(56))
      cf(64) = Ug(4)*cf(47) 
      phi(3) = phi(1)*Ug(3) + Ug(5)*normal(2)
      phi(4) = phi(1)*Ug(4) + Ug(5)*normal(3)
      phi(5) = phi(1)*H + Ug(5)*vitno
      cf(65) = (H*cf(40) + Ug(5)*cf(17))
      cf(66) = (H*cf(41) + Ug(5)*cf(18))
      cf(67) = (H*cf(42) + Ug(5)*cf(19))
      cf(68) = (H*cf(47) + Ug(5)*invnorm)
      cf(69) = (H*cf(43) + cf(56)*cf(35))
      cf(70) = (H*cf(44) + cf(56)*cf(36))
      cf(71) = (H*cf(45) + cf(56)*cf(37))
      cf(72) = (H*cf(46) + cf(56)*cf(38))
      cf(73) = cf(56)*cf(39)
c
      VdotN = Ud(2)*normal(1)+Ud(3)*normal(2)+Ud(4)*normal(3)
      VdotN = VdotN - vitno
      cf(74) = (cf(23)-cf(17))
      cf(75) = (cf(24)-cf(18))
      cf(76) = (cf(25)-cf(19))
      cf(77) = 2.0d0*Ud(2)
      cf(78) = 2.0d0*Ud(3)
      cf(79) = 2.0d0*Ud(4)
      vitd2 = Ud(2)*Ud(2) + Ud(3)*Ud(3) + Ud(4)*Ud(4)
      H = gam*(Ud(5)+pstiff)+0.5*gam1*Ud(1)*vitd2
      cf(80) = 0.5*gam1*vitd2
      cf(81) = 0.5*gam1*Ud(1)*cf(77)
      cf(82) = 0.5*gam1*Ud(1)*cf(78)
      cf(83) = 0.5*gam1*Ud(1)*cf(79)
      cf(84) = gam1*Ud(1)/((gam1*Ud(1))*(gam1*Ud(1)))
      cf(85) = -H*gam1/((gam1*Ud(1))*(gam1*Ud(1)))
      H = H/(gam1*Ud(1))
      phi(1) = phi(1) + Ud(1)*VdotN
      cf(86) = (cf(40)+Ud(1)*cf(74))
      cf(87) = (cf(41)+Ud(1)*cf(75))
      cf(88) = (cf(42)+Ud(1)*cf(76))
      cf(89) = (cf(47)-Ud(1)*invnorm)
      cf(90) = Ud(1)*normal(1)
      cf(91) = Ud(1)*normal(2)
      cf(92) = Ud(1)*normal(3)
      cf(93) = VdotN
      phi(2) = phi(2) + Ud(1)*Ud(2)*cf(93) +Ud(5)*normal(1)
      cf(94) = (cf(48)+Ud(1)*Ud(2)*cf(74)+Ud(5)*cf(8))
      cf(95) = (cf(49)+Ud(1)*Ud(2)*cf(75)+Ud(5)*cf(9))
      cf(96) = (cf(50)+Ud(1)*Ud(2)*cf(76)+Ud(5)*cf(10)) 
      cf(97) = (cf(55)-Ud(1)*Ud(2)*invnorm)
      cf(98) = Ud(2)*cf(93)
      cf(99) = cf(93)*Ud(1)
      cf(100) = Ud(1)*Ud(2)*normal(1)
      cf(101) = Ud(1)*Ud(2)*normal(2)
      cf(102) = Ud(1)*Ud(2)*normal(3)
      phi(3) = phi(3) + Ud(1)*Ud(3)*cf(93) +Ud(5)*normal(2)
      cf(103) = ((Ug(3)*cf(40) + Ug(5)*cf(11)) + Ud(1)*Ud(3)*cf(74) +
     &          Ud(5)*cf(11))
      cf(104) = ((Ug(3)*cf(41) + Ug(5)*cf(12)) + Ud(1)*Ud(3)*cf(75) +
     &          Ud(5)*cf(12))
      cf(105) = ((Ug(3)*cf(42) + Ug(5)*cf(13)) + Ud(1)*Ud(3)*cf(76) +
     &          Ud(5)*cf(13))
      cf(106) = (Ug(3)*cf(47) - Ud(1)*Ud(3)*invnorm)
      cf(107) = Ug(3)*cf(43)
      cf(108) = Ug(3)*cf(44)
      cf(109) = (Ug(3)*cf(45) + cf(56))
      cf(110) = Ug(3)*cf(46)
      cf(111) = Ud(3)*cf(93)
      cf(112) = Ud(1)*Ud(3)*normal(1)
      cf(113) = (cf(93)*Ud(1) + Ud(1)*Ud(3)*normal(2))
      cf(114) = Ud(1)*Ud(3)*normal(3)
      phi(4) = phi(4) + Ud(1)*Ud(4)*cf(93) +Ud(5)*normal(3)
      cf(115) = (cf(57)+Ud(5)*cf(14)+Ud(1)*Ud(4)*cf(74))
      cf(116) = (cf(58)+Ud(1)*Ud(4)*cf(75)+Ud(5)*cf(15))
      cf(117) = (cf(59)+Ud(1)*Ud(4)*cf(76)+Ud(5)*cf(16))
      cf(118) = (cf(64)-Ud(1)*Ud(4)*invnorm)
      cf(119) = Ud(4)*cf(93)
      cf(120) = Ud(1)*Ud(4)*normal(1)
      cf(121) = Ud(1)*Ud(4)*normal(2)
      cf(122) = (Ud(1)*cf(93)+Ud(1)*Ud(4)*normal(3))
      cf(123) = (cf(93)*H+Ud(1)*cf(93)*cf(84)*cf(80)
     &        + Ud(1)*cf(93)*cf(85))
      cf(124) = (Ud(1)*H*normal(1)+Ud(1)*cf(93)*cf(84)*cf(81))
      cf(125) = (Ud(1)*H*normal(2)+Ud(1)*cf(93)*cf(84)*cf(82))
      cf(126) = (Ud(1)*H*normal(3)+Ud(1)*cf(93)*cf(84)*cf(83))
      phi(5) = phi(5) + Ud(1)*cf(93)*H +Ud(5)*vitno
      cf(127) = (cf(65)+Ud(1)*H*cf(74)+Ud(5)*cf(17))
      cf(128) = (cf(66)+Ud(1)*H*cf(75)+Ud(5)*cf(18))
      cf(129) = (cf(67)+Ud(1)*H*cf(76)+Ud(5)*cf(19))
      cf(130) = (cf(68)+Ud(5)*invnorm-Ud(1)*H*invnorm)
      cf(131) = (cf(73)+Ud(1)*cf(93)*cf(84)*gam)
      cf(132) = (cf(73)+vitno)
      cf(133) = (Ud(1)*cf(93)*cf(84)*gam+vitno)
c
c Computation of the Roe-averaged state
c
      squsr1 = SQRT(Ug(1))
      cf(134) = 1.0d0/(2.0d0*SQRT(Ug(1)))
      squsr2 = SQRT(Ud(1))
      cf(135) = 1.0d0/(2.0d0*SQRT(Ud(1)))
c     
      cf(136) = 1.0d0/gam1
      ener1 = (Ug(5)+gam*pstiff)*cf(136) + 0.5d0*Ug(1)*vitg2
      cf(137) = cf(136)*gam
      cf(138) = 0.5d0*vitg2
      cf(139) = 0.5d0*Ug(1)*cf(26)
      cf(140) = 0.5d0*Ug(1)*cf(27)
      cf(141) = 0.5d0*Ug(1)*cf(28)
c
      ener2 = cf(136)*(Ud(5)+gam*pstiff) + 0.5d0*Ud(1)*vitd2
      cf(142) = cf(136)*gam
      cf(143) = 0.5d0*vitd2
      cf(144) = 0.5d0*Ud(1)*cf(77)
      cf(145) = 0.5d0*Ud(1)*cf(78)
      cf(146) = 0.5d0*Ud(1)*cf(79)
c
      usro  = 1.0d0/(squsr1 + squsr2)
      cf(147) = -1.0d0/((squsr1 + squsr2)*(squsr1 + squsr2)) 
      cf(148) = cf(147)*cf(134)
      cf(149) = cf(147)*cf(135)
c
      uar1 = (squsr1*Ug(1) + squsr2*Ud(1))*usro
c      duar1 = (cf(134)*dUg(1)*Ug(1) 
c     &      + squsr1*dUg(1) + cf(135)*dUd(1)*Ud(1) + 
c     &        squsr2*dUd(1))*usro + (squsr1*Ug(1) + 
c     &        squsr2*Ud(1))*(cf(148)*dUg(1)+cf(149)*dUd(1))
c
      uar2 = (squsr1*Ug(2) + squsr2*Ud(2))*usro
      cf(150) = (squsr1*Ug(2) + squsr2*Ud(2))
      cf(151) = (usro*cf(134)*Ug(2)+cf(150)*cf(148))
      cf(152) = usro*squsr1
      cf(153) = (usro*cf(135)*Ud(2)+cf(150)*cf(149))
      cf(154) = usro*squsr2
      cf(155) = (squsr1*Ug(3) + squsr2*Ud(3))
c
      uar3 = (squsr1*Ug(3) + squsr2*Ud(3))*usro
      cf(155) = (squsr1*Ug(3) + squsr2*Ud(3))
      cf(156) = (usro*cf(134)*Ug(3)+cf(155)*cf(148))
      cf(157) = usro*squsr1
      cf(158) = (cf(155)*cf(149)+usro*cf(135)*Ud(3))
      cf(159) = usro*squsr2
c
      uar4 = (squsr1*Ug(4) + squsr2*Ud(4))*usro
      cf(160) = (squsr1*Ug(4) + squsr2*Ud(4))
      cf(161) = (usro*cf(134)*Ug(4)+cf(160)*cf(148))
      cf(162) = usro*squsr1
      cf(163) = (cf(160)*cf(149)+usro*cf(135)*Ud(4))
      cf(164) = usro*squsr2
      cf(165) = 1.0d0/(squsr1*squsr1)
c
      uar5 = ((ener1 + Ug(5))/squsr1 + (ener2 + Ud(5))/squsr2)*usro
      cf(166) = 1.0d0/(squsr2*squsr2)
      cf(167) = ((ener1 + Ug(5))/squsr1+(ener2 + Ud(5))/squsr2)
      cf(168) = (usro*cf(165)*squsr1*cf(137)
     &        + usro*cf(166)*squsr2*cf(142))
      cf(169) = (usro*cf(165)*squsr1*cf(138)
     &        - usro*cf(165)*(ener1 + Ug(5))*cf(134)+cf(167)*cf(148))
      cf(170) = (usro*cf(165)*squsr1*cf(136)+usro*cf(165)*squsr1)
      cf(171) = usro*cf(165)*squsr1*cf(139)
      cf(172) = usro*cf(165)*squsr1*cf(140)
      cf(173) = usro*cf(165)*squsr1*cf(141)
      cf(174) = (usro*cf(166)*squsr2*cf(143)-usro*cf(166)*(ener2 +
     &          Ud(5))*cf(135)+cf(167)*cf(149))
      cf(175) = (usro*cf(166)*squsr2*cf(136)+usro*cf(166)*squsr2)
      cf(176) = usro*cf(166)*squsr2*cf(144)
      cf(177) = usro*cf(166)*squsr2*cf(145)
      cf(178) = usro*cf(166)*squsr2*cf(146)
c
c Computation of the dissipation term
c
      VdotN = normal(1)*uar2 + normal(2)*uar3 + normal(3)*uar4
      cf(179) = (uar2*cf(8)+uar3*cf(11)+uar4*cf(14))
      cf(180) = (uar2*cf(9)+uar4*cf(15)+uar3*cf(12))
      cf(181) = (uar2*cf(10)+uar3*cf(13)+uar4*cf(16))
      cf(182) = (normal(1)*cf(151)+normal(2)*cf(156)+normal(3)*cf(161))
      cf(183) = (normal(1)*cf(153)+normal(2)*cf(158)+normal(3)*cf(163))
      cf(184) = normal(1)*cf(152)
      cf(185) = normal(2)*cf(157)
      cf(186) = normal(3)*cf(162)
      cf(187) = normal(1)*cf(154)
      cf(188) = normal(2)*cf(159)
      cf(189) = normal(3)*cf(164)
c
      qir = 0.5d0*(uar2*uar2 + uar3*uar3 + uar4*uar4)
      cf(190) = (uar2*cf(151)+uar4*cf(161)+uar3*cf(156))
      cf(191) = (uar2*cf(153)+uar4*cf(163)+uar3*cf(158))
      cf(192) = uar2*cf(152)
      cf(193) = uar3*cf(157)
      cf(194) = uar4*cf(162)
      cf(195) = uar2*cf(154)
      cf(196) = uar3*cf(159)
      cf(197) = uar4*cf(164)
c
      tet1 = normal(3)*uar3 - normal(2)*uar4
      cf(198) = (uar3*cf(14)-uar4*cf(11))
      cf(199) = (uar3*cf(15)-uar4*cf(12))
      cf(200) = (uar3*cf(16)-uar4*cf(13))
      cf(201) = (normal(3)*cf(156)-normal(2)*cf(161))
      cf(202) = normal(3)*cf(157)
      cf(203) = -normal(2)*cf(162)
      cf(204) = (normal(3)*cf(158)-normal(2)*cf(163))
      cf(205) = normal(3)*cf(159)
      cf(206) = -normal(2)*cf(164)
      tet2 = normal(1)*uar4 - normal(3)*uar2
      cf(207) = (uar4*cf(8)-uar2*cf(14))
      cf(208) = (uar4*cf(9)-uar2*cf(15))
      cf(209) = (uar4*cf(10)-uar2*cf(16))
      cf(210) = (normal(1)*cf(161)-normal(3)*cf(151))
      cf(211) = - normal(3)*cf(152)
      cf(212) = normal(1)*cf(162)
      cf(213) = (normal(1)*cf(163)-normal(3)*cf(153))
      cf(214) = - normal(3)*cf(154)
      cf(215) =  normal(1)*cf(164)
      tet3 = normal(2)*uar2 - normal(1)*uar3
      cf(216) = (uar2*cf(11)-uar3*cf(8))
      cf(217) = (uar2*cf(12)-uar3*cf(9))
      cf(218) = (uar2*cf(13)-uar3*cf(10))
      cf(219) = (normal(2)*cf(151)-normal(1)*cf(156))
      cf(220) = normal(2)*cf(152)
      cf(221) = -normal(1)*cf(157)
      cf(222) = (normal(2)*cf(153)-normal(1)*cf(158))
      cf(223) = normal(2)*cf(154)
      cf(224) = -normal(1)*cf(159)
c
      cr2 = gam1*(uar5 - qir)
      cf(225) = gam1*cf(168)
      cf(226) = gam1*(cf(169)-cf(190))
      cf(227) = gam1*(cf(171)-cf(192))
      cf(228) = gam1*(cf(172)-cf(193))
      cf(229) = gam1*(cf(173)-cf(194))
      cf(230) = gam1*cf(170)
      cf(231) = gam1*(cf(174)-cf(191))
      cf(232) = gam1*(cf(176)-cf(195))
      cf(233) = gam1*(cf(177)-cf(196))
      cf(234) = gam1*(cf(178)-cf(197))
      cf(235) = gam1*cf(175)
      cr = SQRT(cr2)
      cf(236) = 1.0d0/(2.0d0*SQRT(cr2))
      cf(237) = cf(236)*cf(225)
      cf(238) = cf(236)*cf(226)
      cf(239) = cf(236)*cf(227)
      cf(240) = cf(236)*cf(228)
      cf(241) = cf(236)*cf(229)
      cf(242) = cf(236)*cf(230)
      cf(243) = cf(236)*cf(231)
      cf(244) = cf(236)*cf(232)
      cf(245) = cf(236)*cf(233)
      cf(246) = cf(236)*cf(234)
      cf(247) = cf(236)*cf(235)
      cf(248) = -1.0d0/(cr2*cr2)
      cf(249) = cf(248)*cf(225)
      cf(250) = cf(248)*cf(226)
      cf(251) = cf(248)*cf(227)
      cf(252) = cf(248)*cf(228)
      cf(253) = cf(248)*cf(229)
      cf(254) = cf(248)*cf(230)
      cf(255) = cf(248)*cf(231)
      cf(256) = cf(248)*cf(232)
      cf(257) = cf(248)*cf(233)
      cf(258) = cf(248)*cf(234)
      cf(259) = cf(248)*cf(235)
      cr2 = 1.0/cr2
c
      energ = (Ugr(5)+gam*pstiff)*cf(136) 
     &       + 0.5d0*Ugr(1)*(Ugr(2)*Ugr(2) + 
     &        Ugr(3)*Ugr(3) + Ugr(4)*Ugr(4))
      cf(260) = 0.5d0*(Ugr(2)*Ugr(2) + Ugr(3)*Ugr(3) + Ugr(4)*Ugr(4))
      cf(261) = cf(136)*gam
      cf(262) = Ugr(1)*Ugr(2)
      cf(263) = Ugr(1)*Ugr(3)
      cf(264) = Ugr(1)*Ugr(4)
      enerd = (Udr(5)+gam*pstiff)*cf(136) 
     &        + 0.5d0*Udr(1)*(Udr(2)*Udr(2) + 
     &        Udr(3)*Udr(3) + Udr(4)*Udr(4))
      cf(265) = 0.5d0*Udr(2)*Udr(2)+
     &          0.5d0*Udr(3)*Udr(3)+0.5d0*Udr(4)*Udr(4)
      cf(266) = Udr(1)*Udr(2)
      cf(267) = Udr(1)*Udr(3)
      cf(268) = Udr(1)*Udr(4)
      cf(269) = cf(136)*gam
      dif1 = - Ugr(1) + Udr(1)
      dif2 = - Ugr(1)*Ugr(2) + Udr(1)*Udr(2)
      dif3 = - Ugr(1)*Ugr(3) + Udr(1)*Udr(3)
      dif4 = - Ugr(1)*Ugr(4) + Udr(1)*Udr(4)
      dif5 = - energ + enerd
c
      vp1 = gamma * ABS(VdotN - vitno)
      if (ABS(VdotN - vitno) .eq. 0.0d0) then
        dvp1 = 0.0d0
      else
        cf(270) = gamma * (VdotN - vitno)/ABS(VdotN - vitno)
        cf(271) = cf(270)*(cf(179)-cf(17))
        cf(272) = cf(270)*(cf(180)-cf(18))
        cf(273) = cf(270)*(cf(181)-cf(19))
        cf(274) = -cf(270)*invnorm
        cf(275) = cf(270)*cf(182)
        cf(276) = cf(270)*cf(184)
        cf(277) = cf(270)*cf(185)
        cf(278) = cf(270)*cf(186)
        cf(279) = cf(270)*cf(183)
        cf(280) = cf(270)*cf(187)
        cf(281) = cf(270)*cf(188)
        cf(282) = cf(270)*cf(189) 
        dvp1 = 
     &         cf(271)*denormal(1)
     &       + cf(272)*denormal(2)
     &       + cf(273)*denormal(3)
     &       - cf(274)*devitno 
     &       + cf(275)*dUg(1)
     &       + cf(276)*dUg(2)
     &       + cf(277)*dUg(3)
     &       + cf(278)*dUg(4)
     &       + cf(279)*dUd(1)
     &       + cf(280)*dUd(2)
     &       + cf(281)*dUd(3)
     &       + cf(282)*dUd(4) 
      end if
      vp4 = gamma * ABS((VdotN + cr) - vitno)
      if (ABS((VdotN + cr) - vitno) .eq. 0.0d0) then
        dvp4 = 0.0d0
      else
        cf(283) = gamma * ((VdotN + cr) 
     &    - vitno)/ABS((VdotN + cr) - vitno)
        cf(284) = cf(283)*(cf(179)-cf(17))
        cf(285) = cf(283)*(cf(180)-cf(18))
        cf(286) = cf(283)*(cf(181)-cf(19))
        cf(287) = cf(283)*cf(237)
        cf(288) = cf(283)*(cf(182)+cf(238))
        cf(289) = cf(283)*(cf(184)+cf(239))
        cf(290) = cf(283)*(cf(185)+cf(240))
        cf(291) = cf(283)*(cf(186)+cf(241))
        cf(292) = cf(283)*cf(242)
        cf(293) = cf(283)*(cf(183)+cf(243))
        cf(294) = cf(283)*(cf(187)+cf(244))
        cf(295) = cf(283)*(cf(188)+cf(245))
        cf(296) = cf(283)*(cf(189)+cf(246))
        cf(297) = cf(283)*cf(247)
        cf(298) = - cf(283)*invnorm 
        dvp4 = 
     &         cf(284)*denormal(1)
     &       + cf(285)*denormal(2)
     &       + cf(286)*denormal(3)
     &       + cf(287)*dpstiff
     &       + cf(288)*dUg(1)
     &       + cf(289)*dUg(2)
     &       + cf(290)*dUg(3)
     &       + cf(291)*dUg(4)
     &       + cf(292)*dUg(5)
     &       + cf(293)*dUd(1)
     &       + cf(294)*dUd(2)
     &       + cf(295)*dUd(3)
     &       + cf(296)*dUd(4)
     &      + cf(297)*dUd(5) 
     &      + cf(298)*devitno
      end if
      vp5 = gamma * ABS((VdotN - cr) - vitno)
      if (ABS((VdotN - cr) - vitno) .eq. 0.0d0) then
        dvp5 = 0.0d0
      else
        cf(299) = gamma * ((VdotN - cr) 
     &            - vitno)/ABS((VdotN - cr) - vitno)
        cf(300) = cf(299)*(cf(179)-cf(17))
        cf(301) = cf(299)*(cf(180)-cf(18))
        cf(302) = cf(299)*(cf(181)-cf(19))
        cf(303) = -cf(299)*cf(237)
        cf(304) = cf(299)*(cf(182)-cf(238))
        cf(305) = cf(299)*(cf(184)-cf(239))
        cf(306) = cf(299)*(cf(185)-cf(240))
        cf(307) = cf(299)*(cf(186)-cf(241))
        cf(308) = -cf(299)*cf(242)
        cf(309) = cf(299)*(cf(183)-cf(243))
        cf(310) = cf(299)*(cf(187)-cf(244))
        cf(311) = cf(299)*(cf(188)-cf(245))
        cf(312) = cf(299)*(cf(189)-cf(246))
        cf(313) = -cf(299)*cf(247)
        cf(314) = -cf(299)*invnorm
        dvp5 = 
     &         cf(300)*denormal(1)
     &       + cf(301)*denormal(2)
     &       + cf(302)*denormal(3)
     &       + cf(303)*dpstiff
     &       + cf(304)*dUg(1)
     &       + cf(305)*dUg(2)
     &       + cf(306)*dUg(3)
     &       + cf(307)*dUg(4)
     &       + cf(308)*dUg(5)
     &       + cf(309)*dUd(1)
     &       + cf(310)*dUd(2)
     &       + cf(311)*dUd(3)
     &       + cf(312)*dUd(4)
     &       + cf(313)*dUd(5) 
     &       + cf(314)*devitno
      end if
c
      flur1 = vp1*
     &        ((normal(1)*(1.0d0 - gam1*qir*cr2) - tet1)*dif1 +
     &        (normal(1)*gam1*uar2*cr2)*dif2 +
     &        (normal(3) + (normal(1)*gam1*uar3*cr2))*dif3 +
     &        (-normal(2) + (normal(1)*gam1*uar4*cr2))*dif4 -
     &        (normal(1)*gam1*cr2)*dif5)
      cf(315) = ((normal(1)*(1.0d0 - gam1*qir*cr2) - tet1)*dif1 +
     &         (normal(1)*gam1*uar2*cr2)*dif2 +
     &         (normal(3) + (normal(1)*gam1*uar3*cr2))*dif3 +
     &         (-normal(2) + (normal(1)*gam1*uar4*cr2))*dif4 -
     &         (normal(1)*gam1*cr2)*dif5)
      cf(316) = (1.0d0 - gam1*qir*cr2)
      cf(317) = (normal(1)*(1.0d0-gam1*qir*cr2)-tet1)
      cf(318) = gam1*uar2*cr2
      cf(319) = dif1*(cf(316)*cf(8)-cf(198))
      cf(320) = dif1*(cf(316)*cf(9)-cf(199))
      cf(321) = dif1*(cf(316)*cf(10)-cf(200))
      cf(322) = -dif1*normal(1)*gam1*qir*cf(249)
      cf(323) = - dif1*(normal(1)*(gam1*cr2*cf(190)
     &  +gam1*qir*cf(250))+cf(201))
      cf(324) = - dif1*normal(1)*(gam1*cr2*cf(192)+gam1*qir*cf(251))
      cf(325) = - dif1*(normal(1)*(gam1*cr2*cf(193)
     & +gam1*qir*cf(252))+cf(202))
      cf(326) = - dif1*(normal(1)*(gam1*cr2*cf(194)
     &  +gam1*qir*cf(253))+cf(203))
      cf(327) = - dif1*normal(1)*gam1*qir*cf(254)
      cf(328) = - dif1*(normal(1)*(gam1*cr2*cf(191)
     &  +gam1*qir*cf(255))+cf(204))
      cf(329) = - dif1*normal(1)*(gam1*cr2*cf(195)
     &  +gam1*qir*cf(256))
      cf(330) = - dif1*(normal(1)*(gam1*cr2*cf(196)
     &  +gam1*qir*cf(257))+cf(205))
      cf(331) = - dif1*(normal(1)*(gam1*cr2*cf(197)
     & +gam1*qir*cf(258))+cf(206))
      cf(332) = - dif1*normal(1)*gam1*qir*cf(259)
      cf(333) = normal(1)*gam1*cr2
      cf(334) = normal(1)*gam1*uar2
      cf(335) = (normal(1)*gam1*uar2*cr2)
      cf(336) = gam1*uar3*cr2
      cf(337) = normal(1)*gam1*uar3
      cf(338) = (cf(319)+dif2*cf(318)*cf(8)+dif3*(cf(14)+cf(336)*cf(8)))
      cf(339) = (cf(320)+dif3*(cf(15)+cf(336)*cf(9))+dif2*cf(318)*cf(9))
      cf(340)=(cf(321)+dif2*cf(318)*cf(10)+dif3*(cf(16)+cf(336)*cf(10)))
      cf(341) = ((dif2*(cf(334)*cf(250)
     &      +cf(333)*cf(151))+cf(323))-cf(317)
     &      +dif3*(cf(337)*cf(250)+cf(333)*cf(156)))
      cf(342) = (dif3*cf(337)*cf(251)
     &      +(dif2*(cf(334)*cf(251)+cf(333)*cf(152))+cf(324)))
      cf(343) = (dif3*(cf(333)*cf(157)
     &      +cf(337)*cf(252))+(cf(325)+dif2*cf(334)*cf(252)))
      cf(344) = (dif3*cf(337)*cf(253)+(cf(326)+dif2*cf(334)*cf(253)))
      cf(345) = (dif3*cf(337)*cf(254)+(cf(327)+dif2*cf(334)*cf(254))) 
      cf(346) = ((cf(328)+dif2*(cf(334)*cf(255)+cf(333)*cf(153)))
     &       +dif3*(cf(333)*cf(158)+cf(337)*cf(255)))
      cf(347) = (dif3*cf(337)*cf(256)+(cf(329)
     &          +dif2*(cf(333)*cf(154)+cf(334)*cf(256))))
      cf(348) = (dif3*(cf(333)*cf(159)+cf(337)*cf(257))+(cf(330)+
     &      dif2*cf(334)*cf(257)))
      cf(349) = (dif3*cf(337)*cf(258)+(cf(331)+dif2*cf(334)*cf(258)))
      cf(350) = (dif3*cf(337)*cf(259)+(cf(332)+dif2*cf(334)*cf(259)))
      cf(351) = (cf(322)+dif3*cf(337)*cf(249)+dif2*cf(334)*cf(249))
      cf(352) = (normal(3) + (normal(1)*gam1*uar3*cr2))
      cf(353) = gam1*uar4*cr2
      cf(354) = normal(1)*gam1*cr2
      cf(355) = normal(1)*gam1*uar4
      cf(356) = (cf(353)*cf(8)-cf(11))
      cf(357) = (cf(353)*cf(9)-cf(12))
      cf(358) = (cf(353)*cf(10)-cf(13))
      cf(359) = (-normal(2) + (normal(1)*gam1*uar4*cr2))
      cf(360) = normal(1)*gam1
      cf(361) = (cf(338)+ dif4*cf(356)-dif5*gam1*cr2*cf(8))
      cf(362) = (cf(339)+dif4*cf(357)-dif5*gam1*cr2*cf(9))
      cf(363) = (cf(340)+dif4*cf(358)-dif5*gam1*cr2*cf(10))
      cf(364) = (cf(351)+dif4*cf(355)*cf(249)-dif5*cf(360)*cf(249))
      cf(365) = ((cf(335)*Udr(2)+cf(317))+cf(359)*Udr(4)+cf(352)*Udr(3))
      cf(366) = (cf(341)+dif4*(cf(354)*cf(161)
     &        + cf(355)*cf(250))-dif5*cf(360)*cf(250))
      cf(367) = (cf(342)+dif4*cf(355)*cf(251)-dif5*cf(360)*cf(251)) 
      cf(368) = (cf(343)+dif4*cf(355)*cf(252)-dif5*cf(360)*cf(252)) 
      cf(369) = (cf(344)+dif4*(cf(354)*cf(162)
     &        + cf(355)*cf(253))-dif5*cf(360)*cf(253))
      cf(370) = (cf(345)+dif4*cf(355)*cf(254)-dif5*cf(360)*cf(254))
      cf(371) = (cf(346)+dif4*(cf(354)*cf(163)
     &         +cf(355)*cf(255))-dif5*cf(360)*cf(255))
      cf(372) = (cf(347)+dif4*cf(355)*cf(256)-dif5*cf(360)*cf(256))
      cf(373) = (cf(348)+dif4*cf(355)*cf(257)-dif5*cf(360)*cf(257))
      cf(374) = (cf(349)+dif4*(cf(354)*cf(164)
     &         +cf(355)*cf(258))-dif5*cf(360)*cf(258))
      cf(375) = (cf(350)+dif4*cf(355)*cf(259)-dif5*cf(360)*cf(259))
      cf(376) = - (normal(1)*gam1*cr2)
      cf(377) = (cf(364)+cf(376)*(cf(269)-cf(261))) 
      cf(378) = -((cf(352)*Ugr(3)+cf(335)*Ugr(2)
     &         +cf(359)*Ugr(4))+cf(376)*cf(260))
      cf(379) = vp1*cf(361)
      cf(380) = vp1*cf(362)
      cf(381) = vp1*cf(363)
      cf(382) = vp1*cf(377)
      cf(383) = vp1*cf(366)
      cf(384) = vp1*cf(367)
      cf(385) = vp1*cf(368)
      cf(386) = vp1*cf(369)
      cf(387) = vp1*cf(370)
      cf(388) = vp1*cf(371)
      cf(389) = vp1*cf(372)
      cf(390) = vp1*cf(373)
      cf(391) = vp1*cf(374)
      cf(392) = vp1*cf(375)
      cf(393) = vp1*(cf(365)+cf(376)*cf(265))
      cf(394) = vp1*(cf(335)*Udr(1)+cf(376)*cf(266))
      cf(395) = vp1*(cf(352)*Udr(1)+cf(376)*cf(267)) 
      cf(396) = vp1*(cf(359)*Udr(1)+cf(376)*cf(268))
      cf(397) = vp1*cf(376)*cf(136)
      cf(398) = vp1*cf(378) 
      cf(399) = -vp1*(cf(335)*Ugr(1)+cf(376)*cf(262))
      cf(400) = -vp1*(cf(352)*Ugr(1)+cf(376)*cf(263))
      cf(401) = - vp1*(cf(359)*Ugr(1)+cf(376)*cf(264))
      cf(402) = - vp1*cf(376)*cf(136)
c     &         + cf(400)*dUgr(3)
c     &         + cf(401)*dUgr(4)
c     &         + cf(402)*dUgr(5)
c
      flur2 = vp1*
     &        ((normal(2)*(1.0d0 - gam1*qir*cr2) - tet2)*dif1 +
     &        (-normal(3) + (normal(2)*gam1*uar2*cr2))*dif2 +
     &        (normal(2)*gam1*uar3*cr2)*dif3 +
     &        (normal(1)  + (normal(2)*gam1*uar4*cr2))*dif4 -
     &        (normal(2)*gam1*cr2)*dif5)
      cf(403) = ((normal(2)*(1.0d0 - gam1*qir*cr2) - tet2)*dif1 +
     &         (-normal(3) + (normal(2)*gam1*uar2*cr2))*dif2 +
     &         (normal(2)*gam1*uar3*cr2)*dif3 +
     &         (normal(1)  + (normal(2)*gam1*uar4*cr2))*dif4 -
     &         (normal(2)*gam1*cr2)*dif5)
      cf(404) = (1.0d0 - gam1*qir*cr2)
      cf(405) = gam1*cr2
      cf(406) = gam1*qir
      cf(407) = - normal(2)*(cf(405)*cf(190)+cf(406)*cf(250))
      cf(408) = - normal(2)*(cf(405)*cf(192)+cf(406)*cf(251))
      cf(409) = - normal(2)*(cf(405)*cf(193)+cf(406)*cf(252))
      cf(410) = - normal(2)*(cf(405)*cf(194)+cf(406)*cf(253))
      cf(411) = - normal(2)*cf(406)*cf(254)
      cf(412) = - normal(2)*(cf(405)*cf(191)+cf(406)*cf(255))
      cf(413) = - normal(2)*(cf(405)*cf(195)+cf(406)*cf(256))
      cf(414) = - normal(2)*(cf(405)*cf(196)+cf(406)*cf(257))
      cf(415) = - normal(2)*(cf(405)*cf(197)+cf(406)*cf(258))
      cf(416) = (normal(2)*(1.0d0 - gam1*qir*cr2) - tet2)
      cf(417) = gam1*uar2*cr2
      cf(418) = normal(2)*gam1*cr2
      cf(419) = normal(2)*gam1*uar2
      cf(420) = (-normal(3) + (normal(2)*gam1*uar2*cr2))
      cf(421) = gam1*uar3*cr2
      cf(422) = normal(2)*gam1*cr2
      cf(423) = normal(2)*gam1*uar3
      cf(424) = (normal(2)*gam1*uar3*cr2)
      cf(425) = gam1*uar4*cr2
      cf(426) = normal(2)*gam1*uar4
      cf(427) =(dif1*(cf(404)*cf(11)-cf(207))+dif2*(cf(417)*cf(11)
     &    -cf(14))+dif4*(cf(8)+cf(425)*cf(11))+dif3*cf(421)*cf(11))
      cf(428) =(dif1*(cf(404)*cf(12)-cf(208))+dif4*(cf(425)*cf(12)
     *    +cf(9))+dif3*cf(421)*cf(12)+dif2*(cf(417)*cf(12)-cf(15)))
      cf(429) =(dif1*(cf(404)*cf(13)-cf(209))+dif2*(cf(417)*cf(13)
     &    -cf(16))+dif3*cf(421)*cf(13)+dif4*(cf(10)+cf(425)*cf(13)))
      cf(430) =(dif4*(cf(422)*cf(161)+cf(426)*cf(250))+
     &    dif3*(cf(422)*cf(156)+cf(423)*cf(250))+dif2*(cf(419)*cf(250)
     &    +cf(418)*cf(151))+dif1*(cf(407)-cf(210)))
      cf(431) = (dif2*(cf(419)*cf(251)+cf(418)*cf(152))+dif1*(cf(408)
     &    -cf(211))+dif3*cf(423)*cf(251)+dif4*cf(426)*cf(251))
      cf(432) = (dif2*cf(419)*cf(252)+dif1*cf(409)+dif4*cf(426)*cf(252)
     &     +dif3*(cf(422)*cf(157)+cf(423)*cf(252)))
      cf(433) =(dif1*(cf(410)-cf(212))+dif2*cf(419)*cf(253)
     &     +dif3*cf(423)*cf(253)+
     &      dif4*(cf(422)*cf(162)+cf(426)*cf(253)))
      cf(434) =(dif1*cf(411)+dif3*cf(423)*cf(254)
     &      +dif2*cf(419)*cf(254)+dif4*cf(426)*cf(254))
      cf(435) =(dif4*(cf(422)*cf(163)+cf(426)*cf(255))+dif1*(cf(412)
     &      -cf(213))+dif3*(cf(422)*cf(158)+cf(423)*cf(255))+
     &        dif2*(cf(419)*cf(255)+cf(418)*cf(153)))
      cf(436) =(dif2*(cf(419)*cf(256)+cf(418)*cf(154))
     &   +dif1*(cf(413)-cf(214))+dif3*cf(423)*cf(256)
     &   +dif4*cf(426)*cf(256))
      cf(437) =(dif1*cf(414)+dif2*cf(419)*cf(257)+dif3*cf(422)*cf(159)
     &      +dif3*cf(423)*cf(257)+dif4*cf(426)*cf(257))
      cf(438)=(dif1*(cf(415)-cf(215))+dif3*cf(423)*cf(258)+
     &      dif4*(cf(422)*cf(164)+cf(426)*cf(258))+dif2*cf(419)*cf(258))
      cf(439)=(dif2*cf(419)*cf(259)-dif1*normal(2)*cf(406)*cf(259)
     &      +dif3*cf(423)*cf(259)+dif4*cf(426)*cf(259))
      cf(440)=(dif2*cf(419)*cf(249)-dif1*normal(2)*cf(406)*cf(249)+
     &       dif3*cf(423)*cf(249)+dif4*cf(426)*cf(249))
      cf(441)= - (cf(416)+cf(420)*Ugr(2)+cf(424)*Ugr(3))
      cf(442)=(cf(416)+cf(420)*Udr(2)+cf(424)*Udr(3)) 
      cf(443)=(normal(1) + (normal(2)*gam1*uar4*cr2))
      cf(444)=normal(2)*gam1
      cf(445)=-(normal(2)*gam1*cr2)
      cf(446)=vp1*(cf(440)-dif5*cf(444)*cf(249)
     &      +cf(445)*(cf(269)-cf(261)))
      cf(447)=vp1*(cf(427)-dif5*gam1*cr2*cf(11))
      cf(448) = vp1*(cf(428)-dif5*gam1*cr2*cf(12))
      cf(449) = vp1*(cf(429)-dif5*gam1*cr2*cf(13))
      cf(450) = vp1*(cf(430)-dif5*cf(444)*cf(250))
      cf(451) = vp1*(cf(431)-dif5*cf(444)*cf(251))
      cf(452) = vp1*(cf(432)-dif5*cf(444)*cf(252))
      cf(453) = vp1*(cf(433)-dif5*cf(444)*cf(253))
      cf(454) = vp1*(cf(434)-dif5*cf(444)*cf(254))
      cf(455) = vp1*(cf(435)-dif5*cf(444)*cf(255))
      cf(456) = vp1*(cf(436)-dif5*cf(444)*cf(256))
      cf(457) = vp1*(cf(437)-dif5*cf(444)*cf(257))
      cf(458) = vp1*(cf(438)-dif5*cf(444)*cf(258))
      cf(459) = vp1*(cf(439)-dif5*cf(444)*cf(259))
      cf(460) = vp1*(cf(441)-cf(443)*Ugr(4)-cf(445)*cf(260))
      cf(461) = -vp1*(cf(420)*Ugr(1)+cf(445)*cf(262))
      cf(462) = -vp1*(cf(424)*Ugr(1)+cf(445)*cf(263))
      cf(463) = -vp1*(cf(443)*Ugr(1)+cf(445)*cf(264))
      cf(464) = -vp1*cf(445)*cf(136)
      cf(465) = + vp1*(cf(442)+cf(443)*Udr(4)+cf(445)*cf(265))
      cf(466) = vp1*(cf(420)*Udr(1)+cf(445)*cf(266))
      cf(467) = vp1*(cf(424)*Udr(1)+cf(445)*cf(267))
      cf(468) = vp1*(cf(443)*Udr(1)+cf(445)*cf(268))
      cf(469) = vp1*cf(445)*cf(136)
c
      flur3 = vp1*
     &        ((normal(3)*(1.0d0 - cf(406)*cr2) - tet3)*dif1 +
     &        (normal(2)  + (normal(3)*gam1*uar2*cr2))*dif2 +
     &        (-normal(1) + (normal(3)*gam1*uar3*cr2))*dif3 +
     &        (normal(3)*gam1*uar4*cr2)*dif4 -
     &        (normal(3)*gam1*cr2)*dif5)
      cf(470) =((normal(3)*(1.0d0 - cf(406)*cr2) - tet3)*dif1 +
     &         (normal(2)  + (normal(3)*gam1*uar2*cr2))*dif2 +
     &         (-normal(1) + (normal(3)*gam1*uar3*cr2))*dif3 +
     &         (normal(3)*gam1*uar4*cr2)*dif4 -
     &         (normal(3)*gam1*cr2)*dif5)
      cf(471) = (1.0d0 - cf(406)*cr2)
      cf(472) = (cf(405)*cf(190)+cf(406)*cf(250))
      cf(473) = (cf(405)*cf(192)+cf(406)*cf(251))
      cf(474) = (cf(405)*cf(193)+cf(406)*cf(252))
      cf(475) = (cf(405)*cf(194)+cf(406)*cf(253))
      cf(476) = cf(406)*cf(254)
      cf(477) = (cf(405)*cf(191)+cf(406)*cf(255))
      cf(478) = (cf(405)*cf(195)+cf(406)*cf(256))
      cf(479) = (cf(405)*cf(196)+cf(406)*cf(257))
      cf(480) = (cf(405)*cf(197)+cf(406)*cf(258))
      cf(481) = cf(406)*cf(259)
      cf(482) = cf(406)*cf(249) 
      cf(483) = (normal(3)*(1.0d0 - cf(406)*cr2) - tet3)
      cf(484) = normal(3)*cf(405)
      cf(485) = normal(3)*gam1*uar2
      cf(486) = (normal(2) + (normal(3)*cf(417)))
      cf(487) = gam1*uar3*cr2
      cf(488) = normal(3)*cf(405)
      cf(489) = normal(3)*gam1*uar3
      cf(490) = (-normal(1) + (normal(3)*gam1*uar3*cr2))
      cf(491) = gam1*uar4*cr2
      cf(492) = normal(3)*cf(405)
      cf(493) = normal(3)*gam1*uar4
      cf(494) = (normal(3)*gam1*uar4*cr2)
      cf(495) = normal(3)*gam1
      cf(496) = (normal(3)*cf(405))
      cf(497) = vp1*(dif1*(cf(471)*cf(14)-cf(216))-dif5*cf(405)*cf(14)
     &          +dif4*cf(491)*cf(14)+dif2*(cf(11)+cf(417)*cf(14))
     &          +dif3*(cf(487)*cf(14)-cf(8)))
      cf(498) = vp1*(dif1*(cf(471)*cf(15)-cf(217))+dif2*(cf(12)
     &          +cf(417)*cf(15))+dif3*(cf(487)*cf(15)-cf(9))
     &          +dif4*cf(491)*cf(15)-dif5*cf(405)*cf(15))
      cf(499) = vp1*(dif1*(cf(471)*cf(16)-cf(218))+dif2*(cf(13)
     &          +cf(417)*cf(16))+dif3*(cf(487)*cf(16)-cf(10))
     &          +dif4*cf(491)*cf(16)-dif5*cf(405)*cf(16))
      cf(500) =vp1*(dif3*(cf(489)*cf(250)+cf(488)*cf(156))
     &          +dif2*(cf(484)*cf(151)+cf(485)*cf(250))
     &    +dif4*(cf(493)*cf(250)+cf(492)*cf(161))-dif5*cf(495)*cf(250)-
     &    dif1*(normal(3)*cf(472)+cf(219)))
      cf(501) =vp1*(dif2*(cf(484)*cf(152)+cf(485)*cf(251))
     &    +dif4*cf(493)*cf(251)
     &    -dif5*cf(495)*cf(251)+dif3*cf(489)*cf(251)
     &    -dif1*(normal(3)*cf(473)+cf(220)))
      cf(502) =vp1*(dif3*(cf(488)*cf(157)+cf(489)*cf(252))+dif4*cf(493)
     &   *cf(252)-dif5*cf(495)*cf(252)+dif2*cf(485)*cf(252)
     &   -dif1*(normal(3)*cf(474)+cf(221)))
      cf(503) =vp1*(dif3*cf(489)*cf(253)+dif2*cf(485)*cf(253)
     &  +dif4*(cf(492)*cf(162)+cf(493)*cf(253))
     &  -dif5*cf(495)*cf(253)-dif1*normal(3)*cf(475))
      cf(504) = vp1*(dif3*cf(489)*cf(254)+dif2*cf(485)*cf(254)
     &   +dif4*cf(493)*cf(254)-dif5*cf(495)*cf(254)
     &    -dif1*normal(3)*cf(476))
      cf(505) =vp1*(dif3*(cf(488)*cf(158)+cf(489)*cf(255))
     &    +dif2*(cf(484)*cf(153)+cf(485)*cf(255))
     &    +dif4*(cf(492)*cf(163)+cf(493)*cf(255))
     &    -dif5*cf(495)*cf(255)-dif1*(normal(3)*cf(477)+cf(222))) 
      cf(506) = vp1*(dif3*cf(489)*cf(256)+dif2*(cf(484)*cf(154)
     &    +cf(485)*cf(256))+dif4*cf(493)*cf(256)-dif5*cf(495)
     &    *cf(256)-dif1*(normal(3)*cf(478)+cf(223)))
      cf(507) =vp1*(dif3*(cf(488)*cf(159)+cf(489)*cf(257))
     &    +dif2*cf(485)*cf(257)+dif4*cf(493)*cf(257)
     &    -dif5*cf(495)*cf(257)-dif1*(normal(3)*cf(479)+cf(224)))
      cf(508) =vp1*(dif3*cf(489)*cf(258)+dif2*cf(485)*cf(258)
     &    +dif4*(cf(492)*cf(164)+cf(493)*cf(258))
     &    -dif5*cf(495)*cf(258)-dif1*normal(3)*cf(480))
      cf(509) = vp1*(dif2*cf(485)*cf(259)-dif1*normal(3)*cf(481)
     &   +dif4*cf(493)*cf(259)-dif5*cf(495)*cf(259)
     &    +dif3*cf(489)*cf(259))
      cf(510) = vp1*(dif3*cf(489)*cf(249)+dif2*cf(485)*cf(249)
     &  +dif4*cf(493)*cf(249)-dif5*cf(495)*cf(249)
     &   -cf(496)*(cf(269)-cf(261))-dif1*normal(3)*cf(482))
      cf(511) = vp1*(cf(496)*cf(260)-cf(483)-cf(490)
     &    *Ugr(3)-cf(494)*Ugr(4)-cf(486)*Ugr(2))
      cf(512) = vp1*(cf(496)*cf(262)-cf(486)*Ugr(1))
      cf(513) = vp1*(cf(496)*cf(263)-cf(490)*Ugr(1))
      cf(514) = vp1*(cf(496)*cf(264)-cf(494)*Ugr(1))
      cf(515) = vp1*(cf(483)+cf(486)*Udr(2)
     & +cf(490)*Udr(3)+cf(494)*Udr(4)-cf(496)*cf(265))
      cf(516) = vp1*(cf(486)*Udr(1)-cf(496)*cf(266))
      cf(517) = vp1*(cf(490)*Udr(1)-cf(496)*cf(267))
      cf(518) = vp1*(cf(494)*Udr(1)-cf(496)*cf(268))
      cf(519) = vp1*cf(496)*cf(136)
      cf(520) = - vp1*cf(496)*cf(136)
c
      flur4 = vp4*
     &        ((-cr*VdotN   + cf(406))*dif1 +
     &        ( cr*normal(1) - gam1*uar2)*dif2 +
     &        ( cr*normal(2) - gam1*uar3)*dif3 +
     &        ( cr*normal(3) - gam1*uar4)*dif4 +
     &        gam1*dif5)
      cf(521) = ((-cr*VdotN   + cf(406))*dif1 +
     &         ( cr*normal(1) - gam1*uar2)*dif2 +
     &         ( cr*normal(2) - gam1*uar3)*dif3 +
     &         ( cr*normal(3) - gam1*uar4)*dif4 +
     &         gam1*dif5)
      cf(522) = (gam1*cf(190)-VdotN*cf(238)-cr*cf(182))
      cf(523) = (gam1*cf(192)-cr*cf(184)-VdotN*cf(239))
      cf(524) = (gam1*cf(193)-cr*cf(185)-VdotN*cf(240))
      cf(525) = (gam1*cf(194)-cr*cf(186)-VdotN*cf(241))
      cf(526) = (gam1*cf(191)-cr*cf(183)-VdotN*cf(243)) 
      cf(527) = (gam1*cf(195)-cr*cf(187)-VdotN*cf(244))
      cf(528) = (gam1*cf(196)-cr*cf(188)-VdotN*cf(245))
      cf(529) = (gam1*cf(197)-cr*cf(189)-VdotN*cf(246))
      cf(530) = (-cr*VdotN+cf(406))
      cf(531) = (normal(1)*cf(237)*dpstiff
     &      + normal(1)*cf(238)-gam1*cf(151))
      cf(532) = (normal(1)*cf(239)-gam1*cf(152))
      cf(533) = (cr*normal(1)-gam1*uar2)
      cf(554) = (cr*normal(2) - gam1*uar3)
      cf(555) = (cr*normal(3) - gam1*uar4)
      cf(556) =vp4*(dif2*cr*cf(8)-dif1*cr*cf(179)
     &        +dif3*cr*cf(11)+dif4*cr*cf(14))
      cf(557) =vp4*(dif4*cr*cf(15)+dif2*cr*cf(9)
     &        -dif1*cr*cf(180)+dif3*cr*cf(12))
      cf(558) =vp4*(dif3*cr*cf(13)+dif4*cr*cf(16)
     &         +dif2*cr*cf(10)-dif1*cr*cf(181))
      cf(559) =vp4*(dif3*normal(2)*cf(237)+dif4*normal(3)
     &    *cf(237)+gam1*(cf(269)-cf(261))-dif1*VdotN*cf(237))
      cf(560) =vp4*(dif1*cf(522)+dif3*(normal(2)*cf(238)
     &   -gam1*cf(156))+dif4*(normal(3)*cf(238)
     &   -gam1*cf(161))+dif2*cf(531))
      cf(561) =vp4*(dif2*cf(532)+dif3*normal(2)*cf(239)
     &    +dif4*normal(3)*cf(239)+dif1*cf(523))
      cf(562) =vp4*(dif2*normal(1)*cf(240)+dif3*(normal(2)
     &   *cf(240)-gam1*cf(157))+dif4*normal(3)
     &   *cf(240)+dif1*cf(524))
      cf(563) =vp4*(dif2*normal(1)*cf(241)+dif3*normal(2)
     &   *cf(241)+dif4*(normal(3)*cf(241)-gam1*cf(162))
     &   +dif1*cf(525))
      cf(564) =vp4*(dif2*normal(1)*cf(242)+dif3*normal(2)
     &  *cf(242)+dif4*normal(3)*cf(242)-dif1*VdotN*cf(242))
      cf(565) =vp4*(dif2*(normal(1)*cf(243)-gam1*cf(153))+dif3
     &   *(normal(2)*cf(243)-gam1*cf(158))+dif4*(normal(3)
     &   *cf(243)-gam1*cf(163))+dif1*cf(526))
      cf(566) =vp4*(dif1*cf(527)+dif3*normal(2)*cf(244)+dif4
     & *normal(3)*cf(244)+dif2*(normal(1)*cf(244)-gam1*cf(154)))
      cf(567) =vp4*(dif1*cf(528)+dif2*normal(1)*cf(245)+dif4
     & *normal(3)*cf(245)+dif3*(normal(2)*cf(245)-gam1*cf(159)))
      cf(568) =vp4*(dif2*normal(1)*cf(246)+dif3*normal(2)*cf(246)
     &   +dif4*(normal(3)*cf(246)-gam1*cf(164))+dif1*cf(529))
      cf(569) =vp4*(dif2*normal(1)*cf(247)+dif3*normal(2)*cf(247)
     &  +dif4*normal(3)*cf(247)-dif1*VdotN*cf(247))
      cf(570) =-vp4*(cf(530)+cf(533)*Ugr(2)+cf(554)*Ugr(3)
     &   +gam1*cf(260)+cf(555)*Ugr(4)) 
      cf(571) = -vp4*(gam1*cf(262)+cf(533)*Ugr(1))
      cf(572) = -vp4*(cf(554)*Ugr(1)+gam1*cf(263))
      cf(573) = -vp4*(cf(555)*Ugr(1)+gam1*cf(264))
      cf(574) =vp4*(cf(530)+cf(533)*Udr(2)+cf(554)
     &   *Udr(3)+cf(555)*Udr(4)+gam1*cf(265))
      cf(575) = vp4*(cf(533)*Udr(1)+gam1*cf(266))
      cf(576) = vp4*(cf(554)*Udr(1)+gam1*cf(267))
      cf(577) = vp4*(cf(555)*Udr(1)+gam1*cf(268))
      cf(578) = -vp4*gam1*cf(136) 
      cf(579) = vp4*gam1*cf(136)
c
      flur5 = vp5*
     &        (( cr*VdotN  + gam1* qir)*dif1 +
     &        (-cr*normal(1) - gam1*uar2)*dif2 +
     &        (-cr*normal(2) - gam1*uar3)*dif3 +
     &        (-cr*normal(3) - gam1*uar4)*dif4 +
     &        gam1*dif5)
      cf(580) = (( cr*VdotN  + gam1* qir)*dif1 +
     &         (-cr*normal(1) - gam1*uar2)*dif2 +
     &         (-cr*normal(2) - gam1*uar3)*dif3 +
     &         (-cr*normal(3) - gam1*uar4)*dif4 +
     &         gam1*dif5)
      cf(581) = (VdotN*cf(238)+cr*cf(182)+gam1*cf(190))
      cf(582) = (VdotN*cf(239)+cr*cf(184)+gam1*cf(192))
      cf(583) = (VdotN*cf(240)+cr*cf(185)+gam1*cf(193))
      cf(584) = (VdotN*cf(241)+gam1*cf(194)+cr*cf(186))
      cf(585) = (gam1*cf(191)+cr*cf(183)+VdotN*cf(243))
      cf(586) = (cr*cf(187)+gam1*cf(195)+VdotN*cf(244))
      cf(587) = (VdotN*cf(245)+cr*cf(188)+gam1*cf(196))
      cf(588) = (VdotN*cf(246)+cr*cf(189)+gam1*cf(197))
      cf(589) = (cr*VdotN  + gam1* qir)
      cf(590) = (-cr*normal(1) - gam1*uar2)
      cf(591) = (-cr*normal(2) - gam1*uar3)
      cf(592) = (-cr*normal(3) - gam1*uar4)
      cf(593) =vp5*(dif1*cr*cf(179)-dif2*cr*cf(8)-dif3
     &      *cr*cf(11)-dif4*cr*cf(14))
      cf(594) =vp5*(dif1*cr*cf(180)-dif2*cr*cf(9)-dif3
     &      *cr*cf(12)-dif4*cr*cf(15))
      cf(595) =vp5*(dif1*cr*cf(181)-dif2*cr*cf(10)-dif3
     &      *cr*cf(13)-dif4*cr*cf(16))
      cf(596) =vp5*(dif1*VdotN*cf(237)-dif2*normal(1)*cf(237)
     &   *dpstiff-dif3*normal(2)*cf(237)-dif4*normal(3)*cf(237)
     &   +gam1*(cf(269)-cf(261)))
      cf(597)=vp5*(dif1*cf(581)-dif2*(normal(1)*cf(238)+gam1*cf(151))
     &   -dif3*(normal(2)*cf(238)+gam1*cf(156))-dif4*(normal(3)
     &   *cf(238)+gam1*cf(161))) 
      cf(598) =vp5*(dif1*cf(582)-dif2*(normal(1)*cf(239)+gam1
     &   *cf(152))-dif4*normal(3)*cf(239)-dif3*normal(2)*cf(239))
      cf(599) =-vp5*(dif2*normal(1)*cf(240)+dif3*(normal(2)*cf(240)
     &   +gam1*cf(157))+dif4*normal(3)*cf(240)-dif1*cf(583)) 
      cf(600) = -vp5*(dif2*normal(1)*cf(241)+dif3*normal(2)*cf(241)
     &   +dif4*(normal(3)*cf(241)+gam1*cf(162))-dif1*cf(584))
      cf(601) =vp5*(dif1*VdotN*cf(242)-dif2*normal(1)*cf(242)-dif3
     &   *normal(2)*cf(242)-dif4*normal(3)*cf(242))
      cf(602)=vp5*(dif1*cf(585)-dif2*(normal(1)*cf(243)+gam1*cf(153))
     &   -dif3*(normal(2)*cf(243)+gam1*cf(158))-dif4*(normal(3)*cf(243)
     &    +gam1*cf(163)))
      cf(603) =vp5*(dif1*cf(586)-dif4*normal(3)*cf(244)-dif3*normal(2)
     &   *cf(244)-dif2*(normal(1)*cf(244)+gam1*cf(154))) 
      cf(604) =vp5*(dif1*cf(587)-dif2*normal(1)*cf(245)-dif3*(normal(2)
     &   *cf(245)+gam1*cf(159))-dif4*normal(3)*cf(245))
      cf(605) =vp5*(dif1*cf(588)-dif2*normal(1)*cf(246)-dif3*normal(2)
     &   *cf(246)-dif4*(normal(3)*cf(246)+gam1*cf(164)))
      cf(606) =vp5*(dif1*VdotN*cf(247)-dif2*normal(1)*cf(247)-dif4
     &   *normal(3)*cf(247)-dif3*normal(2)*cf(247))
      cf(607) =-vp5*(cf(590)*Ugr(2)+cf(591)*Ugr(3)+cf(592)*Ugr(4)
     &   +gam1*cf(260)+ cf(589)) 
      cf(608) = -vp5*(gam1*cf(262)+cf(590)*Ugr(1))
      cf(609) = - vp5*(gam1*cf(263)+cf(591)*Ugr(1))
      cf(610) = - vp5*(cf(592)*Ugr(1)+gam1*cf(264))
      cf(611) =vp5*(cf(589)+cf(590)*Udr(2)+cf(591)*Udr(3)
     &   +cf(592)*Udr(4)+gam1*cf(265))
      cf(612) = vp5*(cf(590)*Udr(1)+gam1*cf(266))
      cf(613) = vp5*(cf(591)*Udr(1)+gam1*cf(267))
      cf(614) = vp5*(cf(592)*Udr(1)+gam1*cf(268))
      cf(615) = - vp5*gam1*cf(136)
      cf(616) = vp5*gam1*cf(136) 
c
      phi(1) = phi(1) - (normal(1)*flur1 +
     &         normal(2)*flur2 + normal(3)*flur3 +
     &         0.5d0*(flur4 + flur5)*cr2)
      cf(617) = 0.5d0*(flur4 + flur5)
      cf(618) =(0.5d0*cr2*(cf(556)+cf(593))+flur1*cf(8)+flur2*cf(11)
     &   +normal(1)*cf(379)+flur3*cf(14)+normal(3)*cf(497)
     &   +normal(2)*cf(447))
      cf(619) = (normal(1)*cf(315)+normal(2)*cf(403)+normal(3)*cf(470))
      cf(620) =(flur1*cf(9)+flur2*cf(12)+normal(2)*cf(448)+flur3*cf(15)
     &   +normal(3)*cf(498)+0.5d0*cr2*(cf(557)+cf(594))
     &   +normal(1)*cf(380))
      cf(621) =(flur1*cf(10)+flur2*cf(13)+normal(2)*cf(449)+flur3
     &   *cf(16)+normal(3)*cf(499)+0.5d0*cr2*(cf(558)+cf(595))
     &   +normal(1)*cf(381))
      cf(622)=(normal(2)*cf(446)+normal(3)*cf(510)+0.5d0*cr2*(cf(596)
     &   +cf(559))+cf(617)*cf(249)+normal(1)*cf(382))
      cf(623) =(normal(2)*cf(450)+normal(3)*cf(500)+0.5d0*cr2
     &   *(cf(560)+cf(597))+cf(617)*cf(250)+normal(1)*cf(383)) 
      cf(624) =(cf(617)*cf(251)+normal(2)*cf(451)+normal(3)*cf(501)
     &   +0.5d0*cr2*(cf(561)+cf(598))+normal(1)*cf(384))
      cf(625) =(normal(2)*cf(452)+normal(3)*cf(502)+0.5d0*cr2
     &   *(cf(562)+cf(599))+cf(617)*cf(252)+normal(1)*cf(385))
      cf(626) =(normal(2)*cf(453)+normal(3)*cf(503)+0.5d0*cr2
     &   *(cf(563)+cf(600))+cf(617)*cf(253)+normal(1)*cf(386))
      cf(627) =(normal(2)*cf(454)+normal(3)*cf(504)+0.5d0*cr2
     &   *(cf(564)+cf(601))+cf(617)*cf(254)+normal(1)*cf(387))
      cf(628) =(normal(2)*cf(455)+normal(3)*cf(505)+0.5d0*cr2
     &   *(cf(565)+cf(602))+cf(617)*cf(255)+normal(1)*cf(388))
      cf(629)=(normal(2)*cf(456)+normal(3)*cf(506)+0.5d0*cr2
     &   *(cf(603)+cf(566))+cf(617)*cf(256)+normal(1)*cf(389)) 
      cf(630)=(normal(2)*cf(457)+normal(3)*cf(507)+0.5d0*cr2*(cf(604)
     &   +cf(567))+cf(617)*cf(257)+normal(1)*cf(390))
      cf(631)=(normal(2)*cf(458)+normal(3)*cf(508)+0.5d0*cr2*(cf(605)
     &   +cf(568))+cf(617)*cf(258)+normal(1)*cf(391))
      cf(632)=(normal(2)*cf(459)+normal(3)*cf(509)+0.5d0*cr2
     &   *(cf(569)+cf(606))+cf(617)*cf(259)+normal(1)*cf(392))
      cf(633)=(normal(2)*cf(465)+normal(3)*cf(515)+0.5d0*cr2*(cf(574)
     &    +cf(611))+normal(1)*cf(393))
      cf(634)=(normal(1)*cf(394)+normal(3)*cf(516)+0.5d0*cr2*(cf(575)
     &    +cf(612))+normal(2)*cf(466))
      cf(635)=(normal(2)*cf(467)+normal(3)*cf(517)+0.5d0*cr2*(cf(576)
     &    +cf(613))+normal(1)*cf(395))
      cf(636)=(normal(2)*cf(468)+normal(3)*cf(518)+0.5d0*cr2*(cf(577)
     &    +cf(614))+normal(1)*cf(396))
      cf(637)=(normal(2)*cf(469)+normal(3)*cf(520)+0.5d0*cr2*(cf(579)
     &    +cf(616))+normal(1)*cf(397))
      cf(638)=(normal(2)*cf(460)+normal(3)*cf(511)+0.5d0*cr2*(cf(607)
     &    +cf(570))+normal(1)*cf(398))
      cf(639)=(normal(2)*cf(461)+normal(3)*cf(512)+0.5d0*cr2*(cf(571)
     &    +cf(608))+normal(1)*cf(399))
      cf(640)=(normal(1)*cf(400)+normal(3)*cf(513)+0.5d0*cr2*(cf(572)
     &    +cf(609))+normal(2)*cf(462))
      cf(641)=(0.5d0*cr2*(cf(573)+cf(610))+normal(2)*cf(463)+
     &    normal(3)*cf(514)+normal(1)*cf(401))
      cf(642)=(normal(1)*cf(402)+normal(2)*cf(464)+normal(3)*cf(519)
     &    +0.5d0*cr2*(cf(578)+cf(615)))
c
      phi(2) = phi(2) - (  
     &         (uar2*normal(1))*flur1 +
     &         (uar2*normal(2) - normal(3))*flur2 +
     &         (uar2*normal(3) + normal(2))*flur3 +
     &         0.5d0*normal(1)*(flur4 - flur5)/cr +
     &         0.5d0*uar2*(flur4 + flur5)*cr2)
      cf(643)=(uar2*normal(1))
      cf(644)=(uar2*normal(2) - normal(3))
      cf(645)=(uar2*normal(3) + normal(2))
      cf(646)=0.5d0*(flur4 - flur5)/cr
      cf(647)=0.5d0*normal(1)/cr
      cf(648)=-0.5d0*normal(1)*(flur4 - flur5)/(cr*cr)
      cf(649)=0.5d0*(flur4 + flur5)*cr2
      cf(650)=0.5d0*uar2*cr2
      cf(651)=0.5d0*uar2*(flur4 + flur5)
      cf(652) =(flur1*uar2*cf(8)+cf(643)*cf(379)+flur2*uar2*cf(11)
     %     -flur2*cf(14)+cf(644)*cf(447)+flur3*uar2*cf(14)+flur3*cf(11)
     &     +cf(645)*cf(497)+cf(646)*cf(8)+cf(647)*(cf(556)-cf(593))
     &     +cf(650)*(cf(556)+cf(593)))
      cf(653)=(flur1*uar2*cf(9)+cf(643)*cf(380)+flur2*uar2*cf(12)
     &     -flur2*cf(15)+cf(644)*cf(448)+flur3*uar2*cf(15)+flur3
     &     *cf(12)+cf(645)*cf(498)+cf(646)*cf(9)+cf(647)*(cf(557)
     &     -cf(594))+cf(650)*(cf(557)+cf(594)))
      cf(654)=(flur1*uar2*cf(10)+cf(643)*cf(381)+flur2*uar2*cf(13)
     &     -flur2*cf(16)+cf(644)*cf(449)+flur3*uar2*cf(16)+flur3
     &    *cf(13)+cf(645)*cf(499)+cf(646)*cf(10)+cf(647)*(cf(558)
     &    -cf(595))+ cf(650)*(cf(595)+cf(558)))
      cf(655)=(flur1*normal(1)*cf(151)+flur2*normal(2)*cf(151)
     &    +cf(644)*cf(450)+flur3*normal(3)*cf(151)+cf(645)*cf(500)
     &    +cf(647)*(cf(560)-cf(597))+cf(648)*cf(238)+cf(649)*cf(151)
     &    +cf(650)*(cf(597)+cf(560))+cf(651)*cf(250)+cf(643)*cf(383))
      cf(656)=(cf(643)*cf(384)+flur2*normal(2)*cf(152)+cf(644)*cf(451)
     &    +flur3*normal(3)*cf(152)+cf(645)*cf(501)+cf(647)*(cf(561)
     &    -cf(598))+cf(648)*cf(239)+cf(649)*cf(152)+cf(650)*(cf(598)
     &    +cf(561))+cf(651)*cf(251)+flur1*normal(1)*cf(152))
      cf(657)=(cf(644)*cf(452)+cf(645)*cf(502)+cf(647)*(cf(562)-cf(599))
     &    +cf(648)*cf(240)+cf(650)*(cf(599)+cf(562))+cf(651)*cf(252)
     &    +cf(643)*cf(385))
      cf(658)=(cf(644)*cf(453)+cf(645)*cf(503)+cf(647)*(cf(563)-cf(600))
     &    +cf(648)*cf(241)+cf(650)*(cf(600)+cf(563))+cf(651)*cf(253)
     &    +cf(643)*cf(386))
      cf(659)=(cf(644)*cf(454)+cf(645)*cf(504)+cf(647)*(cf(564)-cf(601))
     &    +cf(648)*cf(242)+cf(650)*(cf(601)+cf(564))+cf(651)*cf(254)
     &    +cf(643)*cf(387))
      cf(660)=(cf(643)*cf(388)+flur2*normal(2)*cf(153)+cf(644)*cf(455)
     &    +flur3*normal(3)*cf(153)+cf(645)*cf(505)+cf(647)*(cf(565)
     &    -cf(602))+cf(648)*cf(243)+cf(649)*cf(153)+cf(650)*(cf(602)
     &    +cf(565))+cf(651)*cf(255)+flur1*normal(1)*cf(153))
      cf(661)=(cf(643)*cf(389)+flur2*normal(2)*cf(154)+cf(644)*cf(456)
     &    +flur3*normal(3)*cf(154)+cf(645)*cf(506)+cf(647)*(cf(566)
     &    -cf(603))+cf(648)*cf(244)+cf(649)*cf(154))
      cf(662)=(cf(661)+cf(650)*(cf(603)+cf(566))+cf(651)*cf(256)+flur1
     &    *normal(1)*cf(154))
      cf(663)=(cf(644)*cf(457)+cf(645)*cf(507)+cf(647)*(cf(567)-cf(604))
     &    +cf(648)*cf(245)+cf(650)*(cf(604)+cf(567))+cf(651)*cf(257)
     &    +cf(643)*cf(390))
      cf(664)=(cf(644)*cf(458)+cf(645)*cf(508)+cf(647)*(cf(568)-cf(605))
     &    +cf(648)*cf(246)+cf(650)*(cf(605)+cf(568))+cf(651)*cf(258)
     &    +cf(643)*cf(391))
      cf(665)=(cf(644)*cf(459)+cf(645)*cf(509)+cf(647)*(cf(569)-cf(606))
     &    +cf(648)*cf(247)+cf(650)*(cf(606)+cf(569))+cf(651)*cf(259)
     &    +cf(643)*cf(392))
      cf(666)=(cf(643)*cf(315)+cf(645)*cf(470)+cf(644)*cf(403))
      cf(667)=(cf(644)*cf(446)+cf(645)*cf(510)+cf(647)*(cf(559)-cf(596))
     &    +cf(648)*cf(237)+cf(650)*(cf(559)+cf(596))+cf(651)*cf(249)
     &    +cf(643)*cf(382))
      cf(668)=(cf(644)*cf(465)+cf(645)*cf(515)+cf(647)*(cf(574)-cf(611))
     &    +cf(650)*(cf(611)+cf(574))+cf(643)*cf(393))
      cf(669)=(cf(644)*cf(466)+cf(645)*cf(516)+cf(647)*(cf(575)-cf(612))
     &    +cf(650)*(cf(612)+cf(575))+cf(643)*cf(394))
      cf(670)=(cf(644)*cf(467)+cf(645)*cf(517)+cf(647)*(cf(576)-cf(613))
     &    +cf(650)*(cf(576)+cf(613))+cf(643)*cf(395))
      cf(671)=(cf(644)*cf(468)+cf(645)*cf(518)+cf(647)*(cf(577)-cf(614))
     &    +cf(650)*(cf(577)+cf(614))+cf(643)*cf(396))
      cf(672)=(cf(644)*cf(469)+cf(645)*cf(520)+cf(647)*(cf(579)-cf(616))
     &    +cf(650)*(cf(579)+cf(616))+cf(643)*cf(397))
      cf(673)=(cf(643)*cf(398)+cf(645)*cf(511)+cf(647)*(cf(570)-cf(607))
     &    +cf(650)*(cf(607)+cf(570))+cf(644)*cf(460))
      cf(674)=(cf(643)*cf(399)+cf(645)*cf(512)+cf(647)*(cf(571)-cf(608))
     &    +cf(650)*(cf(608)+cf(571))+cf(644)*cf(461))
      cf(675)=(cf(644)*cf(462)+cf(645)*cf(513)+cf(647)*(cf(572)-cf(609))
     &    +cf(650)*(cf(609)+cf(572))+cf(643)*cf(400))
      cf(676)=(cf(644)*cf(463)+cf(645)*cf(514)+cf(647)*(cf(573)-cf(610))
     &    +cf(650)*(cf(610)+cf(573))+cf(643)*cf(401))
      cf(677)=(cf(643)*cf(402)+cf(644)*cf(464)+cf(647)*(cf(578)-cf(615))
     &    +cf(645)*cf(519)+cf(650)*(cf(615)+cf(578)))
      cf(678)=(cf(647)*cf(521)+cf(650)*cf(521))
      cf(679)=(cf(650)*cf(580)-cf(647)*cf(580))
      cf(680)=(cf(94)-cf(652))
      cf(681)=(cf(95)-cf(653))
      cf(682)=(cf(96)-cf(654))
      cf(683)=(cf(51)-cf(655))
      cf(684)=(cf(52)+cf(56)-cf(656))
      cf(685)=(cf(53)-cf(657))
      cf(686)=(cf(54)-cf(658))
      cf(687)=(normal(1)-cf(659))
      cf(688)=(cf(98)-cf(660))
      cf(689)=(cf(99)-cf(662)+cf(100))
      cf(690)=(cf(101)-cf(663))
      cf(691)=(cf(102)-cf(664))
      cf(692)=(normal(1)-cf(665))
      cf(693)=(uar3*normal(1) + normal(3))
      cf(694)=(uar3*normal(2))
c
      phi(3) = phi(3) - (
     &         (uar3*normal(1) + normal(3))*flur1 +
     &         (uar3*normal(2))*flur2 +
     &         (uar3*normal(3) - normal(1))*flur3 +
     &         0.5d0*normal(2)*(flur4 - flur5)/cr +
     &         0.5d0*uar3*(flur4 + flur5)*cr2)
      cf(695) = (uar3*normal(3) - normal(1))
      cf(696) = 0.5d0*(flur4 - flur5)/cr
      cf(697) = 0.5d0*normal(2)/cr
      cf(698) = - 0.5d0*normal(2)*(flur4 - flur5)/(cr*cr)
      cf(699) =(flur1*uar3*cf(8)+flur1*cf(14)+flur2*uar3*cf(11)
     &     +cf(694)*cf(447)+flur3*uar3*cf(14)-flur3*cf(8)+cf(695)
     &     *cf(497)+cf(696)*cf(11)+cf(697)*(cf(556)-cf(593))
     &     +cf(693)*cf(379))
      cf(700) = 0.5d0*cr2*(flur4 + flur5)
      cf(701)=0.5d0*uar3*cr2
      cf(702)=0.5d0*uar3*(flur4+flur5)
      cf(703)=(cf(701)*(cf(593)+cf(556))+cf(699)) 
      cf(704)=(flur1*cf(15)+cf(693)*cf(380)+flur2*uar3*cf(12)
     &     +cf(694)*cf(448)+flur3*uar3*cf(15)-flur3*cf(9)
     &     +cf(695)*cf(498)+cf(696)*cf(12)+cf(697)*(cf(557)
     &     -cf(594))+cf(701)*(cf(594)+cf(557))+flur1*uar3*cf(9))
      cf(705)=(flur1*uar3*cf(10)+flur1*cf(16)+flur2*uar3*cf(13)
     &     +cf(694)*cf(449)+flur3*uar3*cf(16)-flur3*cf(10)
     &     +cf(695)*cf(499)+cf(696)*cf(13)+cf(697)*(cf(558)
     &     -cf(595))+cf(701)*(cf(595)+cf(558))+cf(693)*cf(381))
      cf(706)=(cf(694)*cf(403)+cf(695)*cf(470)+cf(693)*cf(315))
      cf(707)=(cf(694)*cf(446)+cf(695)*cf(510)+cf(697)*(cf(559)-cf(596))
     &     +cf(698)*cf(237)+cf(701)*(cf(596)+cf(559))+cf(702)*cf(249)
     &     +cf(693)*cf(382))
      cf(708)=(flur1*normal(1)*cf(156)+flur2*normal(2)*cf(156)
     &     +cf(694)*cf(450)+flur3*normal(3)*cf(156)+cf(695)*cf(500)
     &     +cf(697)*(cf(560)-cf(597))+cf(698)*cf(238)+cf(700)*cf(156)
     &     +cf(702)*cf(250)+cf(701)*(cf(597)+cf(560))+cf(693)*cf(383))
      cf(709)=(cf(694)*cf(451)+cf(695)*cf(501)+cf(697)*(cf(561)-cf(598))
     &     +cf(698)*cf(239)+cf(702)*cf(251)+cf(701)*(cf(598)+cf(561))
     &     +cf(693)*cf(384))
      cf(710)=(flur1*normal(1)*cf(157)+flur2*normal(2)*cf(157)+cf(694)
     &     *cf(452)+flur3*normal(3)*cf(157)+cf(695)*cf(502)+cf(697)
     &     *(cf(562)-cf(599))+cf(698)*cf(240)+cf(701)*(cf(599)
     &     +cf(562))+cf(702)*cf(252)+cf(700)*cf(157)+cf(693)*cf(385))
      cf(711)=(cf(694)*cf(453)+cf(695)*cf(503)+cf(697)*(cf(563)-cf(600))
     &     +cf(698)*cf(241)+cf(702)*cf(253)+cf(701)*(cf(600)+cf(563))
     &     +cf(693)*cf(386))
      cf(712)=(cf(694)*cf(454)+cf(695)*cf(504)+cf(697)*(cf(564)-cf(601))
     &     +cf(698)*cf(242)+cf(702)*cf(254)+cf(701)*(cf(601)+cf(564))
     &     +cf(693)*cf(387))
      cf(713)=(cf(701)*(cf(602)+cf(565))+cf(700)*cf(158)+flur2*normal(2)
     &     *cf(158)+flur1*normal(1)*cf(158)+cf(694)*cf(455)+flur3
     &     *normal(3)*cf(158)+cf(695)*cf(505)+cf(697)*(cf(565)-cf(602))
     &     +cf(698)*cf(243)+cf(693)*cf(388)+cf(702)*cf(255))
      cf(714)=(cf(694)*cf(456)+cf(695)*cf(506)+cf(697)*(cf(566)-cf(603))
     &     +cf(698)*cf(244)+cf(702)*cf(256)+cf(701)*(cf(603)+cf(566))
     &     +cf(693)*cf(389))
      cf(715)=(flur1*normal(1)*cf(159)+cf(694)*cf(457)+flur2*normal(2)
     &     *cf(159)+flur3*normal(3)*cf(159)+cf(695)*cf(507)+cf(697)
     &     *(cf(567)-cf(604))+cf(698)*cf(245)+cf(702)*cf(257)+cf(701)
     &     *(cf(604)+cf(567))+cf(700)*cf(159)+cf(693)*cf(390))
      cf(716)=(cf(694)*cf(458)+cf(695)*cf(508)+cf(693)*cf(391)+cf(698)
     &     *cf(246)+cf(697)*(cf(568)-cf(605))+cf(702)*cf(258)+cf(701)
     &     *(cf(605)+cf(568)))
      cf(717)=(cf(694)*cf(459)+cf(695)*cf(509)+cf(697)*(cf(569)-cf(606))
     &     +cf(698)*cf(247)+cf(701)*(cf(606)+cf(569))+cf(702)*cf(259)
     &     +cf(693)*cf(392))
      cf(718)=(cf(694)*cf(465)+cf(695)*cf(515)+cf(697)*(cf(574)-cf(611))
     &     +cf(701)*(cf(611)+cf(574))+cf(693)*cf(393))
      cf(719)=(cf(694)*cf(466)+cf(695)*cf(516)+cf(697)*(cf(575)-cf(612))
     &     +cf(693)*cf(394)+cf(701)*(cf(612)+cf(575)))
      cf(720)=(cf(694)*cf(467)+cf(695)*cf(517)+cf(697)*(cf(576)-cf(613))
     &     +cf(693)*cf(395)+cf(701)*(cf(613)+cf(576)))
      cf(721)=(cf(694)*cf(468)+cf(695)*cf(518)+cf(697)*(cf(577)-cf(614))
     &     +cf(693)*cf(396)+cf(701)*(cf(614)+cf(577)))
      cf(722)=(cf(697)*(cf(579)-cf(616))+cf(695)*cf(520)+cf(694)*cf(469)
     &     +cf(693)*cf(397)+cf(701)*(cf(616)+cf(579)))
      cf(723)=(cf(694)*cf(460)+cf(695)*cf(511)+cf(697)*(cf(570)-cf(607))
     &     +cf(693)*cf(398)+cf(701)*(cf(607)+cf(570)))
      cf(724)=(cf(697)*(cf(571)-cf(608))+cf(694)*cf(461)+cf(695)*cf(512)
     &     +cf(693)*cf(399)+cf(701)*(cf(608)+cf(571)))
      cf(725)=(cf(697)*(cf(572)-cf(609))+cf(695)*cf(513)+cf(694)*cf(462)
     &     +cf(693)*cf(400)+cf(701)*(cf(609)+cf(572)))
      cf(726)=(cf(694)*cf(463)+cf(693)*cf(401)+cf(697)*(cf(573)-cf(610))
     &     +cf(695)*cf(514)+cf(701)*(cf(610)+cf(573)))
      cf(727)=(cf(693)*cf(402)+cf(694)*cf(464)+cf(695)*cf(519)+cf(697)
     &     *(cf(578)-cf(615))+cf(701)*(cf(615)+cf(578)))
c
      phi(4) = phi(4) - ( 
     &         (uar4*normal(1) - normal(2))*flur1 +
     &         (uar4*normal(2) + normal(1))*flur2 +
     &         (uar4*normal(3))*flur3 +
     &         0.5d0*normal(3)*(flur4 - flur5)/cr +
     &         0.5d0*uar4*(flur4 + flur5)*cr2)
      cf(728)=(uar4*normal(1) - normal(2))
      cf(729)=(uar4*normal(2) + normal(1))
      cf(730)=(uar4*normal(3))
      cf(731)=0.5d0*(flur4 - flur5)/cr
      cf(732)=0.5d0*normal(3)/cr
      cf(733)=-0.5d0*normal(3)*(flur4 - flur5)/(cr*cr)
      cf(734)=0.5d0*(flur4 + flur5)*cr2
      cf(735)=0.5d0*uar4*cr2
      cf(736)=0.5d0*uar4*(flur4 + flur5)
      cf(737)=(flur1*normal(1)*cf(161)+cf(728)*cf(383)+flur2
     &   *normal(2)*cf(161)+flur3*normal(3)*cf(161)+cf(729)
     &   *cf(450)+cf(730)*cf(500)+cf(732)*(cf(560)-cf(597))
     &   +cf(733)*cf(238)+cf(734)*cf(161)+cf(735)*(cf(597)
     &   +cf(560))+cf(736)*cf(250))
      cf(738)=(cf(735)*(cf(598)+cf(561))+cf(733)*cf(239)+cf(732)
     &   *(cf(561)-cf(598))+cf(730)*cf(501)+cf(729)*cf(451)
     &   +cf(728)*cf(384)+cf(736)*cf(251))
      cf(739)=(cf(735)*(cf(599)+cf(562))+cf(733)*cf(240)+cf(732)
     &   *(cf(562)-cf(599))+cf(730)*cf(502)+cf(729)*cf(452)
     &   +cf(728)*cf(385)+cf(736)*cf(252))
      cf(740)=(cf(735)*(cf(600)+cf(563))+cf(734)*cf(162)+cf(733)
     &   *cf(241)+cf(732)*(cf(563)-cf(600))+cf(730)*cf(503)+flur3
     &   *normal(3)*cf(162)+cf(729)*cf(453)+flur2*normal(2)
     &   *cf(162)+cf(728)*cf(386)+flur1*normal(1)*cf(162)
     &   +cf(736)*cf(253))
      cf(741)=(cf(736)*cf(254)+cf(735)*(cf(601)+cf(564))+cf(733)
     &   *cf(242)+cf(732)*(cf(564)-cf(601))+cf(730)*cf(504)
     &   +cf(729)*cf(454)+cf(728)*cf(387))
      cf(742)=(flur1*normal(1)*cf(163)+flur2*normal(2)*cf(163)
     &   +cf(728)*cf(388)+cf(729)*cf(455)+flur3*normal(3)*cf(163)
     &   +cf(730)*cf(505)+cf(732)*(cf(565)-cf(602))+cf(733)*cf(243)
     &   +cf(734)*cf(163)+cf(735)*(cf(602)+cf(565))+cf(736)*cf(255))
      cf(743)=(cf(735)*(cf(603)+cf(566))+cf(733)*cf(244)+cf(732)
     &   *(cf(566)-cf(603))+cf(730)*cf(506)+cf(729)*cf(456)+cf(728)
     &   *cf(389)+cf(736)*cf(256))
      cf(744)=(cf(736)*cf(257)+cf(733)*cf(245)+cf(732)*(cf(567)
     &   -cf(604))+cf(730)*cf(507)+cf(729)*cf(457)+cf(728)*cf(390)
     &   +cf(735)*(cf(604)+cf(567)))
      cf(745)=(cf(735)*(cf(605)+cf(568))+cf(734)*cf(164)+cf(733)*cf(246)
     &   +cf(732)*(cf(568)-cf(605))+cf(730)*cf(508)+flur3*normal(3)
     &   *cf(164)+cf(729)*cf(458)+flur2*normal(2)*cf(164)+cf(728)
     &   *cf(391)+flur1*normal(1)*cf(164)+cf(736)*cf(258))
      cf(746)=(cf(736)*cf(259)+cf(735)*(cf(606)+cf(569))+cf(733)*cf(247)
     &   +cf(732)*(cf(569)-cf(606))+cf(730)*cf(509)+cf(729)*cf(459)
     &   +cf(728)*cf(392))
      cf(747)=(cf(729)*cf(446)+cf(730)*cf(510)+cf(732)*(cf(559)-cf(596))
     &   +cf(733)*cf(237)+cf(735)*(cf(596)+cf(559))+cf(736)*cf(249)
     &   +cf(728)*cf(382))
      cf(748)=(cf(730)*cf(515)+cf(729)*cf(465)+cf(732)*(cf(574)-cf(611))
     &   +cf(735)*(cf(574)+cf(611))+cf(728)*cf(393))
      cf(749)=(cf(728)*cf(394)+cf(730)*cf(516)+cf(729)*cf(466)
     &   +cf(732)*(cf(575)-cf(612))+cf(735)*(cf(612)+cf(575)))
      cf(750)=(cf(732)*(cf(576)-cf(613))+cf(730)*cf(517)+cf(729)*cf(467)
     &   +cf(728)*cf(395)+cf(735)*(cf(613)+cf(576)))
      cf(751)=(cf(732)*(cf(577)-cf(614))+cf(730)*cf(518)+cf(729)*cf(468)
     &   +cf(728)*cf(396)+cf(735)*(cf(614)+cf(577)))
      cf(752)=(cf(735)*(cf(616)+cf(579))+cf(732)*(cf(579)-cf(616))
     &   +cf(730)*cf(520)+cf(729)*cf(469)+cf(728)*cf(397))
      cf(753)=(cf(729)*cf(460)+cf(730)*cf(511)+cf(732)*(cf(570)-cf(607))
     &   +cf(735)*(cf(570)+cf(607))+cf(728)*cf(398))
      cf(754)=(cf(729)*cf(461)+cf(730)*cf(512)+cf(732)*(cf(571)-cf(608))
     &   +cf(735)*(cf(608)+cf(571))+cf(728)*cf(399))
      cf(755)=(cf(729)*cf(462)+cf(730)*cf(513)+cf(732)*(cf(572)-cf(609))
     &   +cf(735)*(cf(609)+cf(572))+cf(728)*cf(400))
      cf(756)=(cf(729)*cf(463)+cf(730)*cf(514)+cf(732)*(cf(573)-cf(610))
     &   +cf(735)*(cf(610)+cf(573))+cf(728)*cf(401))
      cf(757)=(cf(728)*cf(402)+cf(729)*cf(464)+cf(730)*cf(519)+cf(732)
     &   *(cf(578)-cf(615))+cf(735)*(cf(615)+cf(578)))
      cf(758)=(flur1*uar4*cf(8)+cf(728)*cf(379)+flur2*uar4*cf(11)+flur2
     &   *cf(8)+cf(729)*cf(447)+flur3*uar4*cf(14)+cf(730)*cf(497)
     &   +cf(731)*cf(14)+cf(732)*(cf(556)-cf(593))+cf(735)*(cf(556)
     &   +cf(593))-flur1*cf(11))
      cf(759)=(flur1*uar4*cf(9)+cf(728)*cf(380)+flur2*uar4*cf(12)+flur2
     &   *cf(9)+cf(729)*cf(448)+flur3*uar4*cf(15)+cf(730)*cf(498)
     &   +cf(731)*cf(15)+cf(732)*(cf(557)-cf(594))+cf(735)*(cf(594)
     &   +cf(557))-flur1*cf(12))
      cf(760)=(flur1*uar4*cf(10)-flur1*cf(13)+cf(728)*cf(381)+flur2
     &   *uar4*cf(13)+flur2*cf(10)+cf(729)*cf(449)+flur3*uar4*cf(16)
     &   +cf(730)*cf(499)+cf(731)*cf(16)+cf(732)*(cf(558)-cf(595))
     &   +cf(735)*(cf(595)+cf(558)))
      cf(761)=(cf(729)*cf(403)+cf(730)*cf(470)+cf(728)*cf(315))
      cf(762)=-(cf(732)*cf(521)+cf(735)*cf(521))
      cf(763)=(cf(732)*cf(580)-cf(735)*cf(580))
      cf(764)=(cf(115)-cf(758))
      cf(765)=(cf(116)-cf(759))
      cf(766)=(cf(117)-cf(760))
c
      phi(5) = phi(5) - (
     &         (qir*normal(1) + tet1)*flur1 +
     &         (qir*normal(2) + tet2)*flur2 +
     &         (qir*normal(3) + tet3)*flur3 +
     &         0.5d0*VdotN*(flur4 - flur5)/cr +
     &         0.5d0*uar5*(flur4 + flur5)*cr2)
      cf(767)=(qir*normal(1) + tet1)
      cf(768)=(qir*normal(2) + tet2)
      cf(769)=(qir*normal(3) + tet3)
      cf(770)=0.5d0*(flur4 - flur5)/cr
      cf(771)=0.5d0*VdotN/cr
      cf(772)=-0.5d0*VdotN*(flur4 - flur5)/(cr*cr)
      cf(773)=0.5d0*(flur4 + flur5)*cr2
      cf(774)=0.5d0*uar5*cr2
      cf(775)=0.5d0*uar5*(flur4 + flur5)
      cf(776)=(flur1*cf(201)+cf(767)*cf(383)+flur2*normal(2)*cf(190)
     &   +flur1*normal(1)*cf(190)+cf(768)*cf(450)+flur3*normal(3)
     &   *cf(190)+flur3*cf(219)+cf(769)*cf(500)+cf(770)*cf(182)
     &   +cf(771)*(cf(560)-cf(597))+cf(772)*cf(238)+cf(773)*cf(169)
     &   +cf(774)*(cf(597)+cf(560))+cf(775)*cf(250)+flur2*cf(210))
      cf(777)=(cf(767)*cf(384)+flur2*normal(2)*cf(192)+flur2*cf(211)
     &   +cf(768)*cf(451)+flur3*normal(3)*cf(192)+flur3*cf(220)
     &   +cf(769)*cf(501)+cf(770)*cf(184)+cf(771)*(cf(561)-cf(598))
     &   +cf(772)*cf(239)+cf(773)*cf(171)+cf(774)*(cf(598)+cf(561))
     &   +cf(775)*cf(251)+flur1*normal(1)*cf(192))
      cf(778)=(flur1*cf(202)+cf(767)*cf(385)+flur2*normal(2)*cf(193)
     &   +cf(768)*cf(452)+flur3*normal(3)*cf(193)+flur3*cf(221)+flur1
     &   *normal(1)*cf(193)+cf(769)*cf(502)+cf(770)*cf(185)+cf(772)
     &   *cf(240)+cf(773)*cf(172)+cf(774)*(cf(599)+cf(562))+cf(775)
     &   *cf(252)+cf(771)*(cf(562)-cf(599)))
      cf(779)=(flur1*cf(203)+cf(767)*cf(386)+flur2*normal(2)*cf(194)
     &   +flur2*cf(212)+cf(768)*cf(453)+flur3*normal(3)*cf(194)+cf(769)
     &   *cf(503)+cf(770)*cf(186)+cf(771)*(cf(563)-cf(600))+cf(772)
     &   *cf(241)+cf(773)*cf(173)+cf(774)*(cf(600)+cf(563))+cf(775)
     &   *cf(253)+flur1*normal(1)*cf(194))
      cf(780)=(cf(768)*cf(454)+cf(769)*cf(504)+cf(771)*(cf(564)-cf(601))
     &   +cf(772)*cf(242)+cf(773)*cf(170)+cf(774)*(cf(601)+cf(564))
     &   +cf(775)*cf(254)+cf(767)*cf(387))
      cf(781)=(flur1*normal(1)*cf(191)+flur1*cf(204)+cf(767)*cf(388)
     &   +flur2*normal(2)*cf(191)+flur2*cf(213)+cf(768)*cf(455)+flur3
     &   *normal(3)*cf(191)+flur3*cf(222)+cf(769)*cf(505)+cf(770)
     &   *cf(183)+cf(771)*(cf(565)-cf(602))+cf(772)*cf(243)+cf(773)
     &   *cf(174)+cf(774)*(cf(602)+cf(565))+cf(775)*cf(255))
      cf(782)=(cf(774)*(cf(603)+cf(566))+cf(773)*cf(176)+cf(772)
     &   *cf(244)+cf(771)*(cf(566)-cf(603))+cf(770)*cf(187)+cf(769)
     &   *cf(506)+flur3*cf(223)+flur3*normal(3)*cf(195)+cf(768)*cf(456)
     &   +flur2*cf(214)+flur2*normal(2)*cf(195)+cf(767)*cf(389)
     &   +flur1*normal(1)*cf(195)+cf(775)*cf(256))
      cf(783)=(cf(774)*(cf(604)+cf(567))+cf(773)*cf(177)+cf(772)
     &   *cf(245)+cf(771)*(cf(567)-cf(604))+cf(770)*cf(188)+cf(769)
     &   *cf(507)+flur3*cf(224)+flur3*normal(3)*cf(196)+cf(768)*cf(457)
     &   +flur2*normal(2)*cf(196)+cf(767)*cf(390)+flur1*cf(205)+flur1
     &   *normal(1)*cf(196)+cf(775)*cf(257))
      cf(784)=(cf(774)*(cf(605)+cf(568))+cf(773)*cf(178)+cf(772)
     &   *cf(246)+cf(771)*(cf(568)-cf(605))+cf(770)*cf(189)+cf(769)
     &   *cf(508)+flur3*normal(3)*cf(197)+cf(768)*cf(458)+flur2*cf(215)
     &   +flur2*normal(2)*cf(197)+cf(767)*cf(391)+flur1*cf(206)+flur1
     &   *normal(1)*cf(197)+cf(775)*cf(258))
      cf(785)=(cf(775)*cf(259)+cf(774)*(cf(606)+cf(569))+cf(773)*cf(175)
     &   +cf(772)*cf(247)+cf(771)*(cf(569)-cf(606))+cf(769)*cf(509)
     &   +cf(768)*cf(459)+cf(767)*cf(392))
      cf(786)=(flur1*cf(198)+cf(767)*cf(379)+flur2*qir*cf(11)+flur2
     &   *cf(207)+cf(768)*cf(447)+flur3*qir*cf(14)+flur3*cf(216)
     &   +cf(769)*cf(497)+cf(770)*cf(179)+cf(771)*(cf(556)-cf(593))
     &   +cf(774)*(cf(593)+cf(556))+flur1*qir*cf(8))
      cf(787)=(flur1*cf(199)+cf(767)*cf(380)+flur2*qir*cf(12)+flur2
     &   *cf(208)+cf(768)*cf(448)+flur3*qir*cf(15)+flur3*cf(217)
     &   +cf(769)*cf(498)+cf(770)*cf(180)+flur1*qir*cf(9)+cf(774)
     &   *(cf(594)+cf(557))+ cf(771)*(cf(557)-cf(594)))
      cf(788)=(flur1*qir*cf(10)+flur1*cf(200)+flur2*qir*cf(13)+flur2
     &   *cf(209)+cf(768)*cf(449)+flur3*qir*cf(16)+flur3*cf(218)
     &   +cf(769)*cf(499)+cf(770)*cf(181)+cf(771)*(cf(558)-cf(595))
     &   +cf(774)*(cf(595)+cf(558))+cf(767)*cf(381))
      cf(789)=(cf(768)*cf(446)+cf(769)*cf(510)+cf(771)*(cf(559)-cf(596))
     &   +cf(773)*cf(168)+cf(774)*(cf(596)+cf(559))+cf(775)*cf(249)
     &   +cf(767)*cf(382))
      cf(790)=(cf(768)*cf(465)+cf(769)*cf(515)+cf(771)*(cf(574)-cf(611))
     &   +cf(774)*(cf(611)+cf(574))+cf(767)*cf(393))
      cf(791)=(cf(768)*cf(466)+cf(769)*cf(516)+cf(771)*(cf(575)-cf(612))
     &   +cf(774)*(cf(612)+cf(575))+cf(767)*cf(394))
      cf(792)=(cf(768)*cf(467)+cf(769)*cf(517)+cf(771)*(cf(576)-cf(613))
     &   +cf(774)*(cf(613)+cf(576))+cf(767)*cf(395))
      cf(793)=(cf(768)*cf(468)+cf(769)*cf(518)+cf(774)*(cf(614)+cf(577))
     &   +cf(771)*(cf(577)-cf(614))+cf(767)*cf(396))
      cf(794)=(cf(768)*cf(469)+cf(769)*cf(520)+cf(771)*(cf(579)-cf(616))
     &   +cf(774)*(cf(579)+cf(616))+cf(767)*cf(397))
      cf(795)=(cf(768)*cf(460)+cf(769)*cf(511)+cf(771)*(cf(570)-cf(607))
     &   +cf(774)*(cf(607)+cf(570))+cf(767)*cf(398))
      cf(796)=(cf(768)*cf(461)+cf(769)*cf(512)+cf(771)*(cf(571)-cf(608))
     &   +cf(774)*(cf(608)+cf(571))+cf(767)*cf(399))
      cf(797)=(cf(768)*cf(462)+cf(769)*cf(513)+cf(771)*(cf(572)-cf(609))
     &   +cf(774)*(cf(609)+cf(572))+cf(767)*cf(400))
      cf(798)=(cf(768)*cf(463)+cf(769)*cf(514)+cf(771)*(cf(573)-cf(610))
     &   +cf(774)*(cf(610)+cf(573))+cf(767)*cf(401))
      cf(799)=(cf(767)*cf(402)+cf(768)*cf(464)+cf(769)*cf(519)+cf(771)
     &   *(cf(578)-cf(615))+cf(774)*(cf(615)+cf(578)))
c
      cf(800)=0.5*rnorm
      cf(801)=phi(1)*0.5
      cf(802)=-cf(800)*cf(619)
      cf(803)=-cf(800)*0.5d0*cr2*cf(521)
      cf(804)=-cf(800)*0.5d0*cr2*cf(580)
      cf(805)=(cf(800)*(cf(86)-cf(618))+cf(801)*cf(2)) 
      cf(806)=(cf(801)*cf(3)+cf(800)*(cf(87)-cf(620)))
      cf(807)=(cf(801)*cf(4)+cf(800)*(cf(88)-cf(621)))
      cf(808)=cf(800)*cf(89)
      cf(809)=-cf(800)*cf(622)
      cf(810)=cf(800)*(cf(43)-cf(623))
      cf(811)=cf(800)*(cf(44)-cf(624))
      cf(812)=cf(800)*(cf(45)-cf(625))
      cf(813)=cf(800)*(cf(46)-cf(626))
      cf(814)=-cf(800)*cf(627)
      cf(815)=cf(800)*(cf(93)-cf(628))
      cf(816)=cf(800)*(cf(90)-cf(629))
      cf(817)=cf(800)*(cf(91)-cf(630))
      cf(818)=cf(800)*(cf(92)-cf(631))
      cf(819)=-cf(800)*cf(632)
      cf(820)=- cf(800)*cf(633)
      cf(821)=- cf(800)*cf(634)
      cf(822)=- cf(800)*cf(635)
      cf(823)=- cf(800)*cf(636)
      cf(824)=- cf(800)*cf(637)
      cf(825)=- cf(800)*cf(638)
      cf(826)=- cf(800)*cf(639)
      cf(827)=- cf(800)*cf(640)
      cf(828)=- cf(800)*cf(641)
      cf(829)=- cf(800)*cf(642)
      dphi(1) = cf(802)*dvp1
     &         + cf(803)*dvp4
     &         + cf(804)*dvp5
     &         + cf(805)*denormal(1)
     &         + cf(806)*denormal(2)
     &         + cf(807)*denormal(3)
     &         + cf(808)*devitno
     &         + cf(809)*dpstiff
     &         + cf(810)*dUg(1)
     &         + cf(811)*dUg(2)
     &         + cf(812)*dUg(3)
     &         + cf(813)*dUg(4)
     &         + cf(814)*dUg(5)
     &         + cf(825)*dUgr(1)
     &         + cf(826)*dUgr(2)
     &         + cf(827)*dUgr(3)
     &         + cf(828)*dUgr(4)
     &         + cf(829)*dUgr(5)
     &         + cf(815)*dUd(1)
     &         + cf(816)*dUd(2)
     &         + cf(817)*dUd(3)
     &         + cf(818)*dUd(4)
     &         + cf(819)*dUd(5)
     &         + cf(820)*dUdr(1)
     &         + cf(821)*dUdr(2)
     &         + cf(822)*dUdr(3)
     &         + cf(823)*dUdr(4)
     &         + cf(824)*dUdr(5)
      phi(1) = phi(1)*0.5*rnorm
      cf(830)=0.5*rnorm
      cf(831)=phi(2)*0.5
      dphi(2) = -cf(830)*cf(666)*dvp1
     &         - cf(830)*cf(678)*dvp4
     &         - cf(830)*cf(679)*dvp5
     &         + (cf(830)*cf(680)+cf(831)*cf(2))*denormal(1)
     &         + (cf(831)*cf(3)+cf(830)*cf(681))*denormal(2)
     &         + (cf(831)*cf(4)+cf(830)*cf(682))*denormal(3)
     %         + cf(830)*cf(97)*devitno
     &         - cf(830)*cf(667)*dpstiff
     %         + cf(830)*cf(683)*dUg(1)
     %         + cf(830)*cf(684)*dUg(2)
     %         + cf(830)*cf(685)*dUg(3)
     %         + cf(830)*cf(686)*dUg(4)
     %         + cf(830)*cf(687)*dUg(5)
     &         - cf(830)*cf(673)*dUgr(1)
     &         - cf(830)*cf(674)*dUgr(2)
     &         - cf(830)*cf(675)*dUgr(3)
     &         - cf(830)*cf(676)*dUgr(4)
     &         - cf(830)*cf(677)*dUgr(5)
     &         + cf(830)*cf(688)*dUd(1)
     &         + cf(830)*cf(689)*dUd(2)
     &         + cf(830)*cf(690)*dUd(3)
     &         + cf(830)*cf(691)*dUd(4)
     &         + cf(830)*cf(692)*dUd(5)
     &         - cf(830)*cf(668)*dUdr(1)
     &         - cf(830)*cf(669)*dUdr(2)
     &         - cf(830)*cf(670)*dUdr(3)
     &         - cf(830)*cf(671)*dUdr(4)
     &         - cf(830)*cf(672)*dUdr(5)
      phi(2) = phi(2)*0.5*rnorm
      cf(831)=(0.5*rnorm*(cf(103)-cf(703))+phi(3)*0.5*cf(2))
      cf(832)=(phi(3)*0.5*cf(3)+0.5*rnorm*(cf(104)-cf(704)))
      cf(833)=(phi(3)*0.5*cf(4)+0.5*rnorm*(cf(105)-cf(705)))
      dphi(3) =-0.5*rnorm*cf(706)*dvp1
     &         - 0.5*rnorm*(cf(697)*cf(521)+cf(701)*cf(521))*dvp4
     &         + 0.5*rnorm*(cf(697)*cf(580)-cf(701)*cf(580))*dvp5
     &         + cf(831)*denormal(1)
     &         + cf(832)*denormal(2)
     &         + cf(833)*denormal(3)
     &         + 0.5*rnorm*cf(106)*devitno
     &         - 0.5*rnorm*cf(707)*dpstiff
     &         + 0.5*rnorm*(cf(107)-cf(708))*dUg(1)
     &         + 0.5*rnorm*(cf(108)-cf(709))*dUg(2)
     &         + 0.5*rnorm*(cf(109)-cf(710))*dUg(3)
     &         + 0.5*rnorm*(cf(110)-cf(711))*dUg(4)
     &         + 0.5*rnorm*(normal(2)-cf(712))*dUg(5)
     &         - 0.5*rnorm*cf(723)*dUgr(1)
     &         - 0.5*rnorm*cf(724)*dUgr(2)
     &         - 0.5*rnorm*cf(725)*dUgr(3)
     &         - 0.5*rnorm*cf(726)*dUgr(4)
     &         - 0.5*rnorm*cf(727)*dUgr(5)
     &         + 0.5*rnorm*(cf(111)-cf(713))*dUd(1)
     &         + 0.5*rnorm*(cf(112)-cf(714))*dUd(2)
     &         + 0.5*rnorm*(cf(113)-cf(715))*dUd(3)
     &         + 0.5*rnorm*(cf(114)-cf(716))*dUd(4)
     &         + 0.5*rnorm*(normal(2)-cf(717))*dUd(5)
     &         - 0.5*rnorm*cf(718)*dUdr(1)
     &         - 0.5*rnorm*cf(719)*dUdr(2)
     &         - 0.5*rnorm*cf(720)*dUdr(3)
     &         - 0.5*rnorm*cf(721)*dUdr(4)
     &         - 0.5*rnorm*cf(722)*dUdr(5)
      phi(3) = phi(3)*0.5*rnorm
      dphi(4) = -0.5*rnorm*cf(761)*dvp1
     &        + 0.5*rnorm*cf(762)*dvp4
     &         + 0.5*rnorm*cf(763)*dvp5
     &      + (0.5*rnorm*cf(764)+phi(4)*0.5*cf(2))*denormal(1)
     &      + (0.5*rnorm*cf(765)+phi(4)*0.5*cf(3))*denormal(2)
     &      + (0.5*rnorm*cf(766)+phi(4)*0.5*cf(4))*denormal(3)
     &         + 0.5*rnorm*cf(118)*devitno
     &         - 0.5*rnorm*cf(747)*dpstiff
     &         + 0.5*rnorm*(cf(60)-cf(737))*dUg(1)
     &         + 0.5*rnorm*(cf(61)-cf(738))*dUg(2)
     &         + 0.5*rnorm*(cf(62)-cf(739))*dUg(3)
     &         + 0.5*rnorm*(cf(63)-cf(740))*dUg(4)
     &         + 0.5*rnorm*(normal(3)-cf(741))*dUg(5)
     &         - 0.5*rnorm*cf(753)*dUgr(1)
     &         - 0.5*rnorm*cf(754)*dUgr(2)
     &         - 0.5*rnorm*cf(755)*dUgr(3)
     &         - 0.5*rnorm*cf(756)*dUgr(4)
     &         - 0.5*rnorm*cf(757)*dUgr(5) 
     &         + 0.5*rnorm*(cf(119)-cf(742))*dUd(1)
     &         + 0.5*rnorm*(cf(120)-cf(743))*dUd(2)
     &         + 0.5*rnorm*(cf(121)-cf(744))*dUd(3)
     &         + 0.5*rnorm*(cf(122)-cf(745))*dUd(4)
     &         + 0.5*rnorm*(normal(3)-cf(746))*dUd(5)
     &         - 0.5*rnorm*cf(748)*dUdr(1)
     &         - 0.5*rnorm*cf(749)*dUdr(2)
     &      - 0.5*rnorm*cf(750)*dUdr(3)
     &      - 0.5*rnorm*cf(751)*dUdr(4)
     &      - 0.5*rnorm*cf(752)*dUdr(5)
      phi(4) = phi(4)*0.5*rnorm
      cf(834)=(-0.5*rnorm*cf(768)*cf(403)-0.5*rnorm*cf(769)*cf(470)
     &     -0.5*rnorm*cf(767)*cf(315))
      cf(835)=(-0.5*rnorm*cf(771)*cf(521)-0.5*rnorm*cf(774)*cf(521))
      cf(836)=(0.5*rnorm*cf(771)*cf(580)-0.5*rnorm*cf(774)*cf(580))
      cf(837)=(0.5*rnorm*(cf(127)-cf(786))+phi(5)*0.5*cf(2))
      cf(838)=(0.5*rnorm*(cf(128)-cf(787))+phi(5)*0.5*cf(3))
      cf(839)=(phi(5)*0.5*cf(4)+0.5*rnorm*(cf(129)-cf(788)))
      dphi(5) = cf(834)*dvp1
     &        + cf(835)*dvp4
     &        + cf(836)*dvp5
     &        + cf(837)*denormal(1)
     &        + cf(838)*denormal(2)
     &        + cf(839)*denormal(3)
     &        + 0.5*rnorm*cf(130)*devitno
     &        + 0.5*rnorm*(cf(131)-cf(789))*dpstiff
     &        + 0.5*rnorm*(cf(69)-cf(776))*dUg(1)
     &        + 0.5*rnorm*(cf(70)-cf(777))*dUg(2)
     &        + 0.5*rnorm*(cf(71)-cf(778))*dUg(3)
     &        + 0.5*rnorm*(cf(72)-cf(779))*dUg(4)
     &        + 0.5*rnorm*(cf(132)-cf(780))*dUg(5)
     &      - 0.5*rnorm*cf(795)*dUgr(1)
     &      - 0.5*rnorm*cf(796)*dUgr(2)
     &      - 0.5*rnorm*cf(797)*dUgr(3)
     &      - 0.5*rnorm*cf(798)*dUgr(4)
     &      - 0.5*rnorm*cf(799)*dUgr(5)
     &        + 0.5*rnorm*(cf(123)-cf(781))*dUd(1)
     &        + 0.5*rnorm*(cf(124)-cf(782))*dUd(2)
     &        + 0.5*rnorm*(cf(125)-cf(783))*dUd(3)
     &        + 0.5*rnorm*(cf(126)-cf(784))*dUd(4)
     &        + 0.5*rnorm*(cf(133)-cf(785))*dUd(5)
     &      - 0.5*rnorm*cf(790)*dUdr(1)
     &      - 0.5*rnorm*cf(791)*dUdr(2)
     &      - 0.5*rnorm*cf(792)*dUdr(3)
     &      - 0.5*rnorm*cf(793)*dUdr(4)
     &      - 0.5*rnorm*cf(794)*dUdr(5)
      phi(5) = phi(5)*0.5*rnorm

      if (type.eq.1) then
        updir = 0.5d0 + dsign(0.5d0, phi(1))
        phi(6) = phi(1) * (updir * Ugr(6) + (1.0d0 - updir) * Udr(6))
        cf(840)=(updir * Ugr(6) + (1.0d0 - updir) * Udr(6))
        dphi(6) = cf(840)*dphi(1)
     &          + phi(1)*updir*dUgr(6) 
     &          + phi(1)*(1.0d0-updir)*dUdr(6)
      else if (type.eq.2) then
        updir = 0.5d0 + dsign(0.5d0, phi(1))
        phi(6) = phi(1) * (updir * Ugr(6) + (1.0d0 - updir) * Udr(6))
        cf(841)=(updir * Ugr(6) + (1.0d0 - updir) * Udr(6))
        dphi(6) = cf(841)*dphi(1)
     &          + phi(1)*updir*dUgr(6) 
     &          + phi(1)*(1.0d0-updir)*dUdr(6)
        phi(7) = phi(1) * (updir * Ugr(7) + (1.0d0 - updir) * Udr(7))
        cf(842)=(updir * Ugr(7) + (1.0d0 - updir) * Udr(7))
        dphi(7) = cf(842)*dphi(1)
     &          + phi(1)*updir*dUgr(7) 
     &          + phi(1)*(1.0d0-updir)*dUdr(7)
      endif
c
      END
