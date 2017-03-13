c     *******************************************************************
c 
      subroutine BallVertex(coor, coor0, nu, kf, expandCoef)
      IMPLICIT NONE
      real*8 coor(3,*), coor0(3,*), expandCoef;
      integer nu(4)
c
c
      real*8 lambda(4), zeta(4), beta(4)
      real*8 B(3), A(3,2), B1(2), A1(2,2)
      real*8 DetA1
c
      real*8 kf(12,12), k1(12,12), k2(12,12)
      real*8 k3(12,12), k4(12,12)
c
      real*8 rx12, ry12, rz12
      real*8 rx13, ry13, rz13
      real*8 rx14, ry14, rz14
      real*8 rx23, ry23, rz23
      real*8 rx24, ry24, rz24
      real*8 rx34, ry34, rz34
      real*8 rx0_12, ry0_12, rz0_12
      real*8 rx0_13, ry0_13, rz0_13
      real*8 rx0_14, ry0_14, rz0_14
      real*8 rx0_23, ry0_23, rz0_23
      real*8 rx0_24, ry0_24, rz0_24
      real*8 rx0_34, ry0_34, rz0_34
      real*8 L12, L13, L14, L23, L24, L34
      real*8 xp1, yp1, zp1
      real*8 xp2, yp2, zp2
      real*8 xp3, yp3, zp3
      real*8 xp4, yp4, zp4
      real*8 L, Lp1, Lp2, Lp3, Lp4, Lp1_0, Lp2_0, Lp3_0, Lp4_0
      real*8 k11, k22, k33, k44
      real*8 ipx, ipy, ipz
      real*8 iqx, iqy, iqz
      real*8 irx, iry, irz
      real*8 isx, isy, isz
      real*8 ipxipx, ipxipy, ipxipz, ipyipy, ipyipz, ipzipz
      real*8 iqxiqx, iqxiqy, iqxiqz, iqyiqy, iqyiqz, iqziqz
      real*8 irxirx, irxiry, irxirz, iryiry, iryirz, irzirz
      real*8 isxisx, isxisy, isxisz, isyisy, isyisz, iszisz
c
      real*8 avx(4), avy(4), avz(4)
      real*8 avx0(4), avy0(4), avz0(4)
c
      integer ro, col

c
c    Initialize current tetra stiffness matrix
c
c          do ro = 1, 12
c           do col = 1, 12
c	     k1(ro,col) = 0.0
c	     k2(ro,col) = 0.0
c	     k3(ro,col) = 0.0
c	     k4(ro,col) = 0.0
c	     kf(ro,col) = 0.0
c          enddo
c          enddo
c

c      Calculation of the length vectors
c
c
          rx12 =  coor(1,nu(2)) - coor(1,nu(1))
          ry12 =  coor(2,nu(2)) - coor(2,nu(1))
          rz12 =  coor(3,nu(2)) - coor(3,nu(1))
c
c
c
          rx13 =  coor(1,nu(3)) - coor(1,nu(1))
          ry13 =  coor(2,nu(3)) - coor(2,nu(1))
          rz13 =  coor(3,nu(3)) - coor(3,nu(1))
c
c
c
          rx14 =  coor(1,nu(4)) - coor(1,nu(1))
          ry14 =  coor(2,nu(4)) - coor(2,nu(1))
          rz14 =  coor(3,nu(4)) - coor(3,nu(1))
c
c
c
          rx23 =  coor(1,nu(3)) - coor(1,nu(2))
          ry23 =  coor(2,nu(3)) - coor(2,nu(2))
          rz23 =  coor(3,nu(3)) - coor(3,nu(2))
c
c
c
          rx24 =  coor(1,nu(4)) - coor(1,nu(2)) 
          ry24 =  coor(2,nu(4)) - coor(2,nu(2))
          rz24 =  coor(3,nu(4)) - coor(3,nu(2))
c
c
c
          rx34 =  coor(1,nu(4)) - coor(1,nu(3))
          ry34 =  coor(2,nu(4)) - coor(2,nu(3))
          rz34 =  coor(3,nu(4)) - coor(3,nu(3))

c      Calculation of the initial length vectors
c
c
          rx0_12 =  coor0(1,nu(2)) - coor0(1,nu(1))
          ry0_12 =  coor0(2,nu(2)) - coor0(2,nu(1))
          rz0_12 =  coor0(3,nu(2)) - coor0(3,nu(1))
c
c
c
          rx0_13 =  coor0(1,nu(3)) - coor0(1,nu(1))
          ry0_13 =  coor0(2,nu(3)) - coor0(2,nu(1))
          rz0_13 =  coor0(3,nu(3)) - coor0(3,nu(1))
c
c
c
          rx0_14 =  coor0(1,nu(4)) - coor0(1,nu(1))
          ry0_14 =  coor0(2,nu(4)) - coor0(2,nu(1))
          rz0_14 =  coor0(3,nu(4)) - coor0(3,nu(1))
c
c
c
          rx0_23 =  coor0(1,nu(3)) - coor0(1,nu(2))
          ry0_23 =  coor0(2,nu(3)) - coor0(2,nu(2))
          rz0_23 =  coor0(3,nu(3)) - coor0(3,nu(2))
c
c
c
          rx0_24 =  coor0(1,nu(4)) - coor0(1,nu(2)) 
          ry0_24 =  coor0(2,nu(4)) - coor0(2,nu(2))
          rz0_24 =  coor0(3,nu(4)) - coor0(3,nu(2))
c
c
c
          rx0_34 =  coor0(1,nu(4)) - coor0(1,nu(3))
          ry0_34 =  coor0(2,nu(4)) - coor0(2,nu(3))
          rz0_34 =  coor0(3,nu(4)) - coor0(3,nu(3))

c     Length of each side
c
c
      L12 = sqrt(rx12**2+ry12**2+rz12**2)
      L13 = sqrt(rx13**2+ry13**2+rz13**2)
      L14 = sqrt(rx14**2+ry14**2+rz14**2)
      L23 = sqrt(rx23**2+ry23**2+rz23**2)
      L24 = sqrt(rx24**2+ry24**2+rz24**2)
      L34 = sqrt(rx34**2+ry34**2+rz34**2)
      

c     Computation of the normal vectors of the 4 tetra faces
c
c     face 1-2-3
c

      avx(1)  =  ry12*rz13 - rz12*ry13
      avy(1)  = -rx12*rz13 + rz12*rx13
      avz(1)  =  rx12*ry13 - ry12*rx13
      L = sqrt(avx(1)**2+avy(1)**2+avz(1)**2)
      avx(1) = avx(1)/L
      avy(1) = avy(1)/L
      avz(1) = avz(1)/L
c
c     face 2-1-4
c
      avx(2)  =  ry14*rz12 - rz14*ry12
      avy(2)  = -rx14*rz12 + rz14*rx12
      avz(2)  =  rx14*ry12 - ry14*rx12
      L = sqrt(avx(2)**2+avy(2)**2+avz(2)**2)
      avx(2) = avx(2)/L
      avy(2) = avy(2)/L
      avz(2) = avz(2)/L
c
c     face 4-1-3
c
      avx(3)  =  ry13*rz14 - rz13*ry14
      avy(3)  = -rx13*rz14 + rz13*rx14
      avz(3)  =  rx13*ry14 - ry13*rx14
      L = sqrt(avx(3)**2+avy(3)**2+avz(3)**2)
      avx(3) = avx(3)/L
      avy(3) = avy(3)/L
      avz(3) = avz(3)/L
c
c     face 2-4-3
c
      avx(4)  =  ry24*rz23 - rz24*ry23
      avy(4)  = -rx24*rz23 + rz24*rx23
      avz(4)  =  rx24*ry23 - ry24*rx23
      L = sqrt(avx(4)**2+avy(4)**2+avz(4)**2)
      avx(4) = avx(4)/L
      avy(4) = avy(4)/L
      avz(4) = avz(4)/L
c     Computation of the original normal vectors of the 4 tetra faces
c
c     face 1-2-3
c

      avx0(1)  =  ry0_12*rz0_13 - rz0_12*ry0_13
      avy0(1)  = -rx0_12*rz0_13 + rz0_12*rx0_13
      avz0(1)  =  rx0_12*ry0_13 - ry0_12*rx0_13
      L = sqrt(avx0(1)**2+avy0(1)**2+avz0(1)**2)
      avx0(1) = avx0(1)/L
      avy0(1) = avy0(1)/L
      avz0(1) = avz0(1)/L
c
c     face 2-1-4
c
      avx0(2)  =  ry0_14*rz0_12 - rz0_14*ry0_12
      avy0(2)  = -rx0_14*rz0_12 + rz0_14*rx0_12
      avz0(2)  =  rx0_14*ry0_12 - ry0_14*rx0_12
      L = sqrt(avx0(2)**2+avy0(2)**2+avz0(2)**2)
      avx0(2) = avx0(2)/L
      avy0(2) = avy0(2)/L
      avz0(2) = avz0(2)/L
c
c     face 4-1-3
c
      avx0(3)  =  ry0_13*rz0_14 - rz0_13*ry0_14
      avy0(3)  = -rx0_13*rz0_14 + rz0_13*rx0_14
      avz0(3)  =  rx0_13*ry0_14 - ry0_13*rx0_14
      L = sqrt(avx0(3)**2+avy0(3)**2+avz0(3)**2)
      avx0(3) = avx0(3)/L
      avy0(3) = avy0(3)/L
      avz0(3) = avz0(3)/L
c
c     face 2-4-3
c
      avx0(4)  =  ry0_24*rz0_23 - rz0_24*ry0_23
      avy0(4)  = -rx0_24*rz0_23 + rz0_24*rx0_23
      avz0(4)  =  rx0_24*ry0_23 - ry0_24*rx0_23
      L = sqrt(avx0(4)**2+avy0(4)**2+avz0(4)**2)
      avx0(4) = avx0(4)/L
      avy0(4) = avy0(4)/L
      avz0(4) = avz0(4)/L
c   
c    
c    Finding the coordinate of the point where the perpendicular from the vertex
c    meets the opposing face

      xp1 = coor(1,nu(1))+(rx12*avx(4)+ry12*avy(4)+rz12*avz(4))*avx(4)
      yp1 = coor(2,nu(1))+(rx12*avx(4)+ry12*avy(4)+rz12*avz(4))*avy(4)
      zp1 = coor(3,nu(1))+(rx12*avx(4)+ry12*avy(4)+rz12*avz(4))*avz(4)

      xp2 = coor(1,nu(2))-(rx12*avx(3)+ry12*avy(3)+rz12*avz(3))*avx(3)
      yp2 = coor(2,nu(2))-(rx12*avx(3)+ry12*avy(3)+rz12*avz(3))*avy(3)
      zp2 = coor(3,nu(2))-(rx12*avx(3)+ry12*avy(3)+rz12*avz(3))*avz(3)

      xp3 = coor(1,nu(3))-(rx13*avx(2)+ry13*avy(2)+rz13*avz(2))*avx(2)
      yp3 = coor(2,nu(3))-(rx13*avx(2)+ry13*avy(2)+rz13*avz(2))*avy(2)
      zp3 = coor(3,nu(3))-(rx13*avx(2)+ry13*avy(2)+rz13*avz(2))*avz(2)

      xp4 = coor(1,nu(4))-(rx34*avx(1)+ry34*avy(1)+rz34*avz(1))*avx(1)
      yp4 = coor(2,nu(4))-(rx34*avx(1)+ry34*avy(1)+rz34*avz(1))*avy(1)
      zp4 = coor(3,nu(4))-(rx34*avx(1)+ry34*avy(1)+rz34*avz(1))*avz(1)

c
c     Computing lambda and zeta for interpolation using least squares

c     lambda for xp1,yp1,zp1

      B(1) = xp1-coor(1,nu(4))
      B(2) = yp1-coor(2,nu(4))
      B(3) = zp1-coor(3,nu(4))

      A(1,1) = coor(1,nu(2))-coor(1,nu(4))
      A(2,1) = coor(2,nu(2))-coor(2,nu(4))
      A(3,1) = coor(3,nu(2))-coor(3,nu(4))
      A(1,2) = coor(1,nu(3))-coor(1,nu(4))
      A(2,2) = coor(2,nu(3))-coor(2,nu(4))
      A(3,2) = coor(3,nu(3))-coor(3,nu(4))

      B1(1) = A(1,1)*B(1)+A(2,1)*B(2)+A(3,1)*B(3)
      B1(2) = A(1,2)*B(1)+A(2,2)*B(2)+A(3,2)*B(3)

      A1(1,1) = A(1,1)*A(1,1)+A(2,1)*A(2,1)+A(3,1)*A(3,1)
      A1(2,1) = A(1,2)*A(1,1)+A(2,2)*A(2,1)+A(3,2)*A(3,1)
      A1(2,2) = A(1,2)*A(1,2)+A(2,2)*A(2,2)+A(3,2)*A(3,2)

      DetA1 = 1.0/(A1(1,1)*A1(2,2)-A1(2,1)*A1(2,1))

      lambda(1) = DetA1*(A1(2,2)*B1(1) - A1(2,1)*B1(2))
      zeta(1) = DetA1*(A1(1,1)*B1(2) - A1(2,1)*B1(1))
      beta(1) = 1.0-lambda(1)-zeta(1)

c     lambda for xp2,yp2,zp2

      B(1) = xp2-coor(1,nu(4))
      B(2) = yp2-coor(2,nu(4))
      B(3) = zp2-coor(3,nu(4))

      A(1,1) = coor(1,nu(1))-coor(1,nu(4))
      A(2,1) = coor(2,nu(1))-coor(2,nu(4))
      A(3,1) = coor(3,nu(1))-coor(3,nu(4))
      A(1,2) = coor(1,nu(3))-coor(1,nu(4))
      A(2,2) = coor(2,nu(3))-coor(2,nu(4))
      A(3,2) = coor(3,nu(3))-coor(3,nu(4))

      B1(1) = A(1,1)*B(1)+A(2,1)*B(2)+A(3,1)*B(3)
      B1(2) = A(1,2)*B(1)+A(2,2)*B(2)+A(3,2)*B(3)

      A1(1,1) = A(1,1)*A(1,1)+A(2,1)*A(2,1)+A(3,1)*A(3,1)
      A1(2,1) = A(1,2)*A(1,1)+A(2,2)*A(2,1)+A(3,2)*A(3,1)
      A1(2,2) = A(1,2)*A(1,2)+A(2,2)*A(2,2)+A(3,2)*A(3,2)

      DetA1 = 1.0/(A1(1,1)*A1(2,2)-A1(2,1)*A1(2,1))

      lambda(2) = DetA1*(A1(2,2)*B1(1) - A1(2,1)*B1(2))
      zeta(2) = DetA1*(A1(1,1)*B1(2) - A1(2,1)*B1(1))
      beta(2) = 1.0-lambda(2)-zeta(2)

c     lambda for xp3,yp3,zp3

      B(1) = xp3-coor(1,nu(4))
      B(2) = yp3-coor(2,nu(4))
      B(3) = zp3-coor(3,nu(4))

      A(1,1) = coor(1,nu(1))-coor(1,nu(4))
      A(2,1) = coor(2,nu(1))-coor(2,nu(4))
      A(3,1) = coor(3,nu(1))-coor(3,nu(4))
      A(1,2) = coor(1,nu(2))-coor(1,nu(4))
      A(2,2) = coor(2,nu(2))-coor(2,nu(4))
      A(3,2) = coor(3,nu(2))-coor(3,nu(4))

      B1(1) = A(1,1)*B(1)+A(2,1)*B(2)+A(3,1)*B(3)
      B1(2) = A(1,2)*B(1)+A(2,2)*B(2)+A(3,2)*B(3)

      A1(1,1) = A(1,1)*A(1,1)+A(2,1)*A(2,1)+A(3,1)*A(3,1)
      A1(2,1) = A(1,2)*A(1,1)+A(2,2)*A(2,1)+A(3,2)*A(3,1)
      A1(2,2) = A(1,2)*A(1,2)+A(2,2)*A(2,2)+A(3,2)*A(3,2)

      DetA1 = 1.0/(A1(1,1)*A1(2,2)-A1(2,1)*A1(2,1))

      lambda(3) = DetA1*(A1(2,2)*B1(1) - A1(2,1)*B1(2))
      zeta(3) = DetA1*(A1(1,1)*B1(2) - A1(2,1)*B1(1))
      beta(3) = 1.0-lambda(3)-zeta(3)

c     lambda for xp4,yp4,zp4

      B(1) = xp4-coor(1,nu(1))
      B(2) = yp4-coor(2,nu(1))
      B(3) = zp4-coor(3,nu(1))

      A(1,1) = coor(1,nu(2))-coor(1,nu(1))
      A(2,1) = coor(2,nu(2))-coor(2,nu(1))
      A(3,1) = coor(3,nu(2))-coor(3,nu(1))
      A(1,2) = coor(1,nu(3))-coor(1,nu(1))
      A(2,2) = coor(2,nu(3))-coor(2,nu(1))
      A(3,2) = coor(3,nu(3))-coor(3,nu(1))

      B1(1) = A(1,1)*B(1)+A(2,1)*B(2)+A(3,1)*B(3)
      B1(2) = A(1,2)*B(1)+A(2,2)*B(2)+A(3,2)*B(3)

      A1(1,1) = A(1,1)*A(1,1)+A(2,1)*A(2,1)+A(3,1)*A(3,1)
      A1(2,1) = A(1,2)*A(1,1)+A(2,2)*A(2,1)+A(3,2)*A(3,1)
      A1(2,2) = A(1,2)*A(1,2)+A(2,2)*A(2,2)+A(3,2)*A(3,2)

      DetA1 = 1.0/(A1(1,1)*A1(2,2)-A1(2,1)*A1(2,1))

      lambda(4) = DetA1*(A1(2,2)*B1(1) - A1(2,1)*B1(2))
      zeta(4) = DetA1*(A1(1,1)*B1(2) - A1(2,1)*B1(1))
      beta(4) = 1.0-lambda(4)-zeta(4)

c     length of the vertex from the point on the face opposite to it

      Lp1 = sqrt((xp1-coor(1,nu(1)))**2 + (yp1-coor(2,nu(1)))**2 
     & + (zp1-coor(3,nu(1)))**2)   
      Lp2 = sqrt((xp2-coor(1,nu(2)))**2 + (yp2-coor(2,nu(2)))**2 
     & + (zp2-coor(3,nu(2)))**2)
      Lp3 = sqrt((xp3-coor(1,nu(3)))**2 + (yp3-coor(2,nu(3)))**2 
     & + (zp3-coor(3,nu(3)))**2)   
      Lp4 = sqrt((xp4-coor(1,nu(4)))**2 + (yp4-coor(2,nu(4)))**2 
     & + (zp4-coor(3,nu(4)))**2)   

      Lp1_0 = (coor0(1,nu(4))-coor0(1,nu(2)))*avx0(1)
     &      + (coor0(2,nu(4))-coor0(2,nu(2)))*avy0(1)
     &      + (coor0(3,nu(4))-coor0(3,nu(2)))*avz0(1)
      Lp2_0 = (coor0(1,nu(3))-coor0(1,nu(1)))*avx0(2)
     &      + (coor0(2,nu(3))-coor0(2,nu(1)))*avy0(2)
     &      + (coor0(3,nu(3))-coor0(3,nu(1)))*avz0(2)
      Lp3_0 = (coor0(1,nu(2))-coor0(1,nu(1)))*avx0(3)
     &      + (coor0(2,nu(2))-coor0(2,nu(1)))*avy0(3)
     &      + (coor0(3,nu(2))-coor0(3,nu(1)))*avz0(3)
      Lp4_0 = (coor0(1,nu(1))-coor0(1,nu(2)))*avx0(4)
     &      + (coor0(2,nu(1))-coor0(2,nu(2)))*avy0(4)
     &      + (coor0(3,nu(1))-coor0(3,nu(2)))*avz0(4)

      ipx = (xp1-coor(1,nu(1)))/Lp1 
      ipy = (yp1-coor(2,nu(1)))/Lp1
      ipz = (zp1-coor(3,nu(1)))/Lp1

      iqx = (xp2-coor(1,nu(2)))/Lp2
      iqy = (yp2-coor(2,nu(2)))/Lp2
      iqz = (zp2-coor(3,nu(2)))/Lp2

      irx = (xp3-coor(1,nu(3)))/Lp3
      iry = (yp3-coor(2,nu(3)))/Lp3
      irz = (zp3-coor(3,nu(3)))/Lp3

      isx = (xp4-coor(1,nu(4)))/Lp4
      isy = (yp4-coor(2,nu(4)))/Lp4
      isz = (zp4-coor(3,nu(4)))/Lp4

      k11 = (1.0/Lp1)+expandCoef*Lp1/Lp1_0**2
      k22 = (1.0/Lp2)+expandCoef*Lp2/Lp2_0**2
      k33 = (1.0/Lp3)+expandCoef*Lp3/Lp3_0**2
      k44 = (1.0/Lp4)+expandCoef*Lp4/Lp4_0**2

      ipxipx = k11*ipx*ipx
      ipxipy = k11*ipx*ipy
      ipxipz = k11*ipx*ipz
      ipyipy = k11*ipy*ipy
      ipyipz = k11*ipy*ipz
      ipzipz = k11*ipz*ipz
   
      iqxiqx = k22*iqx*iqx
      iqxiqy = k22*iqx*iqy
      iqxiqz = k22*iqx*iqz
      iqyiqy = k22*iqy*iqy
      iqyiqz = k22*iqy*iqz
      iqziqz = k22*iqz*iqz

      irxirx = k33*irx*irx
      irxiry = k33*irx*iry
      irxirz = k33*irx*irz
      iryiry = k33*iry*iry
      iryirz = k33*iry*irz
      irzirz = k33*irz*irz

      isxisx = k44*isx*isx
      isxisy = k44*isx*isy
      isxisz = k44*isx*isz
      isyisy = k44*isy*isy
      isyisz = k44*isy*isz
      iszisz = k44*isz*isz


c    stiffness matrix for the face springs

c    for spring system 1-p

      k1(1,1) = ipxipx 
      k1(1,2) = ipxipy 
      k1(1,3) = ipxipz 
      k1(1,4) = -lambda(1)*ipxipx 
      k1(1,5) = -lambda(1)*ipxipy 
      k1(1,6) = -lambda(1)*ipxipz 
      k1(1,7) = -zeta(1)*ipxipx 
      k1(1,8) = -zeta(1)*ipxipy 
      k1(1,9) = -zeta(1)*ipxipz 
      k1(1,10) = -beta(1)*ipxipx 
      k1(1,11) = -beta(1)*ipxipy 
      k1(1,12) = -beta(1)*ipxipz 
      k1(2,2) = ipyipy 
      k1(2,3) = ipyipz 
      k1(2,4) = k1(1,5)
      k1(2,5) = -lambda(1)*ipyipy 
      k1(2,6) = -lambda(1)*ipyipz 
      k1(2,7) = k1(1,8)
      k1(2,8) = -zeta(1)*ipyipy 
      k1(2,9) = -zeta(1)*ipyipz 
      k1(2,10) = k1(1,11)
      k1(2,11) = -beta(1)*ipyipy 
      k1(2,12) = -beta(1)*ipyipz 
      k1(3,3) = ipzipz 
      k1(3,4) = k1(1,6)
      k1(3,5) = k1(2,6)
      k1(3,6) = -lambda(1)*ipzipz 
      k1(3,7) = k1(1,9)
      k1(3,8) = k1(2,9)
      k1(3,9) = -zeta(1)*ipzipz 
      k1(3,10) = k1(1,12)
      k1(3,11) = k1(2,12)
      k1(3,12) = -beta(1)*ipzipz 
      k1(4,4) = lambda(1)**2 * ipxipx 
      k1(4,5) = lambda(1)**2 * ipxipy 
      k1(4,6) = lambda(1)**2 * ipxipz 
      k1(4,7) = lambda(1)*zeta(1) * ipxipx  
      k1(4,8) = lambda(1)*zeta(1) * ipxipy  
      k1(4,9) = lambda(1)*zeta(1) * ipxipz  
      k1(4,10) = lambda(1)*beta(1) * ipxipx  
      k1(4,11) = lambda(1)*beta(1) * ipxipy  
      k1(4,12) = lambda(1)*beta(1) * ipxipz  
      k1(5,5) = lambda(1)**2 * ipyipy 
      k1(5,6) = lambda(1)**2 * ipyipz 
      k1(5,7) = k1(4,8)
      k1(5,8) = lambda(1)*zeta(1)*ipyipy
      k1(5,9) = lambda(1)*zeta(1)*ipyipz
      k1(5,10) = k1(4,11)
      k1(5,11) = lambda(1)*beta(1)*ipyipy 
      k1(5,12) = lambda(1)*beta(1)*ipyipz 
      k1(6,6) = lambda(1)**2 * ipzipz 
      k1(6,7) = k1(4,9)
      k1(6,8) = k1(5,9)
      k1(6,9) = lambda(1)*zeta(1)*ipzipz
      k1(6,10) = k1(4,12)
      k1(6,11) = k1(5,12)
      k1(6,12) = lambda(1)*beta(1)*ipzipz
      k1(7,7) = zeta(1)**2 * ipxipx 
      k1(7,8) = zeta(1)**2 * ipxipy 
      k1(7,9) = zeta(1)**2 * ipxipz
      k1(7,10) = zeta(1)*beta(1)*ipxipx
      k1(7,11) = zeta(1)*beta(1)*ipxipy
      k1(7,12) = zeta(1)*beta(1)*ipxipz
      k1(8,8) = zeta(1)**2 * ipyipy
      k1(8,9) = zeta(1)**2 * ipyipz
      k1(8,10) = k1(7,11)
      k1(8,11) = zeta(1)*beta(1)*ipyipy 
      k1(8,12) = zeta(1)*beta(1)*ipyipz 
      k1(9,9) = zeta(1)**2 * ipzipz 
      k1(9,10) = k1(7,12)
      k1(9,11) = k1(8,12)
      k1(9,12) = zeta(1)*beta(1)*ipzipz
      k1(10,10) = beta(1)**2 * ipxipx 
      k1(10,11) = beta(1)**2 * ipxipy 
      k1(10,12) = beta(1)**2 * ipxipz 
      k1(11,11) = beta(1)**2 * ipyipy 
      k1(11,12) = beta(1)**2 * ipyipz 
      k1(12,12) = beta(1)**2 * ipzipz 
    
c    for spring system 2-q

      k2(1,1) = lambda(2)**2 * iqxiqx 
      k2(1,2) = lambda(2)**2 * iqxiqy 
      k2(1,3) = lambda(2)**2 * iqxiqz 
      k2(1,4) = -lambda(2)*iqxiqx 
      k2(1,5) = -lambda(2)*iqxiqy 
      k2(1,6) = -lambda(2)*iqxiqz 
      k2(1,7) = lambda(2)*zeta(2)*iqxiqx 
      k2(1,8) = lambda(2)*zeta(2)*iqxiqy 
      k2(1,9) = lambda(2)*zeta(2)*iqxiqz 
      k2(1,10) = lambda(2)*beta(2)*iqxiqx 
      k2(1,11) = lambda(2)*beta(2)*iqxiqy 
      k2(1,12) = lambda(2)*beta(2)*iqxiqz 
      k2(2,2) = lambda(2)**2 * iqyiqy 
      k2(2,3) = lambda(2)**2 * iqyiqz 
      k2(2,4) = k2(1,5)
      k2(2,5) = -lambda(2)*iqyiqy 
      k2(2,6) = -lambda(2)*iqyiqz 
      k2(2,7) = k2(1,8)
      k2(2,8) = lambda(2)*zeta(2)*iqyiqy 
      k2(2,9) = lambda(2)*zeta(2)*iqyiqz 
      k2(2,10) = k2(1,11)
      k2(2,11) = lambda(2)*beta(2)*iqyiqy 
      k2(2,12) = lambda(2)*beta(2)*iqyiqz 
      k2(3,3) = lambda(2)**2 * iqziqz 
      k2(3,4) = k2(1,6)
      k2(3,5) = k2(2,6)
      k2(3,6) = -lambda(2)*iqziqz 
      k2(3,7) = k2(1,9)
      k2(3,8) = k2(2,9)
      k2(3,9) = lambda(2)*zeta(2)*iqziqz 
      k2(3,10) = k2(1,12)
      k2(3,11) = k2(2,12)
      k2(3,12) = lambda(2)*beta(2)*iqziqz 
      k2(4,4) = iqxiqx 
      k2(4,5) = iqxiqy 
      k2(4,6) = iqxiqz 
      k2(4,7) = -zeta(2) * iqxiqx  
      k2(4,8) = -zeta(2) * iqxiqy  
      k2(4,9) = -zeta(2) * iqxiqz  
      k2(4,10) = -beta(2) * iqxiqx  
      k2(4,11) = -beta(2) * iqxiqy  
      k2(4,12) = -beta(2) * iqxiqz  
      k2(5,5) = iqyiqy 
      k2(5,6) = iqyiqz 
      k2(5,7) = k2(4,8)
      k2(5,8) = -zeta(2)*iqyiqy
      k2(5,9) = -zeta(2)*iqyiqz
      k2(5,10) = k2(4,11)
      k2(5,11) = -beta(2)*iqyiqy 
      k2(5,12) = -beta(2)*iqyiqz 
      k2(6,6) = iqziqz 
      k2(6,7) = k2(4,9)
      k2(6,8) = k2(5,9)
      k2(6,9) = -zeta(2)*iqziqz
      k2(6,10) = k2(4,12)
      k2(6,11) = k2(5,12)
      k2(6,12) = -beta(2)*iqziqz
      k2(7,7) = zeta(2)**2 * iqxiqx 
      k2(7,8) = zeta(2)**2 * iqxiqy 
      k2(7,9) = zeta(2)**2 * iqxiqz
      k2(7,10) = zeta(2)*beta(2)*iqxiqx
      k2(7,11) = zeta(2)*beta(2)*iqxiqy
      k2(7,12) = zeta(2)*beta(2)*iqxiqz
      k2(8,8) = zeta(2)**2 * iqyiqy
      k2(8,9) = zeta(2)**2 * iqyiqz
      k2(8,10) = k2(7,11)
      k2(8,11) = zeta(2)*beta(2)*iqyiqy 
      k2(8,12) = zeta(2)*beta(2)*iqyiqz 
      k2(9,9) = zeta(2)**2 * iqziqz 
      k2(9,10) = k2(7,12)
      k2(9,11) = k2(8,12)
      k2(9,12) = zeta(2)*beta(2)*iqziqz
      k2(10,10) = beta(2)**2 * iqxiqx 
      k2(10,11) = beta(2)**2 * iqxiqy 
      k2(10,12) = beta(2)**2 * iqxiqz 
      k2(11,11) = beta(2)**2 * iqyiqy 
      k2(11,12) = beta(2)**2 * iqyiqz 
      k2(12,12) = beta(2)**2 * iqziqz 
    
    
c    for spring system 3-r

      k3(1,1) = lambda(3)**2 * irxirx 
      k3(1,2) = lambda(3)**2 * irxiry 
      k3(1,3) = lambda(3)**2 * irxirz 
      k3(1,4) = lambda(3)*zeta(3)*irxirx 
      k3(1,5) = lambda(3)*zeta(3)*irxiry 
      k3(1,6) = lambda(3)*zeta(3)*irxirz 
      k3(1,7) = -lambda(3)*irxirx 
      k3(1,8) = -lambda(3)*irxiry 
      k3(1,9) = -lambda(3)*irxirz 
      k3(1,10) = lambda(3)*beta(3)*irxirx 
      k3(1,11) = lambda(3)*beta(3)*irxiry 
      k3(1,12) = lambda(3)*beta(3)*irxirz 
      k3(2,2) = lambda(3)**2 * iryiry 
      k3(2,3) = lambda(3)**2 * iryirz 
      k3(2,4) = k3(1,5)
      k3(2,5) = lambda(3)*zeta(3)*iryiry 
      k3(2,6) = lambda(3)*zeta(3)*iryirz 
      k3(2,7) = k3(1,8)
      k3(2,8) = -lambda(3)*iryiry 
      k3(2,9) = -lambda(3)*iryirz 
      k3(2,10) = k3(1,11)
      k3(2,11) = lambda(3)*beta(3)*iryiry 
      k3(2,12) = lambda(3)*beta(3)*iryirz 
      k3(3,3) = lambda(3)**2 * irzirz 
      k3(3,4) = k3(1,6)
      k3(3,5) = k3(2,6)
      k3(3,6) = lambda(3)*zeta(3)*irzirz 
      k3(3,7) = k3(1,9)
      k3(3,8) = k3(2,9)
      k3(3,9) = -lambda(3)*irzirz 
      k3(3,10) = k3(1,12)
      k3(3,11) = k3(2,12)
      k3(3,12) = lambda(3)*beta(3)*irzirz 
      k3(4,4) = zeta(3)**2 * irxirx 
      k3(4,5) = zeta(3)**2 * irxiry 
      k3(4,6) = zeta(3)**2 * irxirz 
      k3(4,7) = -zeta(3) * irxirx  
      k3(4,8) = -zeta(3) * irxiry  
      k3(4,9) = -zeta(3) * irxirz  
      k3(4,10) = zeta(3)*beta(3) * irxirx  
      k3(4,11) = zeta(3)*beta(3) * irxiry  
      k3(4,12) = zeta(3)*beta(3) * irxirz  
      k3(5,5) = zeta(3)**2 * iryiry 
      k3(5,6) = zeta(3)**2 * iryirz 
      k3(5,7) = k3(4,8)
      k3(5,8) = -zeta(3)*iryiry
      k3(5,9) = -zeta(3)*iryirz
      k3(5,10) = k3(4,11)
      k3(5,11) = zeta(3)*beta(3)*iryiry 
      k3(5,12) = zeta(3)*beta(3)*iryirz 
      k3(6,6) = zeta(3)**2 * irzirz 
      k3(6,7) = k3(4,9)
      k3(6,8) = k3(5,9)
      k3(6,9) = -zeta(3)*irzirz
      k3(6,10) = k3(4,12)
      k3(6,11) = k3(5,12)
      k3(6,12) = zeta(3)*beta(3)*irzirz
      k3(7,7) = irxirx 
      k3(7,8) = irxiry 
      k3(7,9) = irxirz
      k3(7,10) = -beta(3)*irxirx
      k3(7,11) = -beta(3)*irxiry
      k3(7,12) = -beta(3)*irxirz
      k3(8,8) = iryiry
      k3(8,9) = iryirz
      k3(8,10) = k3(7,11)
      k3(8,11) = -beta(3)*iryiry 
      k3(8,12) = -beta(3)*iryirz 
      k3(9,9) = irzirz 
      k3(9,10) = k3(7,12)
      k3(9,11) = k3(8,12)
      k3(9,12) = -beta(3)*irzirz
      k3(10,10) = beta(3)**2 * irxirx 
      k3(10,11) = beta(3)**2 * irxiry 
      k3(10,12) = beta(3)**2 * irxirz 
      k3(11,11) = beta(3)**2 * iryiry 
      k3(11,12) = beta(3)**2 * iryirz 
      k3(12,12) = beta(3)**2 * irzirz 
    
    
c    for spring system 4-s

      k4(1,1) = beta(4)**2 * isxisx 
      k4(1,2) = beta(4)**2 * isxisy 
      k4(1,3) = beta(4)**2 * isxisz 
      k4(1,4) = beta(4)*lambda(4)*isxisx 
      k4(1,5) = beta(4)*lambda(4)*isxisy 
      k4(1,6) = beta(4)*lambda(4)*isxisz 
      k4(1,7) = beta(4)*zeta(4)*isxisx 
      k4(1,8) = beta(4)*zeta(4)*isxisy 
      k4(1,9) = beta(4)*zeta(4)*isxisz 
      k4(1,10) = -beta(4)*isxisx 
      k4(1,11) = -beta(4)*isxisy 
      k4(1,12) = -beta(4)*isxisz 
      k4(2,2) = beta(4)**2 * isyisy 
      k4(2,3) = beta(4)**2 * isyisz 
      k4(2,4) = k4(1,5)
      k4(2,5) = beta(4)*lambda(4)*isyisy 
      k4(2,6) = beta(4)*lambda(4)*isyisz 
      k4(2,7) = k4(1,8)
      k4(2,8) = beta(4)*zeta(4)*isyisy 
      k4(2,9) = beta(4)*zeta(4)*isyisz 
      k4(2,10) = k4(1,11)
      k4(2,11) = -beta(4)*isyisy 
      k4(2,12) = -beta(4)*isyisz 
      k4(3,3) = beta(4)**2 * iszisz 
      k4(3,4) = k4(1,6)
      k4(3,5) = k4(2,6)
      k4(3,6) = beta(4)*lambda(4)*iszisz 
      k4(3,7) = k4(1,9)
      k4(3,8) = k4(2,9)
      k4(3,9) = beta(4)*zeta(4)*iszisz 
      k4(3,10) = k4(1,12)
      k4(3,11) = k4(2,12)
      k4(3,12) = -beta(4)*iszisz 
      k4(4,4) = lambda(4)**2 * isxisx 
      k4(4,5) = lambda(4)**2 * isxisy 
      k4(4,6) = lambda(4)**2 * isxisz 
      k4(4,7) = lambda(4)*zeta(4) * isxisx  
      k4(4,8) = lambda(4)*zeta(4) * isxisy  
      k4(4,9) = lambda(4)*zeta(4) * isxisz  
      k4(4,10) = -lambda(4) * isxisx  
      k4(4,11) = -lambda(4) * isxisy  
      k4(4,12) = -lambda(4) * isxisz  
      k4(5,5) = lambda(4)**2 * isyisy 
      k4(5,6) = lambda(4)**2 * isyisz 
      k4(5,7) = k4(4,8)
      k4(5,8) = lambda(4)*zeta(4)*isyisy
      k4(5,9) = lambda(4)*zeta(4)*isyisz
      k4(5,10) = k4(4,11)
      k4(5,11) = -lambda(4)*isyisy 
      k4(5,12) = -lambda(4)*isyisz 
      k4(6,6) = lambda(4)**2 * iszisz 
      k4(6,7) = k4(4,9)
      k4(6,8) = k4(5,9)
      k4(6,9) = lambda(4)*zeta(4)*iszisz
      k4(6,10) = k4(4,12)
      k4(6,11) = k4(5,12)
      k4(6,12) = -lambda(4)*iszisz
      k4(7,7) = zeta(4)**2 * isxisx 
      k4(7,8) = zeta(4)**2 * isxisy 
      k4(7,9) = zeta(4)**2 * isxisz
      k4(7,10) = -zeta(4)*isxisx
      k4(7,11) = -zeta(4)*isxisy
      k4(7,12) = -zeta(4)*isxisz
      k4(8,8) = zeta(4)**2 * isyisy
      k4(8,9) = zeta(4)**2 * isyisz
      k4(8,10) = k4(7,11)
      k4(8,11) = -zeta(4)*isyisy 
      k4(8,12) = -zeta(4)*isyisz 
      k4(9,9) = zeta(4)**2 * iszisz 
      k4(9,10) = k4(7,12)
      k4(9,11) = k4(8,12)
      k4(9,12) = -zeta(4)*iszisz
      k4(10,10) = isxisx 
      k4(10,11) = isxisy 
      k4(10,12) = isxisz 
      k4(11,11) = isyisy 
      k4(11,12) = isyisz 
      k4(12,12) = iszisz 

      do ro = 1, 12
	do col = ro, 12
  	  kf(ro,col) = k1(ro,col) + k2(ro,col) +
     &                  k3(ro,col) + k4(ro,col)
        enddo
      enddo


c   Symmetry imposed

      kf(2,1) = kf(1,2) 
      kf(3,1) = kf(1,3) 
      kf(3,2) = kf(2,3)
      kf(4,1) = kf(1,4)
      kf(4,2) = kf(2,4)
      kf(4,3) = kf(3,4)
      kf(5,1) = kf(1,5) 
      kf(5,2) = kf(2,5)
      kf(5,3) = kf(3,5)
      kf(5,4) = kf(4,5)
      kf(6,1) = kf(1,6) 
      kf(6,2) = kf(2,6)
      kf(6,3) = kf(3,6) 
      kf(6,4) = kf(4,6)
      kf(6,5) = kf(5,6)
      kf(7,1) = kf(1,7) 
      kf(7,2) = kf(2,7) 
      kf(7,3) = kf(3,7) 
      kf(7,4) = kf(4,7) 
      kf(7,5) = kf(5,7) 
      kf(7,6) = kf(6,7) 
      kf(8,1) = kf(1,8) 
      kf(8,2) = kf(2,8) 
      kf(8,3) = kf(3,8) 
      kf(8,4) = kf(4,8) 
      kf(8,5) = kf(5,8) 
      kf(8,6) = kf(6,8) 
      kf(8,7) = kf(7,8) 
      kf(9,1) = kf(1,9) 
      kf(9,2) = kf(2,9) 
      kf(9,3) = kf(3,9) 
      kf(9,4) = kf(4,9) 
      kf(9,5) = kf(5,9) 
      kf(9,6) = kf(6,9) 
      kf(9,7) = kf(7,9) 
      kf(9,8) = kf(8,9) 
      kf(10,1) = kf(1,10) 
      kf(10,2) = kf(2,10)
      kf(10,3) = kf(3,10)
      kf(10,4) = kf(4,10)
      kf(10,5) = kf(5,10)
      kf(10,6) = kf(6,10)
      kf(10,7) = kf(7,10)
      kf(10,8) = kf(8,10)
      kf(10,9) = kf(9,10)
      kf(11,1) = kf(1,11) 
      kf(11,2) = kf(2,11) 
      kf(11,3) = kf(3,11) 
      kf(11,4) = kf(4,11) 
      kf(11,5) = kf(5,11) 
      kf(11,6) = kf(6,11) 
      kf(11,7) = kf(7,11) 
      kf(11,8) = kf(8,11) 
      kf(11,9) = kf(9,11) 
      kf(11,10) = kf(10,11) 
      kf(12,1) = kf(1,12) 
      kf(12,2) = kf(2,12) 
      kf(12,3) = kf(3,12) 
      kf(12,4) = kf(4,12) 
      kf(12,5) = kf(5,12) 
      kf(12,6) = kf(6,12) 
      kf(12,7) = kf(7,12)
      kf(12,8) = kf(8,12)
      kf(12,9) = kf(9,12)
      kf(12,10) = kf(10,12) 
      kf(12,11) = kf(11,12) 


      END
