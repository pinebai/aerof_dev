      subroutine lworksizesvd(ictxt, nprow, npcol, myrow, mycol, m, n, 
     $           rowblock, colblock, lwork)
      integer ictxt, nprow, npcol, myrow, mycol, m, n, rowblock
      integer colblock, locLLD, icpu, rowindex, lwork
c     
      EXTERNAL BLACS_EXIT, BLACS_GRIDEXIT, BLACS_GRIDINFO,
     $         DESCINIT, PDGESVD, SL_INIT
c
      integer mp, nq, mp0, nq0
      integer watobd, wpdlange, wpdgebrd
      integer wpdlared1d, wpdlared2d, wbdtosvd
      integer wdbdsqr, myrowc, nru, sizeq, wpdormbrqln
      integer mycolr, ncvt, sizep, wpdormbrprt, nb
      integer sizeMin
c
      CALL SL_INIT(ictxt, nprow, npcol)
      CALL BLACS_GRIDINFO(ictxt, nprow, npcol, myrow, mycol)
c     
c     Computing the size of work
c
      nb = colblock
      mb = rowblock
      mp = numroc(m, nb, myrow, 0, nprow)
      nq = numroc(n, nb, mycol, 0, npcol)
      mp0 = numroc(m, nb, 0, 0, nprow)
      nq0 = numroc(n, nb, 0, 0, npcol)
c
      wpdlange = nq
c     wpdlange = mp
      wpdgebrd = nb * (mp + nq + 1) + nq
      wpdlared1d = nq0
      wpdlared2d = mq0
c
      sizeMin = min(m,n)
      sizeB = max(m,n)
c
      wdbdsqr = max(1, 4*sizeMin-4)
      myrowc = myrow*npcol + mycol
      nru = numroc(sizeMin, 1, myrowc, 0, nprow*npcol)
      sizeq = numroc(sizeMin, nb, mycol, 0, npcol)
      wpdormbrqln = max(nb*(nb-1)/2, (sizeq+mp)*nb)+nb*nb
      mycolr = myrow*npcol + mycol
      ncvt = numroc(sizeMin, 1, mycolr, 0, nprow*npcol)
      sizep = numroc(sizeMin, nb, myrow, 0, nprow)
      wpdormbrprt = max(mb*(mb-1)/2, (sizep+nq)*mb)+mb*mb
c
      wbdtosvd = sizeMin*(nru+ncvt) 
     $ + max(wdbdsqr, wpdormbrqln, wpdormbrprt)
c
      watobd = max(max(max(wpdlange, wpdgebrd),wpdlared2d),wpdlared1d)
c
      lwork = 2 + 6*sizeB + max(watobd, wbdtosvd)
c
c
      return
      end
c
