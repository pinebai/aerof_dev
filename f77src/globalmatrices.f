      subroutine globalmatrices(nprow, npcol, m_a, n_a, n_b, 
     $           rowblock, colblock, locLLD, icpu, desc_a, desc_b)
      integer ictxt, nprow, npcol, myrow, mycol, m, n, rowblock
      integer icpu,m_a, n_a, n_b, colblock, info, locLLD 
c			===================
c			purpose: initialize process grid and descriptor vectors for global
c							matrices A, B
c			input: nprow, npcol, m_, n_, rowblock, colblock, locLLD
c			output: desc_a, desc_b
c			===================
c
c			===================
c			NOTES
c			info: indicates successful execution
c			===================
      integer desc_a(9), desc_b(9)
c     
      EXTERNAL BLACS_EXIT, BLACS_GRIDEXIT, BLACS_GRIDINFO,
     $         DESCINIT, SL_INIT
c
c			===================
c			INITIALIZE THE PROCESS GRID
c			===================
c		
c			input: nprow, npcol
c			output: ictxt
      CALL SL_INIT(ictxt, nprow, npcol)
c
c			input: ictxt
c			output: myrow, mycol (already have nprow, npcol)
      CALL BLACS_GRIDINFO(ictxt, nprow, npcol, myrow, mycol)
c
c			===================
c			INITIALIZE DESCRIPTOR VECTOR FOR A, B
c			===================
c			
c			input: m_a, n_a, rowblock, colblock, ictxt, locLLD
c			output: desc_a, desc_b, info
c
      CALL DESCINIT(desc_a, m_a, n_a, rowblock, colblock, 0, 0,
     $              ictxt, locLLD, info)
c
      CALL DESCINIT(desc_b, m_a, n_b, rowblock, colblock, 0, 0, 
     $              ictxt, locLLD, info)
c
c			===================
c			ERROR MESSAGE IF FAILURE
c			===================
      IF (info.ne.0) THEN
      WRITE( 6, FMT = 30)icpu, info 
c
30    FORMAT( 'cpu ', I3, ' failed in least squares with info: ', I3  )
c
      end if
      return
      end
c
