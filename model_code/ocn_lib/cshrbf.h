c====================== include file "cshrbf.h" =======================
c
c
c     when multitasking ... only shared variables are allowed to be
c     used outside the parallel region. "ta" & "bufsl" are shared
c     variables with the same name as their private counterparts (which
c     are in "slabs.h") . this approach leads to no increase in memory
c     when not multitasking
c
      common /cshrbf/ ta(imt,km,nvar)
      dimension bufsl(imt,km,nvar)
      equivalence (ta,bufsl)
c
