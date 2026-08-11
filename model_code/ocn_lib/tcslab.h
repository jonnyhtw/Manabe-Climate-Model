c====================== include file "tcslab.h" ========================
c
c     slab block extensions for turbulent closure    
c
#ifndef multitasking
# define task
#endif
cSGI
cc    task common /slabs/
           common /slabs/
     $                   vvc(imt,km,nslabs,ntau)
     $,                  vdc(imt,km,nslabs,ntau)
     $,                  vdqc(imt,km,nslabs,ntau)
     $,                  q2(imt,km,nslabs,ntau)
#  ifdef leq
     $,                  q2l(imt,km,nslabs,ntau)
#  endif
c
c     q2a   = buffer area to hold the updated "n+1" q2 slab data
c     vvca  = buffer area to hold the updated "n+1" km slab data
c     vdca  = buffer area to hold the updated "n+1" kh slab data
c     vdqca = buffer area to hold the updated "n+1" kq slab data
c     q2la  = buffer area to hold the updated "n+1" q2l slab data
c
cSGI
cc    task common /bufout/
           common /bufout/
     $                   vvca(imt,km)
     $,                  vdca(imt,km)
     $,                  vdqca(imt,km)
     $,                  q2a (imt,km)
#ifdef leq
     $,                  q2la(imt,km)
#endif
c
