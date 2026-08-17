c====================== include file "temp.h" ==========================
c
c     temporary scratch areas for use where needed.
c     zzi  = temporary storage for implicit vertical diffusion
c
#ifndef multitasking
# define task
#endif
cSGI
cc    task common /temps/ 
           common /temps/ tempik(imt,km,2)
#if defined implicitvmix || defined isopycmix
     $,                   zzi(imt,km,2)
#endif
c
