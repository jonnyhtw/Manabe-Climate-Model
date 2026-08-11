c====================== include file "cppmix.h" ========================
c
c     variables for vertical mixing scheme of pacanowski & philander 
c
c     inputs: (initialized in include file "blkdta.h")
c
c
c     fricmx = max vertical mixing coefficient
c     wndmix = min vertical mixing in level 1 to simulate wind mixing
c
c
c     derived quantities:
c
c     drhop   = vertical density difference on row j+1  (p for j+1)
c                ("t" grid box bottoms)
c     drhoc   = vertical density difference on row j (c for central row)
c                ("t" grid box bottoms)
c     riuc   = richardson number for row j   ("u" grid box bottoms)
c     rium   = richardson number for row j-1 ("u" grid box bottoms)
c
c
#ifndef multitasking
# define task
#endif
cSGI
      task common /cppmix/ drhop(imt,km), drhoc(imt,km)
           common /cppmix/ drhop(imt,km), drhoc(imt,km)
     $,                    riuc(imt,km), rium(imt,km)
      common /cppmx1/ wndmix, fricmx
c
