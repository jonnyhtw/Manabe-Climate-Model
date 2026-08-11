c=============================== start  camazon.h ======================
c
      common /criver/ kamzn, iamzns, iamzne, jamzns, jamzne
c
c     variables:
c        kamzn  = level thru which p-e+r will be spread (for example,
c                 to mix levels 1 & 2, set kamzn to 2
c        iamzns, iamzne = start and end indices (inclusive) for
c                 rectangular region to be influenced by river mixing
c        jamzns, jamzne = start and end indices (inclusive) for
c                 rectangular region to be influenced by river mixing
c        (these are set in blkdta.F)
c
c=============================== end    camazon.h ======================
