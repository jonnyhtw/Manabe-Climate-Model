c====================== include file "sltdat.h" =========================
c
c  totslt = Total salt content computed in ocn1st
c           from tinit(1,2)
c  sltnrm = The salt normalization term which will be
c           added to the existing salt values in TRACER,
c           It is computed in step.
c  sltsmj = The computed integral of the salt content.
c           It is computed in TRACER and unloaded in STEP
c  sltsum = Total salt content computed in STEP every
c           nmix*5 timesteps.
c
c  use konblk(100) to store sltsum for use by subsequent runs
c
      common / sltdat / totslt, sltnrm, sltsmj
c
