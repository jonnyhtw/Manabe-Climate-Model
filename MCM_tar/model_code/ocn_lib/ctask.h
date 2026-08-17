c====================== include file "ctask.h" =========================
c
c     variables used to control how slabs are distributed in terms of
c     tasks among processors in "step.F"
c
c
c     inputs:  (set through namelist & blkdta.F)
c
c     a task is defined as one or more latitude rows (slabs). with the
c     ifdef "multitasking" directive turned on, each task goes to a 
c     separate processor.
c
c     ntasks = maximum number of tasks. if ntasks=jmt then each row
c              becomes a separate task. if ntasks=1 then all
c              rows become one task. if ntasks=n (n < jmt) then
c              each task is composed of approximately "jmt/n" rows.
c              ideally "ntasks" should equal the number of available 
c              processors in a dedicated environment.
c              
c     note:    when not "multitasking" then  "ntasks" = 1 by default.
c              "ntasks" can be set > 1 as long as the ifdef "diskless"
c              is not chosen.
c              when "ifdef multitasking" is chosen, try to set "ntasks"
c              so that each processor has the same work load. if 
c              it is known a priori that all slabs do not contain
c              similar work, then adjust "jsrowb" & "jerowb" in ocean.F
c              to keep the amount of work similar for all tasks.
c
c     outputs:
c
c     jsrowb  = array of starting j rows for tasks. "jsrowb" only needs
c               to be dimensioned by "nchunk" which is <= "jmt"
c     jerowb  = array of ending j rows for tasks. "jerowb" only needs
c               to be dimensioned by "nchunk" which is <= "jmt"
c     jsrow   = starting j row for the current task (=jsrowb(ntask))
c     jerow   = ending j row for the current task (=jerowb(ntask))
c     nchunk  = used as the tasking loop limit
c
#ifndef multitasking
# define task
#endif
cSGI
cc    task common /ctask/ jsrow, jerow
           common /ctask/ jsrow, jerow
      common /ctask1/ jsrowb(jmt), jerowb(jmt), ntasks, nchunk
c
