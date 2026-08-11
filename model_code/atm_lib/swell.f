      subroutine swell(sml,big,kn)
c
c  swell moves lower resolution spherical harmonics (in sml) into a
c   higher resolution field (in big) while zeroing out the higher
c   harmonic components missing in the lower resolution field.
c   the current parameters in this subroutine are for transferring
c   r15 harmonics to r21 fields for 9 levels.
c
      parameter (mx=16,jx=17,mb=22,jb=23,kx=9)
c
      complex sml(mx,jx,kx),big(mb,jb,kx)
c
      do 2 k=1,kn
      do 2 n=1,jb
      do 2 m=1,mb
  2   big(m,n,k)=(0.0,0.0)
      do 3 k=1,kn
      do 3 n=1,jx
      do 3 m=1,mx
  3   big(m,n,k)=sml(m,n,k)
      return
      end
