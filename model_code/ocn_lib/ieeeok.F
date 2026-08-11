      subroutine ieeeok (data, length, bigstp)
      implicit real*8 (a-h,o-z)
c
      dimension data(length)
c
      logical bigstp
c
c ---------------------------------------------------------------------
c
c   set limits to prevent ieee conversion underflow & overflow
c   on the Cray YMP.
c
c   native Cray 64 bit values smaller the that which can be 
c   successfully represented by the IEEE conversion routines
c   (variable "smli3e" ) are reset to zero.

c   native Cray 64 bit values whose magnitudes are above that which
c   can be successfully represented after conversion to the IEEE 32-bit
c   format (variable "bigi3e") are handled one of two ways...
c   if the user sets variable "bigstp" to false, the data value is
c   reset to the maximum value that can be represented by the IEEE 
c   32-bit words.  however, if the user sets variable "bigstp" to true,
c   then the program will stop upon encountering a value that would
c   result in a overflow upon conversion to the IEEE 32-bit format.
c
c   inputs: 
c      data = a data record of type real and length "length"
c           that is to be prepared for writing out to an IEEE 
c           unformatted file on the YMP.
c      length = number of elements in data vector "data"
c           to check a single number, set    length = 1
c           to check a 2-D array, set    length = idim * jdim, etc.
c      bigstp = logical variable that determines if upon encountering
c           a data value that would result in a overflow upon conversion
c           to the IEEE 32-bit format converted to the IEEE 32-bit 
c           format the program will terminate (if true) or will reset
c           the data value to the largest allowable for the IEEE 32-bit
c           format (if false).
c
c   output: 
c      data = after processing, the data vector contains numbers that
c           are within the ranges that can be handled by the YMP's
c           IEEE conversion process without underflow or overflow.
c
c ---------------------------------------------------------------------
c
      smli3e = 1.00e-37
      bigi3e = 0.33e+39
c
      do 100 n=1,length
c
        if (abs(data(n)) .lt.  smli3e) data(n) =  0.0
c
        if (    data(n)  .gt.  bigi3e) then
          if (bigstp) then
            stop 'bigpos'
          else
            data(n) =  bigi3e
          endif
        endif
c
        if (    data(n)  .lt. -bigi3e)  then
          if (bigstp) then 
            stop 'bigneg'
          else
            data(n) = -bigi3e
          endif
        endif
c
100   continue
c
      return
      end
