      SUBROUTINE GRIDOC(rin,tout)
C
C  GRIDOC STORES R15 GRID FIELD INTO T21 GRID BY USING
C   NEAREST POINT DATA
C
      PARAMETER (IX=48, IL=40, IXX=64, ILL=34)
C
      DIMENSION rin(IX,IL),tout(IXX,ILL),JRT(34)
C
      DATA JRT/1,2,3,5,6,7,8,9,11,12,13,14,15,16,18,19,20,21,22,
     *  23,25,26,27,28,29,30,32,33,34,35,36,38,39,40/
C
      DO 2 JROW=1,ILL
       JR=JRT(JROW)
       DO 2 I=1,IXX,4
        II=.75001*FLOAT(I)+1.0
        tout(I  ,JROW) = rin(II,JR)
        tout(I+1,JROW) = rin(II+1,JR)
        tout(I+2,JROW) = .5*(rin(II+1,JR)+rin(II+2,JR))
        tout(I+3,JROW) = rin(II+2,JR)
  2   CONTINUE
c
      RETURN
      END
