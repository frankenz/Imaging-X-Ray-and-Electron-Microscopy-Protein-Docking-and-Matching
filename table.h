C******************************************
C table.h
C******************************************
c
      CHARACTER*4 codes
      CHARACTER*1 CID_omitA,RID_omitA,CID_omitB,RID_omitB
      CHARACTER*1 CID_sitA,RID_sitA,CID_sitB,RID_sitB
      INTEGER convert_tbl
      INTEGER NP_tbl, rNUM_omitA, rNUM_omitB
C
      COMMON/TABLE/codes(35)
      COMMON/TABLE/convert_tbl(35, 2)
      COMMON/NPTBL/NP_tbl(40)
C
      COMMON/CID/CID_omitA(20),CID_omitB(20),CID_sitA(5,50),
     1CID_sitB(5,50)
C
      COMMON/RID/RID_omitA(20),rNUM_omitA(20),RID_omitB(20),
     1rNUM_omitB(20),RID_sitA(5,50),RID_sitB(5,50)
C
