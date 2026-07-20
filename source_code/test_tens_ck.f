      PROGRAM TEST_TENS_CK
C  =====================================================================
C  STAGE A: изолированная проверка чистой математики (TENS1, CKME) на
C  РЕАЛЬНЫХ THRJ/THREEJ/PARSGN движка MOLSCAT -- до того как это
C  тянется в BAS9IN/CPL9. Эталонные числа посчитаны независимо в
C  Python (tm_tensors_full.py) и уже проверены там дважды (см. чат).
C  =====================================================================
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL ALLOK
      DOUBLE PRECISION TENS1, CKME
      ALLOK=.TRUE.
C
      WRITE(6,*) '===================================================='
      WRITE(6,*) ' STAGE A: TENS1 / CKME vs Python-эталон'
      WRITE(6,*) '===================================================='
C
      CALL CHK('TENS1(K=2,Q=0,J2=7,MR2=7,MC2=7)',
     1         TENS1(2,0,7,7,7), 8.5732140997D0, ALLOK)
      CALL CHK('TENS1(K=2,Q=0,J2=7,MR2=1,MC2=1)',
     1         TENS1(2,0,7,1,1), -6.1237243570D0, ALLOK)
      CALL CHK('TENS1(K=2,Q=1,J2=7,MR2=3,MC2=1)',
     1         TENS1(2,1,7,3,1), -3.8729833462D0, ALLOK)
      CALL CHK('CKME(K=0,Q=0,LR=0,MLR=0,LC=0,MLC=0)',
     1         CKME(0,0,0,0,0,0), 1.0000000000D0, ALLOK)
      CALL CHK('CKME(K=2,Q=0,LR=2,MLR=0,LC=0,MLC=0)',
     1         CKME(2,0,2,0,0,0), 0.4472135955D0, ALLOK)
      CALL CHK('CKME(K=2,Q=1,LR=2,MLR=1,LC=0,MLC=0)',
     1         CKME(2,1,2,1,0,0), 0.4472135955D0, ALLOK)
      CALL CHK('CKME(K=2,Q=0,LR=2,MLR=1,LC=2,MLC=1)',
     1         CKME(2,0,2,1,2,1), 0.1428571429D0, ALLOK)
C
C  Правило отбора: CKME должен быть строго 0 при |lR-lC|>k
      CALL CHK('CKME(K=2,Q=0,LR=0,MLR=0,LC=0,MLC=0) [должен быть 0]',
     1         CKME(2,0,0,0,0,0), 0.0D0, ALLOK)
C
      WRITE(6,*) '===================================================='
      IF (ALLOK) THEN
        WRITE(6,*) ' STAGE A: ВСЕ ПРОВЕРКИ ПРОЙДЕНЫ (PASS)'
      ELSE
        WRITE(6,*) ' STAGE A: ЕСТЬ РАСХОЖДЕНИЯ (FAIL) -- см. выше'
      ENDIF
      WRITE(6,*) '===================================================='
C
      STOP
      END
C
      SUBROUTINE CHK(NAME,GOT,EXPECT,ALLOK)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*(*) NAME
      LOGICAL ALLOK
      DOUBLE PRECISION GOT,EXPECT,ERR
      ERR=ABS(GOT-EXPECT)
      WRITE(6,'(2X,A,/,6X,A,F16.10,A,F16.10,A,1PE10.2)')
     1  NAME,'получено=',GOT,'  ожидалось=',EXPECT,'  ошибка=',ERR
      IF (ERR.GT.1.D-6) THEN
        WRITE(6,*) '    *** РАСХОЖДЕНИЕ ***'
        ALLOK=.FALSE.
      ENDIF
      RETURN
      END
