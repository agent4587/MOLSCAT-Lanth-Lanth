      PROGRAM TEST_POTBLK
C  =====================================================================
C  STAGE B: проверка POTBLK (сборка TENSOR_ME + CKME, сумма по q со знаком
C  (-1)^q) -- в отличие от STAGE A (test_tens_ck.f), который проверял
C  TENS1/CKME ПО ОТДЕЛЬНОСТИ, здесь проверяется их СТЫКОВКА ровно так, как
C  её использует CPL9.  Эталонные числа посчитаны НЕЗАВИСИМЫМ методом --
C  прямым построением матрицы тензора [j(x)j]_Kq через угловые операторы
C  (numpy, kron/couple), а не через ту же замкнутую CG-формулу, что и сам
C  Fortran-код (см. verify_fortran_formulas.py, блок STAGE B).
C
C  ВАЖНО: MLR/MLC -- аргументы POTBLK -- УДВОЕНЫ (2*ml), в отличие от
C  аргументов CKME, где ml уже неудвоен.  Нуклеарно-спиновые аргументы
C  (MIAR,MIBR,MIAC,MIBC) внутри POTBLK не используются (эта проверка
C  сделана на уровне CPL9-вызова), поэтому здесь передаются нули.
C  =====================================================================
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL ALLOK
      DOUBLE PRECISION POTBLK
      ALLOK=.TRUE.
C
      WRITE(6,*) '===================================================='
      WRITE(6,*) ' STAGE B: POTBLK vs независимое матричное построение'
      WRITE(6,*) '===================================================='
C
C  POTBLK(K,MSAR,MIAR,MSBR,MIBR,LR,MLR,MSAC,MIAC,MSBC,MIBC,LC,MLC,J2)
C  MLR/MLC ниже = 2*(физический ml), см. комментарий выше.
C
      CALL CHK('POTBLK mar=7,mbr=7,l=2,ml=0 / mac=7,mbc=7,l=2,ml=0',
     1   POTBLK(2, 7,0, 7,0, 2,0,  7,0, 7,0, 2,0,  7),
     2   4.8989794856D0, ALLOK)
      CALL CHK('POTBLK mar=7,mbr=7,l=2,ml=0 / mac=5,mbc=7,l=2,ml=0',
     1   POTBLK(2, 7,0, 7,0, 2,0,  5,0, 7,0, 2,0,  7),
     2   0.0D0, ALLOK)
      CALL CHK('POTBLK mar=5,mbr=-3,l=2,ml=2 / mac=7,mbc=-5,l=0,ml=0',
     1   POTBLK(2, 5,0,-3,0, 2,4,  7,0,-5,0, 0,0,  7),
     2   0.0D0, ALLOK)
      CALL CHK('POTBLK mar=7,mbr=-7,l=2,ml=0 / mac=-7,mbc=7,l=2,ml=0',
     1   POTBLK(2, 7,0,-7,0, 2,0, -7,0, 7,0, 2,0,  7),
     2   0.0D0, ALLOK)
      CALL CHK('POTBLK mar=1,mbr=3,l=2,ml=-2 / mac=3,mbc=1,l=2,ml=2',
     1   POTBLK(2, 1,0, 3,0, 2,-4, 3,0, 1,0, 2,4,  7),
     2   0.0D0, ALLOK)
      CALL CHK('POTBLK mar=-3,mbr=5,l=2,ml=-2 / mac=-3,mbc=1,l=0,ml=0',
     1   POTBLK(2,-3,0, 5,0, 2,-4,-3,0, 1,0, 0,0,  7),
     2   3.0000000000D0, ALLOK)
      CALL CHK('POTBLK mar=-7,mbr=-1,l=4,ml=3 / mac=-5,mbc=-1,l=2,ml=2',
     1   POTBLK(2,-7,0,-1,0, 4,6, -5,0,-1,0, 2,4,  7),
     2   3.0000000000D0, ALLOK)
      CALL CHK('POTBLK mar=-1,mbr=-3,l=0,ml=0 / mac=-5,mbc=-3,l=2,ml=2',
     1   POTBLK(2,-1,0,-3,0, 0,0, -5,0,-3,0, 2,4,  7),
     2   3.0000000000D0, ALLOK)
      CALL CHK('POTBLK mar=1,mbr=7,l=4,ml=2 / mac=5,mbc=7,l=2,ml=0',
     1   POTBLK(2, 1,0, 7,0, 4,4,  5,0, 7,0, 2,0,  7),
     2   1.6598500055D0, ALLOK)
      CALL CHK('POTBLK mar=3,mbr=5,l=2,ml=1 / mac=5,mbc=5,l=0,ml=0',
     1   POTBLK(2, 3,0, 5,0, 2,2,  5,0, 5,0, 0,0,  7),
     2  -3.0983866770D0, ALLOK)
      CALL CHK('POTBLK mar=-3,mbr=1,l=2,ml=-1 / mac=-3,mbc=-3,l=4,ml=1',
     1   POTBLK(2,-3,0, 1,0, 2,-2,-3,0,-3,0, 4,2,  7),
     2   1.1065666703D0, ALLOK)
C
      WRITE(6,*) '===================================================='
      IF (ALLOK) THEN
        WRITE(6,*) ' STAGE B: VSE PROVERKI PROIDENY (PASS)'
      ELSE
        WRITE(6,*) ' STAGE B: EST RASHOZHDENIYA (FAIL) -- sm. vyshe'
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
     1  NAME,'got=',GOT,'  expected=',EXPECT,'  err=',ERR
      IF (ERR.GT.1.D-6) THEN
        WRITE(6,*) '    *** MISMATCH ***'
        ALLOK=.FALSE.
      ENDIF
      RETURN
      END
