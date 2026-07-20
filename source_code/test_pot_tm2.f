      PROGRAM TEST_POT_TM2
C  =====================================================================
C  STAGE C: проверка pot-Tm2.f (VINIT/VSTAR) -- сплайн Table 5 + хвосты.
C  Эталонные числа посчитаны НЕЗАВИСИМО в Python (собственная реализация
C  PCHIP, построчно сверенная со scipy.interpolate.PchipInterpolator до
C  машинной точности -- см. pchip_check.py в чате), а не тем же самым
C  Fortran-кодом.  Точки покрывают: короткодействующую экспоненциальную
C  стенку (R<7.0), узлы таблицы (должны воспроизводиться ТОЧНО, т.к.
C  сплайн интерполирующий), промежутки между узлами, стыковку с хвостом
C  на Rrel=12.5, и чисто асимптотическую область (R>12.5, вплоть до 50).
C
C  ВСЕ ЭТАЛОННЫЕ R НИЖЕ ЗАДАНЫ В БОРАХ (как в Таблице 5 статьи) и
C  переводятся в ангстремы через тот же модуль physical_constants,
C  который использует сам pot-Tm2.f -- так тест не зависит от того,
C  какая именно версия CODATA у вас слинкована.
C  =====================================================================
      USE physical_constants, ONLY: bohr_to_Angstrom
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL ALLOK
      ALLOK=.TRUE.
C
      WRITE(6,*) '===================================================='
      WRITE(6,*) ' STAGE C: pot-Tm2 (V0,V2) vs Python-etalon'
      WRITE(6,*) '===================================================='
C
      CALL VINIT(1,RM,EPSIL)
      CALL VINIT(2,RM,EPSIL)
C
      CALL CHKR('V0(R=6.00) [korotkodeyst. stenka]',1, 6.00D0,
     1           52791.5441170232D0, ALLOK)
      CALL CHKR('V0(R=6.50) [korotkodeyst. stenka]',1, 6.50D0,
     1            5604.1088040406D0, ALLOK)
      CALL CHKR('V0(R=7.00) [uzel tablicy]',        1, 7.00D0,
     1             594.9065520400D0, ALLOK)
      CALL CHKR('V0(R=7.20) [uzel tablicy]',        1, 7.20D0,
     1             133.5403258700D0, ALLOK)
      CALL CHKR('V0(R=7.80) [uzel tablicy]',        1, 7.80D0,
     1            -618.7180666200D0, ALLOK)
      CALL CHKR('V0(R=8.50) [minimum, uzel]',       1, 8.50D0,
     1            -825.6458504000D0, ALLOK)
      CALL CHKR('V0(R=9.00) [uzel tablicy]',        1, 9.00D0,
     1            -784.3292713500D0, ALLOK)
      CALL CHKR('V0(R=9.60) [uzel tablicy]',        1, 9.60D0,
     1            -664.2169504000D0, ALLOK)
      CALL CHKR('V0(R=10.00) [uzel tablicy]',       1,10.00D0,
     1            -574.3069818500D0, ALLOK)
      CALL CHKR('V0(R=10.40) [uzel, nachalo redkoy zony]',1,10.40D0,
     1            -488.3781175300D0, ALLOK)
      CALL CHKR('V0(R=11.00) [uzel, redkaya zona]', 1,11.00D0,
     1            -375.8561402400D0, ALLOK)
      CALL CHKR('V0(R=11.70) [mezhdu uzlami 11.0/12.0]',1,11.70D0,
     1            -274.8247269440D0, ALLOK)
      CALL CHKR('V0(R=12.00) [posledniy uzel tablicy]',1,12.00D0,
     1            -238.1136071400D0, ALLOK)
      CALL CHKR('V0(R=12.50) [Rrel, styk s khvostom]',1,12.50D0,
     1             -96.2033955501D0, ALLOK)
      CALL CHKR('V0(R=13.00) [chisto asimptotika]', 1,13.00D0,
     1             -76.0309409278D0, ALLOK)
      CALL CHKR('V0(R=15.00) [asimptotika]',        1,15.00D0,
     1             -32.2183225195D0, ALLOK)
      CALL CHKR('V0(R=20.00) [asimptotika]',        1,20.00D0,
     1              -5.7341692179D0, ALLOK)
      CALL CHKR('V0(R=30.00) [asimptotika]',        1,30.00D0,
     1              -0.5034112894D0, ALLOK)
      CALL CHKR('V0(R=50.00) [daleko na khvoste]',  1,50.00D0,
     1              -0.0234871571D0, ALLOK)
C
      CALL CHKR('V2(R=6.00) [ploskoe prodolzhenie]',2, 6.00D0,
     1               2.2468153000D0, ALLOK)
      CALL CHKR('V2(R=6.50) [ploskoe prodolzhenie]',2, 6.50D0,
     1               2.2468153000D0, ALLOK)
      CALL CHKR('V2(R=7.00) [uzel tablicy]',        2, 7.00D0,
     1               2.2468153000D0, ALLOK)
      CALL CHKR('V2(R=7.20) [uzel tablicy]',        2, 7.20D0,
     1               3.0370500000D0, ALLOK)
      CALL CHKR('V2(R=7.80) [maksimum, uzel]',      2, 7.80D0,
     1               3.7584559000D0, ALLOK)
      CALL CHKR('V2(R=8.50) [uzel tablicy]',        2, 8.50D0,
     1               3.0630180000D0, ALLOK)
      CALL CHKR('V2(R=9.00) [uzel tablicy]',        2, 9.00D0,
     1               2.3390504000D0, ALLOK)
      CALL CHKR('V2(R=9.60) [uzel tablicy]',        2, 9.60D0,
     1               1.5651986000D0, ALLOK)
      CALL CHKR('V2(R=10.00) [uzel tablicy]',       2,10.00D0,
     1               1.1570413000D0, ALLOK)
      CALL CHKR('V2(R=10.40) [uzel, redkaya zona]', 2,10.40D0,
     1               0.8366927100D0, ALLOK)
      CALL CHKR('V2(R=11.00) [uzel, redkaya zona]', 2,11.00D0,
     1               0.4963526500D0, ALLOK)
      CALL CHKR('V2(R=11.70) [mezhdu uzlami]',      2,11.70D0,
     1               0.2655104797D0, ALLOK)
      CALL CHKR('V2(R=12.00) [posledniy uzel]',     2,12.00D0,
     1               0.1883851100D0, ALLOK)
      CALL CHKR('V2(R=12.50) [Rrel, styk s khvostom]',2,12.50D0,
     1               0.0453648791D0, ALLOK)
      CALL CHKR('V2(R=13.00) [asimptotika]',        2,13.00D0,
     1               0.0358525229D0, ALLOK)
      CALL CHKR('V2(R=20.00) [asimptotika]',        2,20.00D0,
     1               0.0027039575D0, ALLOK)
      CALL CHKR('V2(R=50.00) [daleko na khvoste]',  2,50.00D0,
     1               0.0000110754D0, ALLOK)
C
C  PROVERKA GLADKOSTI (S1-nepreryvnost) NA STYKE R=Rrel=12.5:
C  znacheniya sleva i sprava dolzhny sovpadat s tochnostyu do 1e-3
C  (eto ne uzel splayna, a granica interpolyaciya/khvost -- nebolshoy
C  skachok dopustim, no ne bolshe).
      CALL VSTAR(1,12.499999D0*bohr_to_Angstrom,VL)
      CALL VSTAR(1,12.500001D0*bohr_to_Angstrom,VR)
      CALL CHKSM('V0 continuity at Rrel=12.5',VL,VR,ALLOK)
      CALL VSTAR(2,12.499999D0*bohr_to_Angstrom,VL)
      CALL VSTAR(2,12.500001D0*bohr_to_Angstrom,VR)
      CALL CHKSM('V2 continuity at Rrel=12.5',VL,VR,ALLOK)
C
C  PROVERKA VSTAR1 (proizvodnaya): sravnenie analiticheskoy proizvodnoy
C  s chislennoy (central'naya raznost' po samoy VSTAR), NAPRYAMUYU v
C  angstremakh -- eto zaodno proveryaet novuyu tsepnuyu konversiyu
C  Angstrom<->bor v VSTAR1, kotoraya ranshe ne testirovalas'.
      DR=1.D-6
      DO ILAM=1,2
        RTEST=4.5D0
        CALL VSTAR(ILAM,RTEST-DR,VL)
        CALL VSTAR(ILAM,RTEST+DR,VR)
        DNUM=(VR-VL)/(2.D0*DR)
        CALL VSTAR1(ILAM,RTEST,DAN)
        ERR=ABS(DAN-DNUM)
        WRITE(6,'(2X,A,I1,A,/,6X,A,F14.6,A,F14.6,A,1PE10.2)')
     1    'VSTAR1 check, lambda=',ILAM,' at R=4.5 Angstrom',
     2    'analytic=',DAN,'  numeric=',DNUM,'  err=',ERR
        IF (ERR.GT.1.D-3*MAX(1.D0,ABS(DNUM))) THEN
          WRITE(6,*) '    *** VSTAR1 MISMATCH ***'
          ALLOK=.FALSE.
        ENDIF
      ENDDO
C
      WRITE(6,*) '===================================================='
      IF (ALLOK) THEN
        WRITE(6,*) ' STAGE C: VSE PROVERKI PROIDENY (PASS)'
      ELSE
        WRITE(6,*) ' STAGE C: EST RASHOZHDENIYA (FAIL) -- sm. vyshe'
      ENDIF
      WRITE(6,*) '===================================================='
C
      STOP
      END
C
      SUBROUTINE CHKR(NAME,ILAM,R,EXPECT,ALLOK)
      USE physical_constants, ONLY: bohr_to_Angstrom
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*(*) NAME
      LOGICAL ALLOK
      DOUBLE PRECISION R,EXPECT,GOT,ERR,RANG
      INTEGER ILAM
C  R IS SUPPLIED IN BOHR (MATCHING TABLE 5); CONVERT TO ANGSTROM, THE
C  UNIT pot-Tm2.f NOW EXPECTS AT ITS VSTAR INTERFACE.
      RANG=R*bohr_to_Angstrom
      CALL VSTAR(ILAM,RANG,GOT)
      ERR=ABS(GOT-EXPECT)
      WRITE(6,'(2X,A,/,6X,A,F16.8,A,F16.8,A,1PE10.2)')
     1  NAME,'got=',GOT,'  expected=',EXPECT,'  err=',ERR
      IF (ERR.GT.1.D-4*MAX(1.D0,ABS(EXPECT))) THEN
        WRITE(6,*) '    *** MISMATCH ***'
        ALLOK=.FALSE.
      ENDIF
      RETURN
      END
C
      SUBROUTINE CHKSM(NAME,VL,VR,ALLOK)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*(*) NAME
      LOGICAL ALLOK
      DOUBLE PRECISION VL,VR,ERR
      ERR=ABS(VL-VR)
      WRITE(6,'(2X,A,/,6X,A,F16.10,A,F16.10,A,1PE10.2)')
     1  NAME,'left=',VL,'  right=',VR,'  jump=',ERR
      IF (ERR.GT.1.D-3) THEN
        WRITE(6,*) '    *** DISCONTINUITY TOO LARGE ***'
        ALLOK=.FALSE.
      ENDIF
      RETURN
      END
