      PROGRAM TEST_GOLDEN
C  =====================================================================
C  STAGE D: zolotoy test (eq. 9 stat'i) na NASTOYASHCHIKH skompilirovannykh
C  POTBLK (base9-Tm_Tm_AI.f) i VSTAR (pot-Tm2.f) -- v otlichie ot Stage A/B/C,
C  kotorye proveryali eti funktsii po otdel'nosti/poparno na proizvol'nykh
C  argumentakh, zdes' stroitsya POLNAYA blochnaya matritsa (fiksirovannyy
C  M_tot, kak eto real'no delaet CPL9) i diagonaliziruetsya cherez LAPACK
C  DSYEV -- rovno tak, kak eto by delal nastoyashchiy raschet.
C
C  Bazis: |mj1,mj2;l,ml>, BEZ yadernogo spina (MIA=MIB=0 vezde -- potentsial
C  ne deystvuet na yadernyy spin, sm. CPL9: MIAR.EQ.MIAC trebovanie), l_max=12
C  (sxodimost' proverena v Python do mashinnoy tochnosti pri etom l_max --
C  sm. golden_test2.py v chate).
C
C  Etalonnye chisla poschitany NEZAVISIMO v Python (golden_test2.py +
C  tm_tensors_full.py: analiticheskaya formula (9) stat'i + otdel'naya
C  Python-realizatsiya POTBLK-podobnoy logiki cherez couple()/ck_me(),
C  NE ispol'zuyushchaya nikakoy Fortran-kod).
C
C  M_tot=1  (doubled MTOTD=2)  -> soderzhit min. sostoyanie |Om1|=|Om2|=1/2
C  M_tot=7  (doubled MTOTD=14) -> soderzhit spin-rastyanutoe |Om1|=|Om2|=7/2
C  =====================================================================
      USE physical_constants, ONLY: bohr_to_Angstrom
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER, PARAMETER :: MAXDIM=400, LMAXV=12, J2V=7, NRPTS=6
      DOUBLE PRECISION RTEST(NRPTS), EMIN_EXP(NRPTS), EMAX_EXP(NRPTS)
      LOGICAL ALLOK
      ALLOK=.TRUE.
C
C  Testovye tochki R (a0) i eq(9)-etalony dlya Om=(1/2,1/2) i (7/2,7/2),
C  poschitany v Python (tm_tensors_full.py, V0_spl/V2_spl + formula 9):
      DATA RTEST   / 8.5382691346D0, 7.5D0,          8.0D0,
     1                9.0D0,          10.0D0,        11.0D0        /
      DATA EMIN_EXP/ -862.8768714971D0, -384.19460928D0, -772.31548007D0,
     1               -812.97667116D0,   -588.47778583D0, -381.93519386D0/
      DATA EMAX_EXP/ -774.4420748749D0, -276.99442925D0, -664.99186299D0,
     1               -744.22291161D0,   -554.46785628D0, -367.34546517D0/
C
      WRITE(6,*) '===================================================='
      WRITE(6,*) ' STAGE D: zolotoy test na nastoyashchikh POTBLK+VSTAR'
      WRITE(6,*) '===================================================='
C
      CALL VINIT(1,RM,EPSIL)
      CALL VINIT(2,RM,EPSIL)
C
      DO IR=1,NRPTS
        RANG = RTEST(IR)*bohr_to_Angstrom
        CALL VSTAR(1,RANG,V0R)
        CALL VSTAR(2,RANG,V2R)
        WRITE(6,'(/,2X,A,F12.6,A)') 'R = ',RTEST(IR),' a0'
        WRITE(6,'(6X,A,F16.8,A,F16.8)') 'V0(R)=',V0R,'   V2(R)=',V2R
C
        CALL GOLDEN_BLOCK(2,LMAXV,J2V,V0R,V2R,MAXDIM,EMIN,EMAXD,NDIM1)
        WRITE(6,'(6X,A,I4,A,F16.8,A,F16.8,A,1PE10.2)')
     1    'M_tot=1 blok, dim=',NDIM1,'  E_min=',EMIN,
     2    '  expected=',EMIN_EXP(IR),'  err=',
     3    ABS(EMIN-EMIN_EXP(IR))
        IF (ABS(EMIN-EMIN_EXP(IR)).GT.1.D-2*MAX(1.D0,ABS(EMIN_EXP(IR))))
     1  THEN
          WRITE(6,*) '    *** MISMATCH (E_min) ***'
          ALLOK=.FALSE.
        ENDIF
C
        CALL GOLDEN_BLOCK(14,LMAXV,J2V,V0R,V2R,MAXDIM,EMIND,EMAXX,NDIM7)
        WRITE(6,'(6X,A,I4,A,F16.8,A,F16.8,A,1PE10.2)')
     1    'M_tot=7 blok, dim=',NDIM7,'  E_max=',EMAXX,
     2    '  expected=',EMAX_EXP(IR),'  err=',
     3    ABS(EMAXX-EMAX_EXP(IR))
        IF (ABS(EMAXX-EMAX_EXP(IR)).GT.1.D-2*MAX(1.D0,ABS(EMAX_EXP(IR))))
     1  THEN
          WRITE(6,*) '    *** MISMATCH (E_max) ***'
          ALLOK=.FALSE.
        ENDIF
      ENDDO
C
      WRITE(6,*) '===================================================='
      IF (ALLOK) THEN
        WRITE(6,*) ' STAGE D: VSE PROVERKI PROIDENY (PASS)'
      ELSE
        WRITE(6,*) ' STAGE D: EST RASHOZHDENIYA (FAIL) -- sm. vyshe'
      ENDIF
      WRITE(6,*) '===================================================='
C
      STOP
      END
C
C  =====================================================================
      SUBROUTINE GOLDEN_BLOCK(MTOTD,LMAXV,J2V,V0R,V2R,MAXDIM,
     1                         EMIN,EMAXD,NDIM)
C  Stroit blochnuyu matritsu potentsiala (K=0,2) dlya fiksirovannogo
C  doubled M_tot=MTOTD, v bazise |MSA,MSB;L,ML> (MIA=MIB=0 vezde),
C  L=0,2,...,LMAXV, i diagonalizuet cherez LAPACK DSYEV (real symmetric).
C  Vozvrashchaet EMIN=min sobstv.znach., EMAXD=max sobstv.znach., NDIM=razmer.
C  =====================================================================
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER, PARAMETER :: MAXLOC=400
      INTEGER MTOTD,LMAXV,J2V,MAXDIM,NDIM
      DOUBLE PRECISION V0R,V2R,EMIN,EMAXD
      INTEGER MSAA(MAXLOC),MSBA(MAXLOC),LA(MAXLOC),MLA(MAXLOC)
      DOUBLE PRECISION VMAT(MAXLOC,MAXLOC),W(MAXLOC)
      DOUBLE PRECISION WORK(4*MAXLOC)
      DOUBLE PRECISION POTBLK,V00,V22
      INTEGER LWORK,INFO
C
C  MAXDIM (dummy arg) sluzhit tol'ko dlya proverki granitsy N<=MAXLOC
C  i dolzhen sovpadat' s MAXLOC (perenositsya iz glavnoy programmy radi
C  edinoobraziya, no fakticheski massivy zdes' fiksirovannogo razmera).
C
      N=0
      DO MSA=-J2V,J2V,2
        DO MSB=-J2V,J2V,2
          MLD=MTOTD-MSA-MSB
          IF (MOD(MLD,2).NE.0) CYCLE
          MLPHYS=MLD/2
          DO LL=0,LMAXV,2
            IF (ABS(MLPHYS).GT.LL) CYCLE
            N=N+1
            IF (N.GT.MAXLOC) THEN
              WRITE(6,*) '*** MAXLOC EXCEEDED IN GOLDEN_BLOCK ***'
              STOP
            ENDIF
            MSAA(N)=MSA
            MSBA(N)=MSB
            LA(N)=LL
            MLA(N)=MLPHYS
          ENDDO
        ENDDO
      ENDDO
      NDIM=N
C
      DO I=1,N
        DO J=I,N
          V00=POTBLK(0,MSAA(I),0,MSBA(I),0,LA(I),2*MLA(I),
     1               MSAA(J),0,MSBA(J),0,LA(J),2*MLA(J),J2V)
          V22=POTBLK(2,MSAA(I),0,MSBA(I),0,LA(I),2*MLA(I),
     1               MSAA(J),0,MSBA(J),0,LA(J),2*MLA(J),J2V)
          VMAT(I,J)=V0R*V00+V2R*V22
          VMAT(J,I)=VMAT(I,J)
        ENDDO
      ENDDO
C
      LWORK=4*MAXLOC
      CALL DSYEV('N','U',N,VMAT,MAXLOC,W,WORK,LWORK,INFO)
      IF (INFO.NE.0) THEN
        WRITE(6,*) '*** DSYEV FAILED, INFO=',INFO
        STOP
      ENDIF
      EMIN=W(1)
      EMAXD=W(N)
      RETURN
      END
