       SUBROUTINE BAS9IN(PRTP,IBOUND,IPRINT)
C  Adapted from base9-alk_alk_ucpld.f (Copyright (C) 2025 J. M. Hutson &
C  C. R. Le Sueur, GNU GPL v3) for Tm2 collisions.
C
C  =====================================================================
C  base9-Tm_Tm.f  --  PLUG-IN BASIS-SET SUITE FOR Tm + Tm (2F7/2 atoms)
C  =====================================================================
C  Physics: Tiesinga, Klos, Li, Petrov, Kotochigova, New J. Phys. 23,
C  085007 (2021).  Interaction potential (their Eq. 1):
C
C    V(R) = Sum_{k=0,2,4} Sum_i V^(i)_k(R) * Sum_q (-1)^q T^(i)_kq * C_{k,-q}(Rhat)
C
C  This version implements the i=1 family in full, plus the magnetic
C  dipole-dipole part of the i=2, k=2 term:
C     k=0, i=1 : T^(1)_00 = identity        (isotropic potential V0(R))
C     k=2, i=1 : T^(1)_2q = [j1(x)j1]_2q + [j2(x)j2]_2q (rank-2 potential V2(R))
C     k=2, i=2 : T^(2)_2q = [j1(x)j2]_2q    (magnetic dipole-dipole,
C                V^(2)_2(R) = D^(2)_2/R^3, Sec 2.3 of the paper).  The
C                C^(2)_2/R^6 piece of the same V^(2)_2(R) is at least
C                two orders of magnitude smaller than D^(2)_2/R^3 over
C                the whole range of interest (paper, Sec 2.3) and is
C                omitted.  The remaining four spin-tensors (k=0,i=2
C                spin-exchange; k=0,2,i=3; k=4,i=1) are not included.
C
C  Basis set (uncoupled, same architecture as base9-alk_alk_ucpld.f):
C     |gammaA> |gammaB> |L ML>,   gamma = |j mj>|i mi>
C  Both atoms identical (169Tm is a BOSON: 69 protons+100 neutrons+69
C  electrons = 238, even -> Bose symmetry, exchange sign +1 for even
C  L+n where n=0 here so effectively symmetric combination for even L
C  after the (msA,miA)<->(msB,miB) exchange test -- IDENTN=.TRUE. as
C  for Rb2).
C
C  IMPORTANT UNITS CONVENTION (same as parent file): all projections
C  (MSA=2*mj_A, MIA=2*mi_A, MTOT=2*Mtot, ML=2*ml) are DOUBLED integers.
C  L itself (orbital angular momentum) is NOT doubled.
C  =====================================================================
      USE efvs, ONLY: EFV, EFVNAM, EFVUNT, MAPEFV, NEFV
      USE potential, ONLY: IREF, NCONST, NEXTMS, NEXTRA, NVLBLK, VCONST
      USE basis_data, ONLY: JHALF
      USE physical_constants, ONLY: Giga_in_SI, bohr_magneton, g_e,
     c                              hartree_in_inv_cm,
     c                              inverse_fine_structure_constant,
     c                              speed_of_light_in_cm
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      CHARACTER(8) PRTP(4),QNAME(10)
      LOGICAL LEVIN,EIN,LCOUNT,IDENTN
      INTEGER M,MP,JJ
      DOUBLE PRECISION E1,E2,EPA,EPB,EMA,EMB,E,BFIELD,EREF
      DIMENSION E(900)
      DIMENSION LREQ(10),MFREQ(10)
      DIMENSION JSTATE(*),VL(*),IV(*),JSINDX(*),L(*),CENT(*),LAM(*)
      DIMENSION DGVL(*)
      DIMENSION MONQN(NQN1)
      DIMENSION NPTS(NDIM),XPT(MXPT,NDIM),XWT(MXPT,NDIM),
     1          X(MX)
C
C  NAMELIST &BASIS9 for Tm2.  ISA/INUCA are DOUBLED (2j, 2i).
C  For 169Tm: j=7/2 -> ISA=7 ; i=1/2 -> INUCA=1  (both defaults below).
C  HFSPLA: hyperfine coupling parameter. Provisionally we treat the
C  context value as the coupling constant zeta (Hhf = h*zeta*I.J) and
C  set HFSPLA = zeta*(j+1/2) in GHz, matching the parent formula's
C  requirement that ANSA reproduces zeta via ANSA=2*HFSPLA/NSFAC.
      NAMELIST /BASIS9/ ISA,ISB,GSA,GSB,INUCA,INUCB,
     1                  HFSPLA,HFSPLB,GA,GB,LMAX,NREQ,LREQ,MFREQ,
     2                  NEXTRA
C
      GS=-g_e
      ALPINV=inverse_fine_structure_constant
      GHZCM=speed_of_light_in_cm/Giga_in_SI
      AUCM=hartree_in_inv_cm
      BM=bohr_magneton
C
      PRTP(1)='TM - TM '
      PRTP(2)='2F7/2 AT'
      PRTP(3)='OMS + MA'
      PRTP(4)='G FIELD '
C
C  SET UP ELEMENTS OF efvs MODULE (one EFV: magnetic field)
      NEFV=1
      EFVNAM(1)='MAGNETIC Z FIELD'
      EFVUNT(1)='GAUSS'
      MAPEFV=2
C  NCONST=2: hyperfine (block MXLAM+1) + Zeeman (block MXLAM+2)
      NCONST=2
C
      VCONST(1)=1.D0/GHZCM
      IBOUND=0
      LMAX=0
C
C  --- Tm2 defaults (169Tm, both atoms identical) ---
      ISA=7
C  Lande g-factor for the 2F7/2 ground term, gJ ~ 1.14119 (NIST atomic spectra database)
      GSA=1.14119D0
      INUCA=1
      GA=0.D0
C  A = -374.1374 MHz (Evans, Sandars & Woodgate, PR 128, 2238 (1962))
C  HFSPLA = A*(j+1/2), j=7/2 -> factor 4
      HFSPLA = -374.1374D-3*4.D0   ! GHz, = -1.4965496 GHz
      ISB=ISA
      GSB=GSA
      INUCB=INUCA
      GB=GA
      HFSPLB=HFSPLA
      IDENTN=.TRUE.
      JHALF=0
      NREQ=0
      DO I=1,10
        LREQ(I)=-1
        MFREQ(I)=999
      ENDDO
C
      NEXTRA=0
      NEXTMS(1)=1
      NEXTMS(2)=1
C
      READ(5,BASIS9)
C
      IF (NEXTRA.GT.2) THEN
        WRITE(6,*) ' *** WARNING: THERE ARE ONLY 2 EXTRA '//
     1             'OPERATORS CODED IN CPL9'
        NEXTRA=2
      ENDIF
C
C  NOTE: unlike the alkali case, there is no "singlet/triplet" total-
C  spin structure here (the potential is expanded directly in tensor
C  ranks k of the ELECTRONIC angular momenta j1,j2, not in the total
C  electron spin).  MXLAM (rung 1-3) = 2 : k=0 (isotropic) and k=2
C  (rank-2).  Set in POTIN9/&POTL, must match this value.
C
      IFMAX=INUCA+ISA
      IFMIN=ABS(INUCA-ISA)
      NSFAC=(IFMAX*(IFMAX+2)-IFMIN*(IFMIN+2))/4
      ANSA=0.D0
      IF (NSFAC.NE.0) ANSA=2.D0*HFSPLA/DBLE(NSFAC)
      ANSB=ANSA
C
      IF (IPRINT.LE.0) RETURN
C
      WRITE(6,'(2X,4A8/)') PRTP
      WRITE(6,601) ISA,GSA,INUCA,GA,HFSPLA,ANSA
  601 FORMAT(/'  Tm ATOM WITH J =',I2,'/2,     MU_J   =',F12.6,
     1       ' MU_B',/,14X,'I =',I2,'/2,     MU_NUC =',F12.6,
     2       ' MU_B',/,14X,'HYPERFINE SPLITTING =',F12.6,
     3       ', COUPLING CONST =',F12.6,' GHZ',/)
      WRITE(6,'(2X,A,I4)') 'L UP TO ',LMAX
      IF (NREQ.GT.0) THEN
        WRITE(6,*) ' ONLY THE FOLLOWING L,MF PAIRS ARE INCLUDED',
     1             ' (NOTE MFREQ >= 999 INCLUDES ALL MF)'
        DO IREQ=1,NREQ
        WRITE(6,'(2X,A,I4,A,I5)') ' L = ',LREQ(IREQ),
     1                            ', MF = ',MFREQ(IREQ),'/2'
        ENDDO
      ENDIF
      WRITE(6,*) ' TENSOR TERMS INCLUDED: k=0,i=1 (isotropic), ',
     1           'k=2,i=1 (rank-2, Table 5), k=2,i=2 (magnetic ',
     2           'dipole-dipole, D(2)_2/R^3 only).'
C
      RETURN
C========================================================== END OF BAS9IN
C
      ENTRY SET9(LEVIN,EIN,NSTATE,JSTATE,NQN,QNAME,MXPAR,NLABV,IPRINT)
C
C  UNCHANGED FROM base9-alk_alk_ucpld.f -- this loop is already fully
C  generic in ISA/INUCA/ISB/INUCB, so it works unmodified for j=7/2.
C
      MXPAR=2
      NLABV=1
      NQN=5
      QNAME(1)='  2*MJA '
      QNAME(2)='  2*MIA '
      QNAME(3)='  2*MJB '
      QNAME(4)='  2*MIB '
C
      DO ILOOP=1,2
        ISTATE=0
        DO MSA=-ISA,ISA,2
        DO MIA=-INUCA,INUCA,2
          MSBMAX=ISB
          IF (IDENTN) MSBMAX=MSA
        DO MSB=-ISB,MSBMAX,2
          MIBMAX=INUCB
          IF (IDENTN .AND. MSA.EQ.MSB) MIBMAX=MIA
        DO MIB=-INUCB,MIBMAX,2
          ISTATE=ISTATE+1
          IF (ILOOP.EQ.2) THEN
            JSTATE(ISTATE         )=MSA
            JSTATE(ISTATE+NSTATE  )=MIA
            JSTATE(ISTATE+NSTATE*2)=MSB
            JSTATE(ISTATE+NSTATE*3)=MIB
            JSTATE(ISTATE+NSTATE*4)=ISTATE
          ENDIF
        ENDDO
        ENDDO
        ENDDO
        ENDDO
        IF (ILOOP.EQ.1) NSTATE=ISTATE
      ENDDO
C
      RETURN
C ========================================================== END OF SET9
C
      ENTRY BASE9(LCOUNT,N,JTOT,IBLOCK,JSTATE,NSTATE,NQN,JSINDX,L,
     1            IPRINT)
C
C  UNCHANGED FROM base9-alk_alk_ucpld.f (see CONTEXT_Tm_CPL9.md section
C  6: "BASE9 -- переносится почти без правок").  JTOT plays the role
C  of 2*Mtot ; IBLOCK selects the parity block.
C
      IF ((MOD(JTOT,2).EQ.0).NEQV.(MOD(ISA+INUCA+ISB+INUCB,2).EQ.0))
     1THEN
        WRITE(6,*) ' *** ERROR. JTOT MUST HAVE THE SAME PARITY AS',
     1             ' ISA + INUCA + ISB + INUCB. STOPPING'
        STOP
      ENDIF
C
      MTOT=JTOT
C
      IF (LCOUNT) THEN
        IF (IPRINT.GE.1) WRITE(6,605) IBLOCK,(-1)**IBLOCK,MTOT
  605   FORMAT('  SYMMETRY BLOCK = ',I3,' SELECTS PARITY',I3,/
     1         '  MTOT =',I3,'/2')
        IF (IDENTN) THEN
          IFA=INUCA+ISA
          IF (2*(IFA/2).EQ.IFA) THEN
            IBOSFR=0
            IF (IPRINT.GE.1) WRITE(6,610) 'BOSONS'
  610         FORMAT(2X,'BASIS SET GENERATED FOR TWO IDENTICAL ',A)
          ELSE
            IBOSFR=1
            IF (IPRINT.GE.1) WRITE(6,610) 'FERMIONS'
          ENDIF
        ENDIF
      ENDIF
C
      N=0
      DO I=1,NSTATE
        MSA=JSTATE(I)
        MIA=JSTATE(NSTATE+I)
        MSB=JSTATE(2*NSTATE+I)
        MIB=JSTATE(3*NSTATE+I)
        MF=MSA+MSB+MIA+MIB
        ML=MTOT-MF
        LSTART=4-2*IBLOCK
        IF (IDENTN .AND. MSA.EQ.MSB .AND. MIA.EQ.MIB .AND.
     1      IBOSFR+IBLOCK.NE.2) LSTART=2*LMAX+1
        DO LL=LSTART,2*LMAX,4
          IF (ABS(ML).GT.LL) CYCLE
          IF (NREQ.EQ.0) GOTO 300
          DO IREQ=1,NREQ
            IF (LREQ(IREQ).GE.0 .AND. LL/2.NE.LREQ(IREQ)) CYCLE
            IF (ABS(MFREQ(IREQ)).LT.999 .AND. MF.NE.MFREQ(IREQ))
     1        CYCLE
            GOTO 300
          ENDDO
          CYCLE
  300     N=N+1
          IF (LCOUNT) CYCLE
          JSINDX(N)=I
          L(N)=LL/2
        ENDDO
      ENDDO
C
      IF (LCOUNT) RETURN
C
      DO I=1,N
        DO JJ=I+1,N
          IF (L(JJ).LT.L(I)) THEN
            LJ=L(I)
            L(I)=L(JJ)
            L(JJ)=LJ
            LJ=JSINDX(I)
            JSINDX(I)=JSINDX(JJ)
            JSINDX(JJ)=LJ
          ENDIF
        ENDDO
      ENDDO
C
      RETURN
C ========================================================= END OF BASE9
C
      ENTRY CPL9(N,IBLOCK,NHAM,LAM,MXLAM,NSTATE,JSTATE,JSINDX,L,JTOT,
     1           VL,IV,CENT,DGVL,IBOUND,IEXCH,IPRINT)
C
C  ===================================================================
C  CPL9 -- THE MAIN REWRITE.  Block layout (MXLAM=3):
C     LL = 1            : isotropic potential, k=0,i=1  (V^(1)_0(R))
C     LL = 2            : rank-2 potential,    k=2,i=1  (V^(1)_2(R))
C     LL = 3            : magnetic dipole-dipole, k=2,i=2 (V^(2)_2(R)
C                         = D^(2)_2/R^3, via DIPBLK/TENSX below)
C     LL = MXLAM+1 = 4  : hyperfine term (unchanged from alkali code,
C                         SDOTI2 is already generic in j,i)
C     LL = MXLAM+2 = 5  : Zeeman term (unchanged, linear in gJ*mJ)
C     LL = MXLAM+NCONST+1,+2 : extra operators (unchanged, mF^2, mF)
C  ===================================================================
C
      DO LL=1,NVLBLK
        NNZ=0
        I=LL
        DO ICOL=1,N
          MSAC=JSTATE(JSINDX(ICOL))
          MIAC=JSTATE(JSINDX(ICOL)+NSTATE)
          MSBC=JSTATE(JSINDX(ICOL)+NSTATE*2)
          MIBC=JSTATE(JSINDX(ICOL)+NSTATE*3)
          LC=L(ICOL)
        DO IROW=1,ICOL
          MSAR=JSTATE(JSINDX(IROW))
          MIAR=JSTATE(JSINDX(IROW)+NSTATE)
          MSBR=JSTATE(JSINDX(IROW)+NSTATE*2)
          MIBR=JSTATE(JSINDX(IROW)+NSTATE*3)
          MLR=MTOT-MSAR-MSBR-MIAR-MIBR
          MLC=MTOT-MSAC-MSBC-MIAC-MIBC
          LR=L(IROW)
          VL(I)=0.D0
          PREFAC=1.D0
          IF (IDENTN .AND. MSAR.EQ.MSBR .AND. MIAR.EQ.MIBR)
     1      PREFAC=PREFAC/SQRT(2.D0)
          IF (IDENTN .AND. MSAC.EQ.MSBC .AND. MIAC.EQ.MIBC)
     1      PREFAC=PREFAC/SQRT(2.D0)
C
          IF (LL.LE.2) THEN
C  POTENTIAL BLOCKS, i=1 FAMILY (k=0 for LL=1, k=2 for LL=2) ============
C  Nuclear-spin projections must be unchanged (potential acts only on
C  the electronic/orbital parts); tensor rank K=2*(LL-1).
            K=2*(LL-1)
            IF (MIAR.EQ.MIAC .AND. MIBR.EQ.MIBC) THEN
              FAC=PREFAC*POTBLK(K,MSAR,MIAR,MSBR,MIBR,LR,MLR,
     1                          MSAC,MIAC,MSBC,MIBC,LC,MLC,ISA)
              VL(I)=FAC
            ENDIF
C  identical-particle exchange term (B<->A swap on the bra)
            IF (IDENTN .AND. MIAR.EQ.MIBC .AND. MIBR.EQ.MIAC) THEN
              FAC2=PREFAC*POTBLK(K,MSBR,MIBR,MSAR,MIAR,LR,MLR,
     1                           MSAC,MIAC,MSBC,MIBC,LC,MLC,ISA)
              IF (MOD(IBOSFR+LR,2).NE.0) FAC2=-FAC2
              VL(I)=VL(I)+FAC2
            ENDIF
C
          ELSEIF (LL.EQ.3 .AND. MXLAM.GE.3) THEN
C  POTENTIAL BLOCK, i=2 FAMILY, k=2: MAGNETIC DIPOLE-DIPOLE ============
C  Tˆ(2)_2q = [j1(x)j2]_2q (paper Eq (6)); radial strength V^(2)_2(R)
C  = D^(2)_2/R^3 supplied as an analytic term via &POTL LAMBDA block 3
C  (NTERM>0, no VSTAR call -- see POTIN9 below).  Same nuclear-spin and
C  exchange handling as the i=1 blocks above.
C  MXLAM.GE.3 GUARD IS LOAD-BEARING: with the OLD MXLAM=2 &POTL (i=1
C  family only, e.g. tm_resonance_pipeline.py's current template),
C  MXLAM+1 also equals 3 -- without this guard this branch would win
C  the ELSEIF race against the HYPERFINE branch below and silently
C  replace hyperfine coupling with a bogus read of an undeclared 3rd
C  LAMBDA block.  Only activates for the NEW MXLAM=3 &POTL.
            IF (MIAR.EQ.MIAC .AND. MIBR.EQ.MIBC) THEN
              FAC=PREFAC*DIPBLK(2,MSAR,MIAR,MSBR,MIBR,LR,MLR,
     1                          MSAC,MIAC,MSBC,MIBC,LC,MLC,ISA)
              VL(I)=FAC
            ENDIF
            IF (IDENTN .AND. MIAR.EQ.MIBC .AND. MIBR.EQ.MIAC) THEN
              FAC2=PREFAC*DIPBLK(2,MSBR,MIBR,MSAR,MIAR,LR,MLR,
     1                           MSAC,MIAC,MSBC,MIBC,LC,MLC,ISA)
              IF (MOD(IBOSFR+LR,2).NE.0) FAC2=-FAC2
              VL(I)=VL(I)+FAC2
            ENDIF
C
          ELSEIF (LL.EQ.MXLAM+1 .AND. LR.EQ.LC .AND. MLR.EQ.MLC) THEN
C  HYPERFINE  ========================================================
C  UNCHANGED from parent file: SDOTI2 is generic in (ISA,INUCA).
            FACA=ANSA*SDOTI2(ISA,MSAR,MSAC,INUCA,MIAR,MIAC)
            FACB=ANSB*SDOTI2(ISB,MSBR,MSBC,INUCB,MIBR,MIBC)
            IF (MSBR.EQ.MSBC .AND. MIBR.EQ.MIBC)
     1        VL(I)=PREFAC*FACA
            IF (MSAR.EQ.MSAC .AND. MIAR.EQ.MIAC)
     1        VL(I)=VL(I)+PREFAC*FACB
            IF (IDENTN) THEN
              IF (MSAR.EQ.MSBC .AND. MIAR.EQ.MIBC) THEN
                FAC2=PREFAC*ANSA*SDOTI2(ISA,MSBR,MSAC,INUCA,MIBR,MIAC)
                IF (MOD(IBOSFR+LR,2).NE.0) FAC2=-FAC2
                VL(I)=VL(I)+FAC2
              ENDIF
              IF (MSBR.EQ.MSAC .AND. MIBR.EQ.MIAC) THEN
                FAC2=PREFAC*ANSA*SDOTI2(ISA,MSAR,MSBC,INUCA,MIAR,MIBC)
                IF (MOD(IBOSFR+LR,2).NE.0) FAC2=-FAC2
                VL(I)=VL(I)+FAC2
              ENDIF
            ENDIF
C
          ELSEIF (LL.EQ.MXLAM+2 .AND. ICOL.EQ.IROW) THEN
C  ZEEMAN, LINEAR IN gJ*mJ (UNCHANGED, generic in ISA) ================
            VL(I)=GSA*DBLE(MSAR)+GSB*DBLE(MSBR)
     1            +GA*DBLE(MIAR)+GB*DBLE(MIBR)
            VL(I)=VL(I)*BM*0.5D0
C
          ELSEIF (LL.EQ.MXLAM+NCONST+1 .AND. ICOL.EQ.IROW) THEN
            VL(I)=DBLE((MSAR+MIAR)**2+(MSBR+MIBR)**2)
          ELSEIF (LL.EQ.MXLAM+NCONST+2 .AND. ICOL.EQ.IROW) THEN
            VL(I)=DBLE(MSAR+MIAR+MSBR+MIBR)
          ENDIF
C
          IF (VL(I).NE.0.D0) NNZ=NNZ+1
          I=I+NVLBLK
        ENDDO
        ENDDO
C
        IF (NNZ.EQ.0) THEN
          WRITE(6,612) JTOT,LL
  612     FORMAT('  * * * NOTE.  FOR JTOT =',I4,',  ALL COUPLING',
     1           ' COEFFICIENTS ARE 0.0 FOR POTENTIAL SYMMETRY',I4)
        ENDIF
      ENDDO
C
      RETURN
C ========================================================= END OF CPL9
C
      ENTRY THRSH9(JREF,MONQN,NQN1,EREF,IPRINT)
C
C  =====================================================================
C  Breit-Rabi-type closed form for the single-atom threshold energy of
C  Tm (electronic J = ISA/2, ARBITRARY, coupled to nuclear I = INUCA/2
C  = 1/2, FIXED).  Re-derived (not copy-pasted) from the same 2-level
C  Hamiltonian H = ANSA*I.J + (GSA*mJ+GA*mI)*mu_B*B that CPL9 builds,
C  because here the roles of the two coupled angular momenta are
C  SWAPPED relative to the parent alkali file:
C     parent (base9-alk_alk_ucpld.f): "spin-1/2 partner" = ELECTRON
C       (fixed s=1/2, g-factor GSA), "arbitrary partner"  = NUCLEUS
C       (I, g-factor GA).  Their x = (GSA-GA)*mu_B*B/HFSPLA.
C     here:  "spin-1/2 partner" = NUCLEUS (I=1/2, g-factor GA),
C            "arbitrary partner" = ELECTRON (J=ISA/2, g-factor GSA).
C       => the linear (outside-sqrt) term now carries GSA, not GA;
C          and BX = (GA-GSA)*mu_B*B/HFSPLA  (note the SIGN FLIP
C          relative to the parent file's BX -- this is not a typo).
C  For a general two-body Hamiltonian H = Lambda*A.B + g_A*muB*B*m_A
C  + g_B*muB*B*m_B with A of fixed spin 1/2 and B arbitrary (spin K),
C  diagonalising the 2x2 block at fixed mF = m_A+m_B gives exactly:
C     E(F=K+-1/2,mF) = -HFSPLA/(2*(2K+1))+g_B*muB*B*mF
C                       +-(HFSPLA/2)*SQRT(1+2*mF*X/(K+1/2)+X**2)
C     X = (g_A-g_B)*muB*B/HFSPLA
C  which is what is coded below with A=nucleus (I=1/2), B=electron
C  (J).  Checked at B=0 against E(F)=(A/2)[F(F+1)-I(I+1)-J(J+1)] (with
C  A=ANSA=2*HFSPLA/NSFAC): reproduces AJ/2 for F=J+1/2 and -A(J+1)/2
C  for F=J-1/2 exactly, for EITHER sign of HFSPLA/ANSA.  The extremal
C  (single-state, |mF|=J+1/2) branch is likewise re-derived with ISA
C  (not INUCA) as the "arbitrary partner" in the first term; the
C  Zeeman part of the extremal formula (GSA*ISA*0.5+GA*INUCA*0.5) is
C  unchanged, since a single fully-stretched product state has no
C  role-dependence there.
C
C     MONQN(1) = 2*F(A) that the desired state correlates with at B->0
C     MONQN(2) = 2*mF(A)
C     MONQN(3) = 2*F(B)
C     MONQN(4) = 2*mF(B)
C  =====================================================================
C
      BFIELD=EFV(1)
      IF (JREF.GT.0) THEN
        WRITE(6,*) ' *** ERROR - THRSH9 CALLED WITH POSITIVE IREF'
        STOP
      ENDIF
C
      IF (MONQN(1).EQ.-99999) THEN
        WRITE(6,*) ' *** ERROR - THRSH9 CALLED WITH MONQN UNSET'
        STOP
      ENDIF
C
      BOHRM=BM*GHZCM
C
C  BREIT-RABI FOR ATOM A ("arbitrary partner" = J = ISA/2 ;
C  "spin-1/2 partner" = I = INUCA/2 = 1/2)
C
      M=MONQN(2)
      IF (ABS(MONQN(1)-ISA).NE.1) THEN
        WRITE(6,*) ' *** THRSH9: INVALID MONQN(1) =',MONQN(1)
        STOP
      ELSEIF (ABS(M).GT.MONQN(1)) THEN
        WRITE(6,*) ' *** THRSH9: MA =',M,' > FA. STOPPING'
        STOP
      ELSEIF (MOD(M+MONQN(1),2).NE.0) THEN
        WRITE(6,*) ' *** THRSH9: INVALID MONQN(1),MONQN(2) PAIR =',
     1             MONQN(1),M
        STOP
      ENDIF
C
      E1=-HFSPLA/(2.D0*DBLE(ISA+1)) + 0.5D0*GSA*BOHRM*DBLE(M)*BFIELD
      BX=BOHRM*BFIELD*(GA-GSA)/HFSPLA
      E2=0.5D0*HFSPLA*SQRT(1.D0+DBLE(M+M)*BX/DBLE(ISA+1)+BX*BX)
C
      IF (ABS(M).EQ.ISA+1) THEN
        EA=DBLE(ISA)*HFSPLA/(2.D0*DBLE(ISA+1))+SIGN(1.D0,DBLE(M))
     1     *BOHRM*BFIELD*(GSA*DBLE(ISA)*0.5D0+GA*DBLE(INUCA)*0.5D0)
      ELSEIF (MONQN(1).EQ.ISA+1) THEN
        EA=E1+E2
      ELSEIF (MONQN(1).EQ.ISA-1) THEN
        EA=E1-E2
      ENDIF
C
C  BREIT-RABI FOR ATOM B (identical structure, using ISB/INUCB/etc)
C
      M=MONQN(4)
      IF (ABS(MONQN(3)-ISB).NE.1) THEN
        WRITE(6,*) ' *** THRSH9: INVALID MONQN(3) =',MONQN(3)
        STOP
      ELSEIF (ABS(M).GT.MONQN(3)) THEN
        WRITE(6,*) ' *** THRSH9: MB =',M,' > FB. STOPPING'
        STOP
      ELSEIF (MOD(M+MONQN(3),2).NE.0) THEN
        WRITE(6,*) ' *** THRSH9: INVALID MONQN(3),MONQN(4) PAIR =',
     1             MONQN(3),M
        STOP
      ENDIF
C
      E1=-HFSPLB/(2.D0*DBLE(ISB+1)) + 0.5D0*GSB*BOHRM*DBLE(M)*BFIELD
      BX=BOHRM*BFIELD*(GB-GSB)/HFSPLB
      E2=0.5D0*HFSPLB*SQRT(1.D0+DBLE(M+M)*BX/DBLE(ISB+1)+BX*BX)
C
      IF (ABS(M).EQ.ISB+1) THEN
        EB=DBLE(ISB)*HFSPLB/(2.D0*DBLE(ISB+1))+SIGN(1.D0,DBLE(M))
     1     *BOHRM*BFIELD*(GSB*DBLE(ISB)*0.5D0+GB*DBLE(INUCB)*0.5D0)
      ELSEIF (MONQN(3).EQ.ISB+1) THEN
        EB=E1+E2
      ELSEIF (MONQN(3).EQ.ISB-1) THEN
        EB=E1-E2
      ENDIF
C
      IF (IPRINT.GE.8) THEN
        WRITE(6,*)
        WRITE(6,667) 'A',MONQN(1),MONQN(2),EA
        WRITE(6,667) 'B',MONQN(3),MONQN(4),EB
  667   FORMAT('  ATOM ',A1,' WITH DOUBLED QUANTUM NOS',2I3,
     1     ' IS AT ENERGY',F12.7,' GHZ')
      ENDIF
      EAB=EA+EB
      EREF=EAB/GHZCM
C
      RETURN
C ======================================================== END OF THRSH9
      ENTRY POTIN9(ITYPE,LAM,MXLAM,NPTS,NDIM,XPT,XWT,
     1             MXPT,IVMIN,IVMAX,L1MAX,L2MAX,
     2             MXLMB,X,MX,IXFAC)
C
C  UNCHANGED IN SPIRIT from the parent file: reuse ITYP=1's generic
C  machinery to read MXLAM radial coefficients v_lambda(R), each either
C  from VSTAR (pot-Tm2.f, NTERM<0) or as analytic power/exponential
C  terms given directly in &POTL (NTERM>0).  MXLAM must equal 3 in
C  &POTL (LMAX=1), ordered lambda(1) -> V^(1)_0(R) [VSTAR I=1],
C  lambda(2) -> V^(1)_2(R) [VSTAR I=2], lambda(3) -> V^(2)_2(R) =
C  D^(2)_2/R^3 [NTERM=1 analytic term, NPOWER=-3, no VSTAR call] --
C  this ORDER MUST MATCH the LL=1,2,3 block numbering used in CPL9.
C  For 169Tm2 (R in the units used by pot-Tm2.f, i.e. Angstrom):
C     D^(2)_2 = -1.38117889454613 cm^-1 Angstrom^3
C  from D^(2)_2 = -sqrt(6)*alpha^2*(g_j/2)^2*E_h*a0^3 (paper, Sec 2.3,
C  with g_j=GSA=1.14119), converted a0^3 -> Angstrom^3 via
C  bohr_to_Angstrom^3, using the SAME physical_constants values linked
C  by the Makefile (physical_constants_module.f, 2022 CODATA/BIPM:
C  inverse_fine_structure_constant=137.035999177, hartree_in_inv_cm=
C  2.1947463136314D5, bohr_in_SI/Angstrom_in_SI=0.529177210544).
C  Example &POTL entry for this 3rd block:
C     MXLAM=3, LAMBDA=0,1,2, NTERM=-1,-1,1,
C     NPOWER(1)=-3, A(1)=-1.38117889454613
      ITYPE=1
C
      RETURN
      END
C ======================================================== END OF POTIN9
      FUNCTION SDOTI2(IS,MS1,MS2,II,MI1,MI2)
C  UNCHANGED FROM base9-alk_alk_ucpld.f.  Generic dot-product matrix
C  element <j m1|J.I|j m2> in the uncoupled basis for ANY j,i (not
C  specific to spin-1/2) -- valid as-is for j=7/2.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SDOTI2=0.D0
      IF (MS1+MI1.EQ.MS2+MI2) THEN
        IF (MS1.EQ.MS2) THEN
          SDOTI2=0.25D0*DBLE(MI1*MS1)
        ELSEIF (ABS(MS1-MS2).EQ.2) THEN
          SDOTI2=0.125D0*SQRT(DBLE(II*(II+2)-MI1*MI2))
     1                  *SQRT(DBLE(IS*(IS+2)-MS1*MS2))
        ENDIF
      ENDIF
      RETURN
      END
C ======================================================= END OF SDOTI2
C  =====================================================================
C  NOTE ON DUPLICATION: POTBLK/TENSOR_ME/TENS1/VECME/CG3/CKME below are
C  intentionally copied verbatim from tens_ck_funcs.f so this file is
C  self-contained for the molscat-Tm2 build.  DO NOT also link
C  tens_ck_funcs.f (or test_tens_ck.f) into the SAME executable as this
C  file -- that will raise "duplicate symbol" errors at link time.
C  tens_ck_funcs.f is for the STAGE A standalone unit test ONLY
C  (test_tens_ck.f); keep the two builds separate.
C  DIPBLK/TENSX further below are NEW (not in tens_ck_funcs.f/
C  test_tens_ck.f) -- they implement the i=2 cross-tensor [j1(x)j2]_kq,
C  reusing VECME/CG3/CKME from the block above.
C  =====================================================================
      FUNCTION POTBLK(K,MSAR,MIAR,MSBR,MIBR,LR,MLR,
     1                 MSAC,MIAC,MSBC,MIBC,LC,MLC,J2)
C  =====================================================================
C  POTBLK -- matrix element of the k-th rank tensor block of the Eq.(1)
C  potential expansion (i=1 family only: T^(1)_kq = [j1(x)j1]_kq +
C  [j2(x)j2]_kq), for the TWO-ATOM uncoupled pair state, SUMMED OVER q:
C
C    <row| V^(1)_k block |col> = Sum_q (-1)^q TENSOR_ME(k,q) * CKME(k,-q)
C
C  Nuclear-spin projections are NOT arguments of the tensor operators
C  themselves; the caller (CPL9) has already imposed MIAR=MIAC,
C  MIBR=MIBC before calling this function.  All spin/electronic
C  projections (MSAR etc.) are DOUBLED integers (2*mj); LR,LC are
C  UNDOUBLED orbital angular momenta; MLR,MLC are DOUBLED (2*ml, always
C  even since ml is an integer).  J2 = 2*j (=ISA=ISB=7 for Tm).
C  =====================================================================
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      POTBLK=0.D0
      IF (K.EQ.0) THEN
C  isotropic block: trivial, diagonal in everything
        IF (MSAR.EQ.MSAC .AND. MSBR.EQ.MSBC .AND. LR.EQ.LC
     1      .AND. MLR.EQ.MLC) POTBLK=1.D0
        RETURN
      ENDIF
C  general k (2, 4, ...): sum over q = -k..k (integer)
      DO IQ=-K,K
C  selection rule on electronic side (checked again inside TENSOR_ME,
C  but pre-filter here to skip the CKME call when possible)
        IF (MSAR+MSBR.NE.MSAC+MSBC+2*IQ) CYCLE
        TME=TENSOR_ME(K,IQ,MSAR,MSBR,MSAC,MSBC,J2)
        IF (TME.EQ.0.D0) CYCLE
C  ml/l orbital element for C_{k,-q};  ml is undoubled here (MLR/2 etc,
C  both guaranteed even)
        CKM=CKME(K,-IQ,LR,MLR/2,LC,MLC/2)
        IF (CKM.EQ.0.D0) CYCLE
        SGN=1.D0
        IF (MOD(ABS(IQ),2).EQ.1) SGN=-1.D0
        POTBLK=POTBLK+SGN*TME*CKM
      ENDDO
      RETURN
      END
C ======================================================== END OF POTBLK
      FUNCTION TENSOR_ME(K,Q,MJAR,MJBR,MJAC,MJBC,J2)
C  =====================================================================
C  <mjAR mjBR| T^(1)_KQ |mjAC mjBC> = delta(mjBR,mjBC)*TENS1(mjA part)
C                                    + delta(mjAR,mjAC)*TENS1(mjB part)
C  where T^(1)_Kq = [j1(x)j1]_Kq + [j2(x)j2]_Kq  (single-atom rank-K
C  tensors, summed).  All Mj arguments DOUBLED; J2=2*j.
C  Verified against Python couple()-kron reference to 3.6e-15 (see
C  step6_fortran_algo_check.py / tm_tensors_full.py).
C  =====================================================================
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C  CRLS/AI FIX: 'Q' falls in the O-Z range, so IMPLICIT above would make
C  it DOUBLE PRECISION by default even though it is a true (undoubled)
C  integer projection index everywhere it is used/called.  Must override
C  explicitly, or every call site passing an INTEGER actual argument to
C  this DOUBLE PRECISION dummy is a silent argument-type mismatch.
      INTEGER Q
      TENSOR_ME=0.D0
      IF (MJBR.EQ.MJBC) TENSOR_ME=TENSOR_ME+TENS1(K,Q,J2,MJAR,MJAC)
      IF (MJAR.EQ.MJAC) TENSOR_ME=TENSOR_ME+TENS1(K,Q,J2,MJBR,MJBC)
      RETURN
      END
C ===================================================== END OF TENSOR_ME
      FUNCTION TENS1(K,Q,J2,MR2,MC2)
C  =====================================================================
C  Single-atom rank-K spherical tensor matrix element
C     <j mR|[j(x)j]_KQ|j mC>
C  built from the rank-1 spherical components of the angular-momentum
C  VECTOR operator itself, coupled to rank K via a Clebsch-Gordan sum:
C     [j(x)j]_KQ = Sum_{q1+q2=Q} <1 q1 1 q2|K Q> J_q1 J_q2
C  Since each J_q is itself a ladder/z operator, only ONE intermediate
C  state contributes per (q1,q2) pair -- see VECME below.
C  All M arguments DOUBLED (2*mR, 2*mC); J2=2*j; K,Q undoubled.
C  =====================================================================
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C  CRLS/AI FIX: see note in TENSOR_ME above -- 'Q' must be forced INTEGER.
      INTEGER Q
      TENS1=0.D0
      IF (MR2-MC2.NE.2*Q) RETURN
      DO IQ1=-1,1
        IQ2=Q-IQ1
        IF (ABS(IQ2).GT.1) CYCLE
        MI2=MC2+2*IQ2
        IF (ABS(MI2).GT.J2) CYCLE
        V2=VECME(IQ2,J2,MI2,MC2)
        IF (V2.EQ.0.D0) CYCLE
        V1=VECME(IQ1,J2,MR2,MI2)
        IF (V1.EQ.0.D0) CYCLE
        CG=CG3(2,2*IQ1,2,2*IQ2,2*K,2*Q)
        TENS1=TENS1+CG*V1*V2
      ENDDO
      RETURN
      END
C ========================================================= END OF TENS1
      FUNCTION VECME(IQ,J2,MOUT2,MIN2)
C  Spherical components of a vector (rank-1) angular-momentum operator,
C  standard Racah normalisation:  J_(+-1) = -+(Jx+-iJy)/sqrt(2), J_0=Jz.
C  All M arguments DOUBLED; J2=2*j.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      VECME=0.D0
      XJ=DBLE(J2)/2.D0
      XMIN=DBLE(MIN2)/2.D0
      IF (IQ.EQ.0) THEN
        IF (MOUT2.EQ.MIN2) VECME=XMIN
      ELSEIF (IQ.EQ.1) THEN
        IF (MOUT2.EQ.MIN2+2)
     1    VECME=-SQRT(XJ*(XJ+1.D0)-XMIN*(XMIN+1.D0))/SQRT(2.D0)
      ELSEIF (IQ.EQ.-1) THEN
        IF (MOUT2.EQ.MIN2-2)
     1    VECME=SQRT(XJ*(XJ+1.D0)-XMIN*(XMIN-1.D0))/SQRT(2.D0)
      ENDIF
      RETURN
      END
C ========================================================= END OF VECME
      FUNCTION CG3(J1X2,M1X2,J2X2,M2X2,J3X2,M3X2)
C  Clebsch-Gordan coefficient <j1 m1 j2 m2|j3 m3>, expressed via the
C  standard MOLSCAT library Wigner-3j function THRJ (real, undoubled
C  arguments), exactly as done for the alkali spin-spin term (SPINSP
C  uses THRJ the same way).  All arguments here DOUBLED integers.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CG3=0.D0
      IF (M1X2+M2X2.NE.M3X2) RETURN
      XJ1=DBLE(J1X2)/2.D0
      XJ2=DBLE(J2X2)/2.D0
      XJ3=DBLE(J3X2)/2.D0
      XM1=DBLE(M1X2)/2.D0
      XM2=DBLE(M2X2)/2.D0
      XM3=DBLE(M3X2)/2.D0
      T=THRJ(XJ1,XJ2,XJ3,XM1,XM2,-XM3)
      IF (T.EQ.0.D0) RETURN
      CG3=PARSGN(NINT(XJ1-XJ2+XM3))*SQRT(XJ3+XJ3+1.D0)*T
      RETURN
      END
C =========================================================== END OF CG3
      FUNCTION CKME(K,Q,LR,MLR,LC,MLC)
C  =====================================================================
C  Orbital matrix element  <LR MLR|C_KQ|LC MLC>  of the unnormalised
C  spherical tensor C_Kq = sqrt(4pi/(2K+1)) Y_Kq, via the standard
C  Gaunt-type closed form (Racah):
C     <lR mR|C_Kq|lC mC> = (-1)^mR sqrt((2lR+1)(2lC+1))
C                           * 3j(lR K lC;0 0 0) * 3j(lR K lC;-mR q mC)
C  THREEJ(L1,L2,L3) is the MOLSCAT library function returning the
C  zero-projection 3j symbol 3j(L1,L2,L3;0,0,0) for integer ranks
C  (used with L2=2 in the parent file's SPINSP; here used generally
C  with L2=K=0,2,4).  All arguments here UNDOUBLED (K,Q,LR,MLR,LC,MLC
C  are true integers: L never doubled, ML passed in already /2 by the
C  caller POTBLK).
C  Verified against sympy.physics.wigner.gaunt to 3.3e-16 in Python
C  (tm_tensors_full.py, step2_ck_me.py) before this translation.
C  =====================================================================
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C  CRLS/AI FIX: see note in TENSOR_ME above -- 'Q' must be forced INTEGER.
      INTEGER Q
      CKME=0.D0
      IF (MLR.NE.MLC+Q) RETURN
      T1=THREEJ(LR,K,LC)
      IF (T1.EQ.0.D0) RETURN
      T2=THRJ(DBLE(LR),DBLE(K),DBLE(LC),-DBLE(MLR),DBLE(Q),DBLE(MLC))
      IF (T2.EQ.0.D0) RETURN
      CKME=PARSGN(MLR)*SQRT(DBLE((2*LR+1)*(2*LC+1)))*T1*T2
      RETURN
      END
C ========================================================== END OF CKME
      FUNCTION DIPBLK(K,MSAR,MIAR,MSBR,MIBR,LR,MLR,
     1                 MSAC,MIAC,MSBC,MIBC,LC,MLC,J2)
C  =====================================================================
C  DIPBLK -- matrix element of the k-th rank CROSS tensor block (i=2
C  family: Tˆ(2)_kq = [j1(x)j2]_kq, Eq (6) of Tiesinga et al. NJP 23,
C  085007 (2021)), for the TWO-ATOM uncoupled pair state, SUMMED OVER q:
C
C    <row| V^(2)_k block |col> = Sum_q (-1)^q TENSX(k,q,...) * CKME(k,-q)
C
C  Structurally identical to POTBLK above, but calling TENSX (the
C  cross-atom [j1(x)j2] tensor) instead of TENSOR_ME (the same-atom
C  [j1(x)j1]+[j2(x)j2] tensor) -- see POTBLK for the argument and unit
C  conventions, which are identical here.  Unlike POTBLK there is no
C  K=0 shortcut: Tˆ(2)_00=[j1(x)j2]_00 is the spin-exchange j1.j2 term,
C  not the identity, but the general Sum_q loop below handles K=0
C  correctly anyway if it is ever needed (not currently called with
C  K=0 -- only K=2 is used, for the magnetic dipole-dipole term).
C  =====================================================================
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIPBLK=0.D0
      DO IQ=-K,K
        IF (MSAR+MSBR.NE.MSAC+MSBC+2*IQ) CYCLE
        TME=TENSX(K,IQ,J2,MSAR,MSAC,J2,MSBR,MSBC)
        IF (TME.EQ.0.D0) CYCLE
        CKM=CKME(K,-IQ,LR,MLR/2,LC,MLC/2)
        IF (CKM.EQ.0.D0) CYCLE
        SGN=1.D0
        IF (MOD(ABS(IQ),2).EQ.1) SGN=-1.D0
        DIPBLK=DIPBLK+SGN*TME*CKM
      ENDDO
      RETURN
      END
C ========================================================= END OF DIPBLK
      FUNCTION TENSX(K,Q,J2A,MAR2,MAC2,J2B,MBR2,MBC2)
C  =====================================================================
C  Cross-atom rank-K spherical tensor matrix element
C     <j1A mAR; j2B mBR|[j1(x)j2]_KQ|j1A mAC; j2B mBC>
C  built exactly like TENS1 (single-atom [j(x)j]_KQ) but with the two
C  rank-1 vector factors drawn from DIFFERENT atoms:
C     [j1(x)j2]_KQ = Sum_{q1+q2=Q} <1 q1 1 q2|K Q> J1_q1 J2_q2
C  Since J1_q1 acts only on atom A and J2_q2 only on atom B (they act
C  on different tensor factors of the pair Hilbert space and trivially
C  commute), the matrix element factorises directly -- no intermediate-
C  state sum over a shared single-particle basis is needed here (unlike
C  TENS1, which composes two ladder operators on the SAME atom).
C  Uses the SAME VECME/CG3 primitives as TENS1, so this is automatically
C  in the identical Racah/Brink&Satchler normalisation used throughout
C  this file for the i=1 tensor family (paper Eq (5)-(7)) -- required
C  so that the long-range coefficients D^(2)_2, C^(2)_2 (paper Table 1,
C  Sec 2.3) can be used directly, with no extra conversion factor
C  between literature conventions.
C  All M arguments DOUBLED (2*mR, 2*mC); J2A,J2B=2*j of each atom
C  (equal, =ISA, for the homonuclear Tm2 case this file is built for);
C  K,Q undoubled.
C  =====================================================================
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C  CRLS/AI FIX: see note in TENSOR_ME above -- 'Q' must be forced INTEGER.
      INTEGER Q
      TENSX=0.D0
      IF (MAR2-MAC2+MBR2-MBC2.NE.2*Q) RETURN
      DO IQ1=-1,1
        IQ2=Q-IQ1
        IF (ABS(IQ2).GT.1) CYCLE
        VA=VECME(IQ1,J2A,MAR2,MAC2)
        IF (VA.EQ.0.D0) CYCLE
        VB=VECME(IQ2,J2B,MBR2,MBC2)
        IF (VB.EQ.0.D0) CYCLE
        CG=CG3(2,2*IQ1,2,2*IQ2,2*K,2*Q)
        TENSX=TENSX+CG*VA*VB
      ENDDO
      RETURN
      END
C ========================================================== END OF TENSX