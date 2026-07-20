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
