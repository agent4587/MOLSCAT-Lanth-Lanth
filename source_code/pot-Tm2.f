      SUBROUTINE VINIT(I,RM,EPSIL)
C  =====================================================================
C  pot-Tm2.f  --  RADIAL STRENGTHS V^(1)_0(R), V^(1)_2(R) FOR Tm2
C  =====================================================================
C  Source: Tiesinga, Klos, Li, Petrov, Kotochigova, New J. Phys. 23,
C  085007 (2021).  Table 5 (ab initio DIRRCI values, 7.0 <= R <= 12.0
C  a0) + Table 1 (long-range van-der-Waals coefficients C^(1)_0,
C  C^(1)_2, in Eh*a0^6) for the i=1 tensor family only (k=0 isotropic,
C  k=2 rank-2).  Used via MOLSCAT's classic ITYP=1 VINIT/VSTAR
C  mechanism (NTERM(I)<0 in &POTL), as set up by POTIN9 in
C  base9-Tm_Tm_AI.f: I=1 -> lambda=0 -> V^(1)_0(R); I=2 -> lambda=1 ->
C  V^(1)_2(R).  MXLAM=2, LMAX=1 required in &POTL.
C
C  METHOD:
C   - R in [7.0, 12.0] a0: table 5 values, exactly reproduced.
C   - Interpolation is done on R^6*V(R) (varies far more smoothly than
C     V(R) itself, exactly as recommended in the paper's Appendix C)
C     using a MONOTONE cubic Hermite spline (Fritsch-Carlson / PCHIP).
C     A plain natural cubic spline was tried first and produced
C     spurious ~10% overshoots on the unevenly-spaced tail of the
C     table (R=10.4,11.0,12.0 a0) -- PCHIP does not have this problem
C     because it is constructed to avoid overshoot by design, and it
C     agrees with the Akima spline the original paper uses to <1% here
C     (checked numerically against scipy.interpolate.Akima1DInterpolator
C     before writing this file; PCHIP was then checked line-by-line
C     against scipy.interpolate.PchipInterpolator, see pchip_check.py).
C   - R in (12.0, Rrel] with Rrel=12.0+0.5=12.5 a0: one extra point
C     (Rrel, C_k/Rrel^6) is appended to the table before building the
C     spline (same "smooth connection" procedure as the paper's
C     Appendix C for V^(1)_2; applied here to V^(1)_0 as well, per the
C     plan in CONTEXT_Tm_CPL9.md section 7).
C   - R > Rrel: pure asymptotic form V(R) = C^(1)_k / R^6.
C   - R < 7.0 a0 (below the table, NOT covered by any ab initio data):
C     V^(1)_0(R) -- dominant, genuinely repulsive-wall behaviour --
C     is extrapolated as a single exponential A*exp(-b*(R-7.0))
C     matched in BOTH value and derivative to the spline at R=7.0
C     (a smooth, C1-continuous steep repulsive wall).
C     V^(1)_2(R) -- an anisotropy STRENGTH, not a repulsive pair
C     potential; its slope at R=7.0 does not point in a direction that
C     would give a sensible exponential wall (see chat -- forcing one
C     gives the wrong sign).  Since V^(1)_0's exponential wall already
C     makes this region completely inaccessible at ultracold collision
C     energies, V^(1)_2(R) is simply held CONSTANT at its R=7.0 value
C     for R<7.0.  This is a deliberate, safe simplification: the exact
C     short-range shape of V^(1)_2 has no measurable effect on
C     scattering observables here.
C
C  *** UNITS: R ARRIVES IN ANGSTROM (CONFIRMED FROM THE MOLSCAT MANUAL
C  -- Herman checked, 2026-07). Table 5 and the C^(1)_k coefficients
C  are natively in a0 (bohr) / cm^-1, so ALL internal computation below
C  (V0FULL, V2FULL, DV0FULL, DV2FULL, the PCHIP spline, RTAB, RMIN,
C  RREL, C0, C2) is deliberately left working in bohr, UNCHANGED from
C  the version already validated in Stage C (test_pot_tm2.f, all PASS)
C  -- only touching this arithmetic again would risk re-introducing a
C  bug in code that is already known to be correct.  Instead, the
C  Angstrom<->bohr conversion is applied ONLY at the two points where
C  this routine talks to the outside world: R is divided by
C  bohr_to_Angstrom on entry to VSTAR/VSTAR1, and the VSTAR1 derivative
C  is divided by bohr_to_Angstrom once more for the chain rule
C  (dV/dR_Angstrom = dV/dR_bohr / bohr_to_Angstrom).  V itself (an
C  energy, cm^-1) needs no conversion either way.  RM is left at 1.D0
C  (see note below) so this conversion is explicit and not implicitly
C  tangled up with however potenl.f may or may not pre-scale R by RM.
C  RM and EPSIL are both set to 1.D0 below: EPSIL=1 because V is
C  already returned directly in cm^-1; RM=1 because the Angstrom
C  conversion above is now handled explicitly in this file instead.
C  RECOMMENDED SANITY CHECK once molscat-Tm2 is built: print V(R) near
C  the potential minimum -- it should be about -825.6 cm^-1 and fall
C  near R=8.5*bohr_to_Angstrom=4.50 Angstrom (NOT near R=8.5 Angstrom,
C  which would indicate the conversion direction got flipped).
C  =====================================================================
      USE physical_constants, ONLY: hartree_in_inv_cm, bohr_to_Angstrom
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      LOGICAL LFIRST
      DIMENSION RTAB(28),V0TAB(28),V2TAB(28)
      DIMENSION REXT(29),F0EXT(29),F2EXT(29),S0EXT(29),S2EXT(29)
      DATA LFIRST/.TRUE./
C
      DATA RTAB / 7.0D0, 7.2D0, 7.4D0, 7.6D0, 7.8D0, 8.0D0, 8.1D0,
     1  8.2D0, 8.3D0, 8.4D0, 8.5D0, 8.6D0, 8.7D0, 8.8D0, 8.9D0, 9.0D0,
     2  9.1D0, 9.2D0, 9.3D0, 9.4D0, 9.5D0, 9.6D0, 9.8D0,10.0D0,10.2D0,
     3 10.4D0,11.0D0,12.0D0 /
      DATA V0TAB /
     1   594.90655204D0,  133.54032587D0, -206.86711784D0,
     2  -450.80090236D0, -618.71806662D0, -727.59730629D0,
     3  -764.32796312D0, -791.25843938D0, -809.69858744D0,
     4  -820.81468219D0, -825.64585040D0, -825.11351391D0,
     5  -820.03654142D0, -811.13357004D0, -799.04967438D0,
     6  -784.32927135D0, -767.46932713D0, -748.89238132D0,
     7  -728.98072835D0, -708.05036358D0, -686.38223486D0,
     8  -664.21695040D0, -619.19702752D0, -574.30698185D0,
     9  -530.48373364D0, -488.37811753D0, -375.85614024D0,
     A  -238.11360714D0 /
      DATA V2TAB /
     1  2.2468153D0, 3.0370500D0, 3.5086367D0, 3.7301851D0,
     2  3.7584559D0, 3.6512236D0, 3.5614030D0, 3.4539311D0,
     3  3.3328684D0, 3.2015754D0, 3.0630180D0, 2.9197678D0,
     4  2.7742074D0, 2.6279605D0, 2.4824668D0, 2.3390504D0,
     5  2.1986898D0, 2.0622174D0, 1.9302292D0, 1.8031947D0,
     6  1.6814415D0, 1.5651986D0, 1.3496968D0, 1.1570413D0,
     7  0.98642598D0, 0.83669271D0, 0.49635265D0, 0.18838511D0 /
C
      RM=1.D0
      EPSIL=1.D0
      IF (I.LT.1 .OR. I.GT.2) THEN
        WRITE(6,*) ' *** VINIT (pot-Tm2): I OUT OF RANGE, I =',I
        STOP
      ENDIF
      IF (.NOT.LFIRST) RETURN
      LFIRST=.FALSE.
C
C  LONG-RANGE VAN DER WAALS COEFFICIENTS, TABLE 1 (i=1 FAMILY ONLY),
C  CONVERTED FROM Eh*a0^6 TO cm^-1*a0^6 USING hartree_in_inv_cm.
      AUCM=hartree_in_inv_cm
      C0=-1672.11503064954D0*AUCM
      C2= 0.788488761D0*AUCM
      RMIN=RTAB(1)
      RREL=RTAB(28)+0.5D0
C
C  BUILD EXTENDED KNOT SET: TABLE + ONE POINT AT RREL WHERE R^6*V(R)
C  IS FORCED TO THE ASYMPTOTIC C_k VALUE (SAME PROCEDURE AS APPENDIX C
C  OF THE PAPER FOR V^(1)_2, APPLIED HERE TO BOTH V^(1)_0 AND V^(1)_2).
      DO K=1,28
        REXT(K)=RTAB(K)
        F0EXT(K)=RTAB(K)**6*V0TAB(K)
        F2EXT(K)=RTAB(K)**6*V2TAB(K)
      ENDDO
      REXT(29)=RREL
      F0EXT(29)=C0
      F2EXT(29)=C2
C
      CALL PCHIPSET(REXT,F0EXT,29,S0EXT)
      CALL PCHIPSET(REXT,F2EXT,29,S2EXT)
C
C  SHORT-RANGE (R<RMIN) PARAMETERS.
C  V0: EXPONENTIAL WALL A*EXP(-B0*(R-RMIN)), MATCHED IN VALUE AND
C  DERIVATIVE TO THE SPLINE AT RMIN.
      F0R=PCHIPEV(REXT,F0EXT,S0EXT,29,RMIN)
      DF0R=PCHIPDV(REXT,F0EXT,S0EXT,29,RMIN)
      V0MIN=F0R/RMIN**6
      DV0MIN=(DF0R-6.D0*F0R/RMIN)/RMIN**6
      B0=-DV0MIN/V0MIN
C  V2: HELD CONSTANT BELOW RMIN (SEE HEADER COMMENT -- THE SIGN OF THE
C  LOCAL SLOPE THERE DOES NOT SUPPORT A PHYSICALLY SENSIBLE WALL, AND
C  V0'S WALL ALREADY MAKES THIS REGION INACCESSIBLE).
      V2MIN=PCHIPEV(REXT,F2EXT,S2EXT,29,RMIN)/RMIN**6
C
      WRITE(6,600)
  600 FORMAT(/'  pot-Tm2: Tm2 radial strengths V^(1)_0(R), V^(1)_2(R)'
     1       ,', Tiesinga et al. NJP 23, 085007 (2021), Table 5.')
      WRITE(6,601) RTAB(1),RTAB(28)
  601 FORMAT('  Table 5 covers R =',F6.2,' to',F6.2,' a0 (PCHIP ',
     1       'spline of R^6*V(R)).')
      WRITE(6,602) RREL,C0,C2
  602 FORMAT('  R >',F6.2,' a0: asymptotic C_k/R^6 tail, ',
     1       'C0 =',1PE14.6,' C2 =',E14.6,' cm^-1 a0^6.')
      WRITE(6,603) RMIN,B0,V2MIN
  603 FORMAT('  R <',F6.2,' a0: V0 exponential wall (rate b0 =',
     1       F8.4,' a0^-1), V2 held constant at',F10.6,' cm^-1.')
      RETURN
C
C============================================================ VSTAR ===
      ENTRY VSTAR(I,R,V)
      RBOHR=R/bohr_to_Angstrom
      IF (I.EQ.1) THEN
        V=V0FULL(RBOHR,RMIN,RREL,C0,V0MIN,B0,REXT,F0EXT,S0EXT)
      ELSEIF (I.EQ.2) THEN
        V=V2FULL(RBOHR,RMIN,RREL,C2,V2MIN,REXT,F2EXT,S2EXT)
      ELSE
        WRITE(6,*) ' *** VSTAR (pot-Tm2): I OUT OF RANGE, I =',I
        STOP
      ENDIF
      RETURN
C
C=========================================================== VSTAR1 ===
      ENTRY VSTAR1(I,R,V)
      RBOHR=R/bohr_to_Angstrom
      IF (I.EQ.1) THEN
        V=DV0FULL(RBOHR,RMIN,RREL,C0,V0MIN,B0,REXT,F0EXT,S0EXT)
      ELSEIF (I.EQ.2) THEN
        V=DV2FULL(RBOHR,RMIN,RREL,C2,REXT,F2EXT,S2EXT)
      ELSE
        WRITE(6,*) ' *** VSTAR1 (pot-Tm2): I OUT OF RANGE, I =',I
        STOP
      ENDIF
C  CHAIN RULE: dV/dR_Angstrom = dV/dR_bohr * dR_bohr/dR_Angstrom
C            = dV/dR_bohr / bohr_to_Angstrom
      V=V/bohr_to_Angstrom
      RETURN
C
C=========================================================== VSTAR2 ===
      ENTRY VSTAR2(I,R,V)
      WRITE(6,*) ' *** VSTAR2 (pot-Tm2): SECOND DERIVATIVES NOT ',
     1           'IMPLEMENTED'
      STOP
      END
C ======================================================= END OF VINIT
      FUNCTION V0FULL(R,RMIN,RREL,C0,V0MIN,B0,REXT,F0EXT,S0EXT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION REXT(29),F0EXT(29),S0EXT(29)
      IF (R.LT.RMIN) THEN
        V0FULL=V0MIN*EXP(-B0*(R-RMIN))
      ELSEIF (R.GT.RREL) THEN
        V0FULL=C0/R**6
      ELSE
        V0FULL=PCHIPEV(REXT,F0EXT,S0EXT,29,R)/R**6
      ENDIF
      RETURN
      END
C ====================================================== END OF V0FULL
      FUNCTION V2FULL(R,RMIN,RREL,C2,V2MIN,REXT,F2EXT,S2EXT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION REXT(29),F2EXT(29),S2EXT(29)
      IF (R.LT.RMIN) THEN
        V2FULL=V2MIN
      ELSEIF (R.GT.RREL) THEN
        V2FULL=C2/R**6
      ELSE
        V2FULL=PCHIPEV(REXT,F2EXT,S2EXT,29,R)/R**6
      ENDIF
      RETURN
      END
C ====================================================== END OF V2FULL
      FUNCTION DV0FULL(R,RMIN,RREL,C0,V0MIN,B0,REXT,F0EXT,S0EXT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION REXT(29),F0EXT(29),S0EXT(29)
      IF (R.LT.RMIN) THEN
        DV0FULL=-B0*V0MIN*EXP(-B0*(R-RMIN))
      ELSEIF (R.GT.RREL) THEN
        DV0FULL=-6.D0*C0/R**7
      ELSE
        FR=PCHIPEV(REXT,F0EXT,S0EXT,29,R)
        DFR=PCHIPDV(REXT,F0EXT,S0EXT,29,R)
        DV0FULL=(DFR-6.D0*FR/R)/R**6
      ENDIF
      RETURN
      END
C ===================================================== END OF DV0FULL
      FUNCTION DV2FULL(R,RMIN,RREL,C2,REXT,F2EXT,S2EXT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION REXT(29),F2EXT(29),S2EXT(29)
      IF (R.LT.RMIN) THEN
        DV2FULL=0.D0
      ELSEIF (R.GT.RREL) THEN
        DV2FULL=-6.D0*C2/R**7
      ELSE
        FR=PCHIPEV(REXT,F2EXT,S2EXT,29,R)
        DFR=PCHIPDV(REXT,F2EXT,S2EXT,29,R)
        DV2FULL=(DFR-6.D0*FR/R)/R**6
      ENDIF
      RETURN
      END
C ===================================================== END OF DV2FULL
      SUBROUTINE PCHIPSET(X,Y,N,S)
C  =====================================================================
C  Fritsch-Carlson monotone cubic Hermite (PCHIP) tangent setup.
C  Given N knots X(1:N) (strictly increasing) and values Y(1:N),
C  returns tangents S(1:N) such that the resulting piecewise cubic
C  Hermite interpolant is monotone on every interval where the data
C  themselves are monotone (no overshoot).  Standard algorithm
C  (Fritsch & Carlson, SIAM J. Numer. Anal. 17, 238 (1980));
C  cross-checked against scipy.interpolate.PchipInterpolator to
C  machine precision on random test data (see pchip_check.py).
C  =====================================================================
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N),Y(N),S(N),H(N-1),D(N-1)
C
      DO K=1,N-1
        H(K)=X(K+1)-X(K)
        D(K)=(Y(K+1)-Y(K))/H(K)
      ENDDO
C
C  INTERIOR TANGENTS
      DO K=2,N-1
        D0=D(K-1)
        D1=D(K)
        IF (D0*D1.LE.0.D0) THEN
          S(K)=0.D0
        ELSE
          W1=2.D0*H(K)+H(K-1)
          W2=H(K)+2.D0*H(K-1)
          S(K)=(W1+W2)/(W1/D0+W2/D1)
        ENDIF
      ENDDO
C
C  LEFT-END TANGENT (SHAPE-PRESERVING 3-POINT ONE-SIDED FORMULA)
      AM=((2.D0*H(1)+H(2))*D(1)-H(1)*D(2))/(H(1)+H(2))
      IF (SIGN(1.D0,AM).NE.SIGN(1.D0,D(1)) .OR. D(1).EQ.0.D0) THEN
        AM=0.D0
      ELSEIF (SIGN(1.D0,D(1)).NE.SIGN(1.D0,D(2))
     1        .AND. ABS(AM).GT.ABS(3.D0*D(1))) THEN
        AM=3.D0*D(1)
      ENDIF
      S(1)=AM
C
C  RIGHT-END TANGENT
      AM=((2.D0*H(N-1)+H(N-2))*D(N-1)-H(N-1)*D(N-2))/(H(N-1)+H(N-2))
      IF (SIGN(1.D0,AM).NE.SIGN(1.D0,D(N-1)) .OR. D(N-1).EQ.0.D0) THEN
        AM=0.D0
      ELSEIF (SIGN(1.D0,D(N-1)).NE.SIGN(1.D0,D(N-2))
     1        .AND. ABS(AM).GT.ABS(3.D0*D(N-1))) THEN
        AM=3.D0*D(N-1)
      ENDIF
      S(N)=AM
C
      RETURN
      END
C ==================================================== END OF PCHIPSET
      FUNCTION PCHIPEV(X,Y,S,N,XQ)
C  Evaluate the PCHIP interpolant (knots/tangents from PCHIPSET) at XQ.
C  XQ outside [X(1),X(N)] is clipped to the nearest end interval
C  (linear search over N-1 intervals; N is small (29) here so this is
C  not a performance concern).
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N),Y(N),S(N)
      IF (XQ.LE.X(1)) THEN
        K=1
      ELSEIF (XQ.GE.X(N)) THEN
        K=N-1
      ELSE
        K=1
        DO WHILE (X(K+1).LT.XQ)
          K=K+1
        ENDDO
      ENDIF
      H=X(K+1)-X(K)
      T=(XQ-X(K))/H
      H00=2.D0*T**3-3.D0*T**2+1.D0
      H10=T**3-2.D0*T**2+T
      H01=-2.D0*T**3+3.D0*T**2
      H11=T**3-T**2
      PCHIPEV=H00*Y(K)+H10*H*S(K)+H01*Y(K+1)+H11*H*S(K+1)
      RETURN
      END
C ===================================================== END OF PCHIPEV
      FUNCTION PCHIPDV(X,Y,S,N,XQ)
C  Derivative (d/dXQ) of the PCHIP interpolant at XQ.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N),Y(N),S(N)
      IF (XQ.LE.X(1)) THEN
        K=1
      ELSEIF (XQ.GE.X(N)) THEN
        K=N-1
      ELSE
        K=1
        DO WHILE (X(K+1).LT.XQ)
          K=K+1
        ENDDO
      ENDIF
      H=X(K+1)-X(K)
      T=(XQ-X(K))/H
      DH00=6.D0*T**2-6.D0*T
      DH10=3.D0*T**2-4.D0*T+1.D0
      DH01=-6.D0*T**2+6.D0*T
      DH11=3.D0*T**2-2.D0*T
      PCHIPDV=(DH00*Y(K)+DH10*H*S(K)+DH01*Y(K+1)+DH11*H*S(K+1))/H
      RETURN
      END
C ===================================================== END OF PCHIPDV
