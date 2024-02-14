      PARAMETER (NX=1244)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION E(0:NX),T(0:NX),F(0:NX),C(0:NX),S(0:NX),R1(0:NX)
      DIMENSION WC(0:NX),FCC1(0:NX),FSS1(0:NX),FCS1(0:NX),E0(0:NX)
      DIMENSION EC(0:NX),FES1(0:NX),FEC1(0:NX),GG(0:NX),PPE(0:NX)

      PI=4*ATAN(1D0)
      NX0=750
      NX1=100
      NX2=100
      NXT=1050
      WRITE(6,100)'NX0 NX1 NX2',NX0,NX1,NX2
      READ(5,*) NX0,NX1,NX2

      OPEN (1,FILE='1')
      DO 1 K=1,NX
        READ(1,*) T(K),F(K)
        WRITE (92,130) T(K),F(K)
    1 CONTINUE
      CLOSE (1)
!--------------------------------------------------------------

      G1= 0.1153D0
      DG= 0.0001D0
      PE1=2.9355D0
      DPE=0.0001D0

      WRITE(6,101)'G1 DG PE1 DPE',G1,DG,PE1,DPE
      READ(5,*) G1,DG,PE1,DPE

      DO 2 I=0,NX1
        G=G1+I*DG
        GG(I)=G
        WRITE(11,*) I,DG,G

        DO 3 J=0,NX2
          PE=PE1+J*DPE
          PPE(J)=PE
          W=2*PI/PE

          FCC=0D0
          FSS=0D0
          FCS=0D0
          FES=0D0
          FEC=0D0

          DO 4 K=NX0,NXT
            WRITE(93,130) T(K),F(K)
            S(K)=SIN(W*T(K))
            C(K)=COS(W*T(K))
            E(K)=EXP(G*T(K))

            FCC=FCC+(E(K)*C(K))**2
            FSS=FSS+(E(K)*S(K))**2
            FCS=FCS+E(K)**2*C(K)*S(K)
            FES=FES+E(K)*S(K)*F(K)
            FEC=FEC+E(K)*C(K)*F(K)
    4     CONTINUE

C          stop "4"

          DET0=FCC*FSS-FCS*FCS
          DETA=FEC*FSS-FES*FCS
          DETB=FCC*FES-FCS*FEC
          A=DETA/DET0
          B=DETB/DET0

          R=0.0D0
          AM=0.0D0

          DO 5 K=NX0,NXT
            S(K)=SIN(W*T(K))
            C(K)=COS(W*T(K))
            E(K)=EXP(G*T(K))
            AMM=(A*C(K)+B*S(K))*E(K)
            R=R+(F(K)-AMM)**2
    5     CONTINUE

          R1(J)=R
          AM=AMM/(NXT-NX0)
          WRITE (3,*) G,PE,R
          WRITE (33,*) A,B
          WRITE (*,130) G,PE,A,B,R,AM
          WRITE (91,*) G,PE,A,B


    3     CONTINUE

      RM=R1(0)
      PM=PPE(0)

      DO 6,M=1,NX2
        IF(R1(M).LE.RM) THEN
          RM=R1(M)
          PM=PPE(M)
        ENDIF
    6 CONTINUE

      WRITE(4,*) GG(I),RM,PM

    2     CONTINUE
      CLOSE(4)

      OPEN (1,FILE='fort.4')
      DO 7 K=0,NX1
        READ(4,*) GG(K),R1(K),PPE(K)
    7 CONTINUE


      RM=R1(0)
      PM=PPE(0)
      GM=GG(0)

      DO 8,M=1,NX1
        IF(R1(M).LE.RM) THEN
          RM=R1(M)
          PM=PPE(M)
          GM=GG(M)
        ENDIF
    8 CONTINUE
      WRITE(98,*) GM,RM,PM


        G=GM
        P=PM
        W=2*PI/P
        A=1.8866215086678465E-004
        B= 2.5719126669956010E-004

      DO 9, K=NX0,NXT
        S(K)=SIN(W*T(K))
        C(K)=COS(W*T(K))
        E(K)=EXP(G*T(K))
        AMM=(A*C(K)+B*S(K))*E(K)
        WRITE (90,*) T(K), AMM
    9 CONTINUE



 100  FORMAT(A,3I4)
 101  FORMAT(A,4F12.6)
 120  FORMAT(4X,4F15.6)
 130  FORMAT(6F12.8)
        END
