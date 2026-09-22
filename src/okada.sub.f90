  recursive subroutine  SPOINT(ALP,X,Y,D,SD,CD,DISL1,DISL2,DISL3, &
  & U1,U2,U3,U11,U12,U21,U22,U31,U32)
!*****
!*****    SURFACE DISPLACEMENT,STRAIN,TILT DUE TO BURIED POINT SOURCE
!*****    IN A SEMIINFINITE MEDIUM     CODED BY  Y.OKADA ... JAN 1985
!*****
!***** INPUT
!*****   ALP   : MEDIUM CONSTANT  MYU/(LAMDA+MYU)
!*****   X,Y   : COORDINATE OF STATION
!*****   D     : SOURCE DEPTH
!*****   SD,CD : SIN,COS OF DIP-ANGLE
!*****          (CD=0.D0, SD=+/-1.D0 SHOULD BE GIVEN FOR VERTICAL FAULT)
!*****   DISL1,DISL2,DISL3 : STRIKE-, DIP- AND TENSILE-DISLOCATION
!*****
!***** OUTPUT
!*****   U1, U2, U3      : DISPLACEMENT ( UNIT= UNIT OF DISL / AREA )
!*****   U11,U12,U21,U22 : STRAIN       ( UNIT= UNIT OF DISL /
!*****   U31,U32         : TILT                 UNIT OF X,Y,D /AREA )
!*****
  use, intrinsic :: iso_fortran_env, only: real64
  implicit none
  real(real64), intent(in) :: ALP, X, Y, D, SD, CD, DISL1, DISL2, DISL3
  real(real64), intent(out) :: U1, U2, U3, U11, U12, U21, U22, U31, U32
  real(real64), parameter :: F0 = 0.D0
  real(real64), parameter :: F1 = 1.D0
  real(real64), parameter :: F2 = 2.D0
  real(real64), parameter :: F3 = 3.D0
  real(real64), parameter :: F4 = 4.D0
  real(real64), parameter :: F5 = 5.D0
  real(real64), parameter :: F8 = 8.D0
  real(real64), parameter :: F9 = 9.D0
  real(real64) :: A1, A2, A3, A4, A5, B1, B2, B3, B4, C1
  real(real64) :: C2, C3, D2, DR, FQ, FS, FX, P, PI2, Q
  real(real64) :: QR, QRP, QRQ, QRX, R, R12, R2, R3, R32, R33
  real(real64) :: R5, R53, R54, RD, S, SDCD, SDSD, UN, X2, XR
  real(real64) :: XY, XYR, Y2, YR
  PI2=6.283185307179586D0
!*****
  P =Y*CD + D*SD
  Q =Y*SD - D*CD
  S =P*SD + Q*CD
  X2=X*X
  Y2=Y*Y
  XY=X*Y
  D2=D*D
  R2=X2 + Y2 + D2
  R =SQRT(R2)
  R3=R *R2
  R5=R3*R2
  QR=F3*Q/R5
  XR =F5*X2/R2
  YR =F5*Y2/R2
  XYR=F5*XY/R2
  DR =F5*D /R2
  RD =R + D
  R12=F1/(R*RD*RD)
  R32=R12*(F2*R + D)/ R2
  R33=R12*(F3*R + D)/(R2*RD)
  R53=R12*(F8*R2 + F9*R*D + F3*D2)/(R2*R2*RD)
  R54=R12*(F5*R2 + F4*R*D +    D2)/R3*R12
!*****
  A1= ALP*Y*(R12-X2*R33)
  A2= ALP*X*(R12-Y2*R33)
  A3= ALP*X/R3 - A2
  A4=-ALP*XY*R32
  A5= ALP*( F1/(R*RD) - X2*R32 )
  B1= ALP*(-F3*XY*R33      + F3*X2*XY*R54)
  B2= ALP*( F1/R3 - F3*R12 + F3*X2*Y2*R54)
  B3= ALP*( F1/R3 - F3*X2/R5) - B2
  B4=-ALP*F3*XY/R5 - B1
  C1=-ALP*Y*(R32 - X2*R53)
  C2=-ALP*X*(R32 - Y2*R53)
  C3=-ALP*F3*X*D/R5 - C2
!*****
  U1 =F0
  U2 =F0
  U3 =F0
  U11=F0
  U12=F0
  U21=F0
  U22=F0
  U31=F0
  U32=F0
!**************************************
!*****                            *****
!*****  STRIKE-SLIP CONTRIBUTION  *****
!*****                            *****
!**************************************
  IF(DISL1.EQ.F0)  GO TO 200
  UN=DISL1/PI2
  QRX=QR*X
  FX=F3*X/R5*SD
  U1 =U1 - UN*( QRX*X + A1*SD )
  U2 =U2 - UN*( QRX*Y + A2*SD )
  U3 =U3 - UN*( QRX*D + A4*SD )
  U11=U11- UN*( QRX* (F2-XR)        + B1*SD )
  U12=U12- UN*(-QRX*XYR      + FX*X + B2*SD )
  U21=U21- UN*( QR*Y*(F1-XR)        + B2*SD )
  U22=U22- UN*( QRX *(F1-YR) + FX*Y + B4*SD )
  U31=U31- UN*( QR*D*(F1-XR)        + C1*SD )
  U32=U32- UN*(-QRX*DR*Y     + FX*D + C2*SD )
!**************************************
!*****                            *****
!*****    DIP-SLIP CONTRIBUTION   *****
!*****                            *****
!**************************************
  200 IF(DISL2.EQ.F0)  GO TO 300
  UN=DISL2/PI2
  SDCD=SD*CD
  QRP=QR*P
  FS=F3*S/R5
  U1 =U1 - UN*( QRP*X - A3*SDCD )
  U2 =U2 - UN*( QRP*Y - A1*SDCD )
  U3 =U3 - UN*( QRP*D - A5*SDCD )
  U11=U11- UN*( QRP*(F1-XR)        - B3*SDCD )
  U12=U12- UN*(-QRP*XYR     + FS*X - B1*SDCD )
  U21=U21- UN*(-QRP*XYR            - B1*SDCD )
  U22=U22- UN*( QRP*(F1-YR) + FS*Y - B2*SDCD )
  U31=U31- UN*(-QRP*DR*X           - C3*SDCD )
  U32=U32- UN*(-QRP*DR*Y    + FS*D - C1*SDCD )
!****************************************
!*****                              *****
!*****  TENSILE-FAULT CONTRIBUTION  *****
!*****                              *****
!****************************************
  300 IF(DISL3.EQ.F0)  GO TO 900
  UN=DISL3/PI2
  SDSD=SD*SD
  QRQ=QR*Q
  FQ=F2*QR*SD
  U1 =U1 + UN*( QRQ*X - A3*SDSD )
  U2 =U2 + UN*( QRQ*Y - A1*SDSD )
  U3 =U3 + UN*( QRQ*D - A5*SDSD )
  U11=U11+ UN*( QRQ*(F1-XR)        - B3*SDSD )
  U12=U12+ UN*(-QRQ*XYR     + FQ*X - B1*SDSD )
  U21=U21+ UN*(-QRQ*XYR            - B1*SDSD )
  U22=U22+ UN*( QRQ*(F1-YR) + FQ*Y - B2*SDSD )
  U31=U31+ UN*(-QRQ*DR*X           - C3*SDSD )
  U32=U32+ UN*(-QRQ*DR*Y    + FQ*D - C1*SDSD )
!*****
  900 RETURN
  end subroutine SPOINT
  recursive subroutine  SRECTF(ALP,X,Y,DEP,AL,AW,SD,CD,DISL1,DISL2,DISL3, &
  & U1,U2,U3,U11,U12,U21,U22,U31,U32)
!*****
!*****   SURFACE DISPLACEMENTS,STRAINS AND TILTS DUE TO RECTANGULAR
!*****   FAULT IN A HALF-SPACE       CODED BY  Y.OKADA ... JAN 1985
!*****
!***** INPUT
!*****   ALP   : MEDIUM CONSTANT  MYU/(LAMDA+MYU)
!*****   X,Y   : COORDINATE OF STATION
!*****   DEP   : SOURCE DEPTH
!*****   AL,AW : LENGTH AND WIDTH OF FAULT
!*****   SD,CD : SIN,COS OF DIP-ANGLE
!*****          (CD=0.D0, SD=+/-1.D0 SHOULD BE GIVEN FOR VERTICAL FAULT)
!*****   DISL1,DISL2,DISL3 : STRIKE-, DIP- AND TENSILE-DISLOCATION
!*****
!***** OUTPUT
!*****   U1, U2, U3      : DISPLACEMENT ( UNIT= UNIT OF DISL     )
!*****   U11,U12,U21,U22 : STRAIN       ( UNIT= UNIT OF DISL /
!*****   U31,U32         : TILT                 UNIT OF X,Y,,,AW )
!*****
!***** SUBROUTINE USED...SRECTG
!*****
  use, intrinsic :: iso_fortran_env, only: real64
  implicit none
  real(real64), intent(in) :: ALP, X, Y, DEP, AL, AW, SD, CD, DISL1, DISL2, DISL3
  real(real64), intent(out) :: U1, U2, U3, U11, U12, U21, U22, U31, U32
  real(real64), parameter :: F0 = 0.D0
  real(real64), parameter :: F1 = 1.D0
  integer :: I, J, JK, K
  real(real64) :: DU(9), ET, P, Q, SIGN, U(9), XI
!*****
  P = Y*CD + DEP*SD
  Q = Y*SD - DEP*CD
!*****
  U = F0
!*****
  do K=1,2
  IF(K.EQ.1)  ET=P
  IF(K.EQ.2)  ET=P-AW
  do J=1,2
  IF(J.EQ.1)  XI=X
  IF(J.EQ.2)  XI=X-AL
  JK=J+K
  IF(JK.NE.3)  SIGN= F1
  IF(JK.EQ.3)  SIGN=-F1
  CALL SRECTG(ALP,XI,ET,Q,SD,CD,DISL1,DISL2,DISL3, &
  & DU(1),DU(2),DU(3),DU(4),DU(5),DU(6),DU(7),DU(8),DU(9))
  do I=1,9
  U(I)=U(I)+SIGN*DU(I)
  end do
  end do
  end do
  U1 =U(1)
  U2 =U(2)
  U3 =U(3)
  U11=U(4)
  U12=U(5)
  U21=U(6)
  U22=U(7)
  U31=U(8)
  U32=U(9)
  RETURN
  end subroutine SRECTF
  recursive subroutine  SRECTG(ALP,XI,ET,Q,SD,CD,DISL1,DISL2,DISL3, &
  & U1,U2,U3,U11,U12,U21,U22,U31,U32)
!*****
!*****   INDEFINITE INTEGRAL OF SURFACE DISPLACEMENTS, STRAINS AND TILTS
!*****   DUE TO FINITE FAULT IN A SEMIINFINITE MEDIUM
!*****                                    CODED BY  Y.OKADA ... JAN 1985
!***** INPUT
!*****   ALP     : MEDIUM CONSTANT  MYU/(LAMDA+MYU)
!*****   XI,ET,Q : FAULT COORDINATE
!*****   SD,CD   : SIN,COS OF DIP-ANGLE
!*****          (CD=0.D0, SD=+/-1.D0 SHOULD BE GIVEN FOR VERTICAL FAULT)
!*****   DISL1,DISL2,DISL3 : STRIKE-, DIP- AND TENSILE-DISLOCATION
!*****
!***** OUTPUT
!*****   U1, U2, U3      : DISPLACEMENT ( UNIT= UNIT OF DISL    )
!*****   U11,U12,U21,U22 : STRAIN       ( UNIT= UNIT OF DISL /
!*****   U31,U32         : TILT                 UNIT OF XI,ET,Q )
!*****
  use, intrinsic :: iso_fortran_env, only: real64
  implicit none
  real(real64), intent(in) :: ALP, XI, ET, Q, SD, CD, DISL1, DISL2, DISL3
  real(real64), intent(out) :: U1, U2, U3, U11, U12, U21, U22, U31, U32
  real(real64), parameter :: F0 = 0.D0
  real(real64), parameter :: F1 = 1.D0
  real(real64), parameter :: F2 = 2.D0
  real(real64) :: A1, A2, A3, A4, A5, AET, AXI, B1, B2, B3
  real(real64) :: B4, C1, C2, C3, D, DLE, ET2, PI2, Q2, R
  real(real64) :: R2, R3, RD, RD2, RE, REQ, RET, RRD, RRE, RRX
  real(real64) :: SDCD, SDSD, TD, TT, UN, X, XI2, Y
  PI2=6.283185307179586D0
!*****
  XI2=XI*XI
  ET2=ET*ET
  Q2=Q*Q
  R2=XI2+ET2+Q2
  R =DSQRT(R2)
  R3=R*R2
  D =ET*SD-Q*CD
  Y =ET*CD+Q*SD
  RET=R+ET
  IF(RET.LT.F0)  RET=F0
  RD =R+D
  RRD=F1/(R*RD)
!*****
  IF( Q .NE.F0)  TT = DATAN( XI*ET/(Q*R) )
  IF( Q .EQ.F0)  TT = F0
  IF(RET.NE.F0)  RE = F1/RET
  IF(RET.EQ.F0)  RE = F0
  IF(RET.NE.F0)  DLE= DLOG(RET)
  IF(RET.EQ.F0)  DLE=-DLOG(R-ET)
! === Modification to prevent zero-division. ===================================
!     RRX=F1/(R*(R+XI))
  IF(ABS(R*(R+XI)).LT.1.0D-6) THEN
  RRX=1.0D6
  ELSE
  RRX=F1/(R*(R+XI))
  ENDIF
! ==============================================================================
  RRE=RE/R
  AXI=(F2*R+XI)*RRX*RRX/R
  AET=(F2*R+ET)*RRE*RRE/R
  IF(CD.EQ.F0)  GO TO 20
!*****
!***** INCLINED FAULT
!*****
  TD=SD/CD
  X =DSQRT(XI2+Q2)
  IF(XI.EQ.F0)  A5=F0
  IF(XI.NE.F0) &
  & A5= ALP*F2/CD*DATAN( (ET*(X+Q*CD)+X*(R+X)*SD) / (XI*(R+X)*CD) )
  A4= ALP/CD*( DLOG(RD) - SD*DLE )
  A3= ALP*(Y/RD/CD - DLE) + TD*A4
  A1=-ALP/CD*XI/RD        - TD*A5
  C1= ALP/CD*XI*(RRD - SD*RRE)
  C3= ALP/CD*(Q*RRE - Y*RRD)
  B1= ALP/CD*(XI2*RRD - F1)/RD - TD*C3
  B2= ALP/CD*XI*Y*RRD/RD       - TD*C1
  GO TO 30
!*****
!***** VERTICAL FAULT
!*****
  20 RD2=RD*RD
  A1=-ALP/F2*XI*Q/RD2
  A3= ALP/F2*( ET/RD + Y*Q/RD2 - DLE )
  A4=-ALP*Q/RD
  A5=-ALP*XI*SD/RD
  B1= ALP/F2*  Q  /RD2*(F2*XI2*RRD - F1)
  B2= ALP/F2*XI*SD/RD2*(F2*Q2 *RRD - F1)
  C1= ALP*XI*Q*RRD/RD
  C3= ALP*SD/RD*(XI2*RRD - F1)
!*****
  30 A2=-ALP*DLE - A3
  B3=-ALP*XI*RRE - B2
  B4=-ALP*( CD/R + Q*SD*RRE ) - B1
  C2= ALP*(-SD/R + Q*CD*RRE ) - C3
!*****
  U1 =F0
  U2 =F0
  U3 =F0
  U11=F0
  U12=F0
  U21=F0
  U22=F0
  U31=F0
  U32=F0
!**************************************
!*****                            *****
!*****  STRIKE-SLIP CONTRIBUTION  *****
!*****                            *****
!**************************************
  IF(DISL1.EQ.F0)  GO TO 200
  UN=DISL1/PI2
  REQ=RRE*Q
  U1 =U1 - UN*( REQ*XI +   TT    + A1*SD )
  U2 =U2 - UN*( REQ*Y  + Q*CD*RE + A2*SD )
  U3 =U3 - UN*( REQ*D  + Q*SD*RE + A4*SD )
  U11=U11+ UN*( XI2*Q*AET - B1*SD )
  U12=U12+ UN*( XI2*XI*( D/(ET2+Q2)/R3 - AET*SD ) - B2*SD )
  U21=U21+ UN*( XI*Q/R3*CD + (XI*Q2*AET - B2)*SD )
  U22=U22+ UN*( Y *Q/R3*CD + (Q*SD*(Q2*AET-F2*RRE) &
  & -(XI2+ET2)/R3*CD - B4)*SD )
  U31=U31+ UN*(-XI*Q2*AET*CD + (XI*Q/R3 - C1)*SD )
  U32=U32+ UN*( D*Q/R3*CD + (XI2*Q*AET*CD - SD/R + Y*Q/R3 - C2)*SD )
!**************************************
!*****                            *****
!*****    DIP-SLIP CONTRIBUTION   *****
!*****                            *****
!**************************************
  200 IF(DISL2.EQ.F0)  GO TO 300
  UN=DISL2/PI2
  SDCD=SD*CD
  U1 =U1 - UN*( Q/R             - A3*SDCD )
  U2 =U2 - UN*( Y*Q*RRX + CD*TT - A1*SDCD )
  U3 =U3 - UN*( D*Q*RRX + SD*TT - A5*SDCD )
  U11=U11+ UN*( XI*Q/R3            + B3*SDCD )
  U12=U12+ UN*( Y *Q/R3 - SD/R     + B1*SDCD )
  U21=U21+ UN*( Y *Q/R3 + Q*CD*RRE + B1*SDCD )
  U22=U22+ UN*( Y*Y*Q*AXI - (F2*Y*RRX + XI*CD*RRE)*SD + B2*SDCD )
  U31=U31+ UN*( D *Q/R3 + Q*SD*RRE + C3*SDCD )
  U32=U32+ UN*( Y*D*Q*AXI - (F2*D*RRX + XI*SD*RRE)*SD + C1*SDCD )
!****************************************
!*****                              *****
!*****  TENSILE-FAULT CONTRIBUTION  *****
!*****                              *****
!****************************************
  300 IF(DISL3.EQ.F0)  GO TO 900
  UN=DISL3/PI2
  SDSD=SD*SD
  U1 =U1 + UN*( Q2*RRE                       - A3*SDSD )
  U2 =U2 + UN*(-D*Q*RRX - SD*(XI*Q*RRE - TT) - A1*SDSD )
  U3 =U3 + UN*( Y*Q*RRX + CD*(XI*Q*RRE - TT) - A5*SDSD )
  U11=U11- UN*( XI*Q2*AET             + B3*SDSD )
  U12=U12- UN*(-D*Q/R3 - XI2*Q*AET*SD + B1*SDSD )
  U21=U21- UN*( Q2*(CD/R3 + Q*AET*SD) + B1*SDSD )
  U22=U22- UN*((Y*CD-D*SD)*Q2*AXI - F2*Q*SD*CD*RRX &
  & - (XI*Q2*AET - B2)*SDSD )
  U31=U31- UN*( Q2*(SD/R3 - Q*AET*CD) + C3*SDSD )
  U32=U32- UN*((Y*SD+D*CD)*Q2*AXI + XI*Q2*AET*SD*CD &
  & - (F2*Q*RRX - C1)*SDSD )
!*****
  900 RETURN
  end subroutine SRECTG
