!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 
!           A planar morphoelastic rod upon an elastic foundation
!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 

 SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP) 
!---------- ---- 

 IMPLICIT NONE
 INTEGER, INTENT(IN) :: NDIM, IJAC, ICP(*)
 DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
 DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
 DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,*), DFDP(NDIM,*)
 DOUBLE PRECISION alpha, gamma, K, L, S, r1, r3, n1, n3, Theta, M 

 ! Define the variables
 S = U(1)
 r1 = U(2)
 r3 = U(3)
 n1 = U(4)
 n3 = U(5)
 Theta = U(6)
 M = U(7)

 ! Define the parameters
 gamma = PAR(1)
 K = PAR(2)
 L = PAR(3)
 alpha = 1.d0 + n3

   ! Define derivatives
   F(1) = L
   F(2) = -L*r3*gamma*M
   F(3) = alpha*gamma*L + L*r1*gamma*M
   F(4) = K*L*(r1 + S*DSIN(Theta)) - n3*gamma*L*M
   F(5) = K*L*alpha*gamma*(r3 - S*DCOS(Theta)) + n1*gamma*L*M
   F(6) = gamma*L*M
   F(7) = -alpha*gamma*L*n1

 END SUBROUTINE FUNC
!----------------------------------------------------------------------

 SUBROUTINE STPNT(NDIM,U,PAR,x) 
!---------- ----- 

 IMPLICIT NONE
 INTEGER, INTENT(IN) :: NDIM
 DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM), PAR(*)
 DOUBLE PRECISION, INTENT(IN) :: x
 DOUBLE PRECISION h, w, L0

 L0 = 0.125
 h = 0.015
 w = 0.01

   PAR(1) = 1.d0 ! Set gamma to one 
   PAR(2) = 0.02  ! Set k
   PAR(3) = 2.d0*(3.d0)**0.5*L0/h   ! Set L    

 END SUBROUTINE STPNT
!----------------------------------------------------------------------

 SUBROUTINE BCND(NDIM,PAR,ICP,NBC,U0,U1,FB,IJAC,DBC) 
!---------- ---- 

 IMPLICIT NONE
 INTEGER, INTENT(IN) :: NDIM, ICP(*), NBC, IJAC
 DOUBLE PRECISION, INTENT(IN) :: PAR(*), U0(NDIM), U1(NDIM)
 DOUBLE PRECISION, INTENT(OUT) :: FB(NBC)
 DOUBLE PRECISION, INTENT(INOUT) :: DBC(NBC,*)
 DOUBLE PRECISION L

 L = PAR(3)

 ! These boundary conditions correspond to the rod being clamped at the boundaries:
 ! S(L) = L, r3(0) = 0, r3(L) = L, r1(0) = r1(L) = 0, Theta(0) = Theta(L) = 0

 FB(1)=U0(1)
 FB(2)=U0(2)
 FB(3)=U1(2)
 FB(4)=U0(3)
 FB(5)=U1(3) - L
 FB(6)=U0(6)
 FB(7)=U1(6)

 END SUBROUTINE BCND

!----------------------------------------------------------------------

 SUBROUTINE PVLS(NDIM,U,PAR)
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: NDIM
 DOUBLE PRECISION, INTENT(IN) :: U(NDIM)
 DOUBLE PRECISION, INTENT(INOUT) :: PAR(*)
 DOUBLE PRECISION, EXTERNAL :: GETP
 DOUBLE PRECISION MAX, MIN, EXTREMUM
 
 MAX = GETP('MAX', 6, U)
 MIN = GETP('MIN', 6, U)

 IF (-MIN > MAX) THEN
   EXTREMUM = -MIN
 ELSE
   EXTREMUM = MAX
 END IF
   
 PAR(4) = EXTREMUM

 END SUBROUTINE PVLS
!----------------------------------------------------------------------
!----------------------------------------------------------------------
 SUBROUTINE ICND
 END SUBROUTINE ICND

 SUBROUTINE FOPT 
 END SUBROUTINE FOPT
!----------------------------------------------------------------------
!----------------------------------------------------------------------