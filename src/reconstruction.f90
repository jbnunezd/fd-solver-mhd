!===============================================================================!
MODULE MOD_Reconstruction
!===============================================================================!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
PRIVATE
!-------------------------------------------------------------------------------!
INTERFACE ReconstructionX
  MODULE PROCEDURE ReconstructionX
END INTERFACE

INTERFACE ReconstructionY
  MODULE PROCEDURE ReconstructionY
END INTERFACE

INTERFACE ReconstructionFixX
  MODULE PROCEDURE ReconstructionFixX
END INTERFACE

INTERFACE ReconstructionFixY
  MODULE PROCEDURE ReconstructionFixY
END INTERFACE
!-------------------------------------------------------------------------------!
PUBLIC :: ReconstructionX
PUBLIC :: ReconstructionY
PUBLIC :: ReconstructionFixX
PUBLIC :: ReconstructionFixY
!-------------------------------------------------------------------------------!
!
!
!
!===============================================================================!
CONTAINS
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE ReconstructionX()
!-------------------------------------------------------------------------------!
USE MOD_FiniteDifference2D_vars,ONLY: V
USE MOD_FiniteDifference2D_vars,ONLY: WM
USE MOD_FiniteDifference2D_vars,ONLY: WP
USE MOD_FiniteDifference2D_vars,ONLY: Fn
USE MOD_FiniteDifference2D_vars,ONLY: Fp
USE MOD_FiniteDifference2D_vars,ONLY: MESH_DX
USE MOD_FiniteDifference2D_vars,ONLY: nVar
USE MOD_FiniteDifference2D_vars,ONLY: nElemsX
USE MOD_FiniteDifference2D_vars,ONLY: nElemsY
USE MOD_FiniteDifference2D_vars,ONLY: nGhosts
USE MOD_FiniteDifference2D_vars,ONLY: Ind
USE MOD_FiniteDifference2D_vars,ONLY: Reconstruction
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! FORMAL ARGUMENTS
!-------------------------------------------------------------------------------!
! LOCAL VARIABLES
!-------------------------------------------------------------------------------!
INTEGER            :: ii, jj
CHARACTER(LEN=255) :: ErrorMessage
!-------------------------------------------------------------------------------!

SELECT CASE (Reconstruction)
  CASE(1)
    DO jj=1,nElemsY
      DO ii=0,nElemsX+1
        WM(1:nVar,ii,jj) = Fn(1:nVar,ii,jj)
        WP(1:nVar,ii,jj) = Fp(1:nVar,ii,jj)
      END DO
    END DO
  CASE(2)
    DO jj=1,nElemsY
      DO ii=0,nElemsX+1
        CALL MUSCL(&
              Fn(1:nVar,-nGhosts+ii:ii+nGhosts,jj),&
              Fp(1:nVar,-nGhosts+ii:ii+nGhosts,jj),&
              WM(1:nVar,ii,jj),&
              WP(1:nVar,ii,jj),&
              Mesh_DX(1))
      END DO
    END DO
  CASE(3,4)
    DO jj=1,nElemsY
      DO ii=0,nElemsX+1
        IF (.NOT. Ind(1,ii,jj)) THEN
          CALL WENO(&
                V(1:nVar,-nGhosts+ii:ii+nGhosts,jj),&
                Fn(1:nVar,-nGhosts+ii:ii+nGhosts,jj),&
                Fp(1:nVar,-nGhosts+ii:ii+nGhosts,jj),&
                WM(1:nVar,ii,jj),&
                WP(1:nVar,ii,jj),&
                Mesh_DX(1),&
                Reconstruction)
        END IF
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Reconstruction not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!-------------------------------------------------------------------------------!
END SUBROUTINE ReconstructionX
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE ReconstructionY()
!-------------------------------------------------------------------------------!
USE MOD_FiniteDifference2D_vars,ONLY: V
USE MOD_FiniteDifference2D_vars,ONLY: WM
USE MOD_FiniteDifference2D_vars,ONLY: WP
USE MOD_FiniteDifference2D_vars,ONLY: Fn
USE MOD_FiniteDifference2D_vars,ONLY: Fp
USE MOD_FiniteDifference2D_vars,ONLY: MESH_DX
USE MOD_FiniteDifference2D_vars,ONLY: nVar
USE MOD_FiniteDifference2D_vars,ONLY: nElemsX
USE MOD_FiniteDifference2D_vars,ONLY: nElemsY
USE MOD_FiniteDifference2D_vars,ONLY: nGhosts
USE MOD_FiniteDifference2D_vars,ONLY: Ind
USE MOD_FiniteDifference2D_vars,ONLY: Reconstruction
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! FORMAL ARGUMENTS
!-------------------------------------------------------------------------------!
! LOCAL VARIABLES
!-------------------------------------------------------------------------------!
INTEGER            :: ii, jj
CHARACTER(LEN=255) :: ErrorMessage
!-------------------------------------------------------------------------------!

SELECT CASE (Reconstruction)
  CASE(1)
    DO jj=0,nElemsY+1
      DO ii=1,nElemsX
        WM(1:nVar,ii,jj) = Fn(1:nVar,ii,jj)
        WP(1:nVar,ii,jj) = Fp(1:nVar,ii,jj)
      END DO
    END DO
  CASE(2)
    DO jj=0,nElemsY+1
      DO ii=1,nElemsX
        CALL MUSCL(&
              Fn(1:nVar,ii,-nGhosts+jj:jj+nGhosts),&
              Fp(1:nVar,ii,-nGhosts+jj:jj+nGhosts),&
              WM(1:nVar,ii,jj),&
              WP(1:nVar,ii,jj),&
              Mesh_DX(2))
      END DO
    END DO
  CASE(3,4)
    DO jj=0,nElemsY+1
      DO ii=1,nElemsX
        IF (.NOT. Ind(2,ii,jj)) THEN
          CALL WENO(&
                V(1:nVar,ii,-nGhosts+jj:jj+nGhosts),&
                Fn(1:nVar,ii,-nGhosts+jj:jj+nGhosts),&
                Fp(1:nVar,ii,-nGhosts+jj:jj+nGhosts),&
                WM(1:nVar,ii,jj),&
                WP(1:nVar,ii,jj),&
                Mesh_DX(2),&
                Reconstruction)
        END IF
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Reconstruction not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!-------------------------------------------------------------------------------!
END SUBROUTINE ReconstructionY
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE ReconstructionFixX()
!-------------------------------------------------------------------------------!
USE MOD_FiniteDifference2D_vars,ONLY: V
USE MOD_FiniteDifference2D_vars,ONLY: WM
USE MOD_FiniteDifference2D_vars,ONLY: WP
USE MOD_FiniteDifference2D_vars,ONLY: Fn
USE MOD_FiniteDifference2D_vars,ONLY: Fp
USE MOD_FiniteDifference2D_vars,ONLY: Ind
USE MOD_FiniteDifference2D_vars,ONLY: MESH_DX
USE MOD_FiniteDifference2D_vars,ONLY: nVar
USE MOD_FiniteDifference2D_vars,ONLY: nElemsX
USE MOD_FiniteDifference2D_vars,ONLY: nElemsY
USE MOD_FiniteDifference2D_vars,ONLY: nGhosts
USE MOD_FiniteDifference2D_vars,ONLY: Reconstruction
USE MOD_FiniteDifference2D_vars,ONLY: ReconstructionFix
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! FORMAL ARGUMENTS
!-------------------------------------------------------------------------------!
! LOCAL VARIABLES
!-------------------------------------------------------------------------------!
INTEGER            :: ii, jj
CHARACTER(LEN=255) :: ErrorMessage
!-------------------------------------------------------------------------------!

IF ((Reconstruction .EQ. 1) .OR. (Reconstruction .EQ. 2)) THEN
  RETURN
END IF

SELECT CASE (ReconstructionFix)
  CASE(1)
    DO jj=1,nElemsY
      DO ii=0,nElemsX+1
        IF (Ind(1,ii,jj)) THEN
          WM(1:nVar,ii,jj) = Fn(1:nVar,ii,jj)
          WP(1:nVar,ii,jj) = Fp(1:nVar,ii,jj)
        END IF
      END DO
    END DO
  CASE(2)
    DO jj=1,nElemsY
      DO ii=0,nElemsX+1
        IF (Ind(1,ii,jj)) THEN
          CALL MUSCL(&
                Fn(1:nVar,-nGhosts+ii:ii+nGhosts,jj),&
                Fp(1:nVar,-nGhosts+ii:ii+nGhosts,jj),&
                WM(1:nVar,ii,jj),&
                WP(1:nVar,ii,jj),&
                Mesh_DX(1))
        END IF
      END DO
    END DO
  CASE(3)
    DO jj=1,nElemsY
      DO ii=0,nElemsX+1
        IF (Ind(1,ii,jj)) THEN
          CALL WENO(&
                V(1:nVar,-nGhosts+ii:ii+nGhosts,jj),&
                Fn(1:nVar,-nGhosts+ii:ii+nGhosts,jj),&
                Fp(1:nVar,-nGhosts+ii:ii+nGhosts,jj),&
                WM(1:nVar,ii,jj),&
                WP(1:nVar,ii,jj),&
                Mesh_DX(1),&
                ReconstructionFix)
        END IF
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "ReconstructionFix not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!-------------------------------------------------------------------------------!
END SUBROUTINE ReconstructionFixX
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE ReconstructionFixY()
!-------------------------------------------------------------------------------!
USE MOD_FiniteDifference2D_vars,ONLY: V
USE MOD_FiniteDifference2D_vars,ONLY: WM
USE MOD_FiniteDifference2D_vars,ONLY: WP
USE MOD_FiniteDifference2D_vars,ONLY: Fn
USE MOD_FiniteDifference2D_vars,ONLY: Fp
USE MOD_FiniteDifference2D_vars,ONLY: Ind
USE MOD_FiniteDifference2D_vars,ONLY: MESH_DX
USE MOD_FiniteDifference2D_vars,ONLY: nVar
USE MOD_FiniteDifference2D_vars,ONLY: nElemsX
USE MOD_FiniteDifference2D_vars,ONLY: nElemsY
USE MOD_FiniteDifference2D_vars,ONLY: nGhosts
USE MOD_FiniteDifference2D_vars,ONLY: Reconstruction
USE MOD_FiniteDifference2D_vars,ONLY: ReconstructionFix
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! FORMAL ARGUMENTS
!-------------------------------------------------------------------------------!
! LOCAL VARIABLES
!-------------------------------------------------------------------------------!
INTEGER            :: ii, jj
CHARACTER(LEN=255) :: ErrorMessage
!-------------------------------------------------------------------------------!

IF ((Reconstruction .EQ. 1) .OR. (Reconstruction .EQ. 2)) THEN
  RETURN
END IF

SELECT CASE (ReconstructionFix)
  CASE(1)
    DO jj=0,nElemsY+1
      DO ii=1,nElemsX
        IF (Ind(2,ii,jj)) THEN
          WM(1:nVar,ii,jj) = Fn(1:nVar,ii,jj)
          WP(1:nVar,ii,jj) = Fp(1:nVar,ii,jj)
        END IF
      END DO
    END DO
  CASE(2)
    DO jj=0,nElemsY+1
      DO ii=1,nElemsX
        IF (Ind(2,ii,jj)) THEN
          CALL MUSCL(&
                Fn(1:nVar,ii,-nGhosts+jj:jj+nGhosts),&
                Fp(1:nVar,ii,-nGhosts+jj:jj+nGhosts),&
                WM(1:nVar,ii,jj),&
                WP(1:nVar,ii,jj),&
                Mesh_DX(2))
        END IF
      END DO
    END DO
  CASE(3)
    DO jj=0,nElemsY+1
      DO ii=1,nElemsX
        IF (Ind(2,ii,jj)) THEN
          CALL WENO(&
                V(1:nVar,ii,-nGhosts+jj:jj+nGhosts),&
                Fn(1:nVar,ii,-nGhosts+jj:jj+nGhosts),&
                Fp(1:nVar,ii,-nGhosts+jj:jj+nGhosts),&
                WM(1:nVar,ii,jj),&
                WP(1:nVar,ii,jj),&
                Mesh_DX(2),&
                ReconstructionFix)
        END IF
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "ReconstructionFix not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!-------------------------------------------------------------------------------!
END SUBROUTINE ReconstructionFixY
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE MUSCL(Fn,Fp,WM,WP,dx)
!-------------------------------------------------------------------------------!
USE MOD_FiniteDifference2D_vars,ONLY: nVar
USE MOD_FiniteDifference2D_vars,ONLY: nGhosts
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! FORMAL ARGUMENTS
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Fn(1:nVar,-nGhosts:nGhosts)
REAL,INTENT(IN)  :: Fp(1:nVar,-nGhosts:nGhosts)
REAL,INTENT(IN)  :: dx
REAL,INTENT(OUT) :: WM(1:nVar)
REAL,INTENT(OUT) :: WP(1:nVar)
!-------------------------------------------------------------------------------!
! LOCAL VARIABLES
!-------------------------------------------------------------------------------!
REAL             :: Q(1:nVar,-nGhosts:nGhosts)
REAL             :: sm, sp, slope
INTEGER          :: iVar
!-------------------------------------------------------------------------------!

Q = Fn
DO iVar = 1, nVar
  sp    = (Q(iVar,+1) - Q(iVar,+0))/dx
  sm    = (Q(iVar,+0) - Q(iVar,-1))/dx
  slope = MINMOD(sm,sp)

  WM(iVar) = Q(iVar,0) - 0.5*slope*dx
END DO

Q = Fp
DO iVar = 1, nVar
  sp    = (Q(iVar,+1) - Q(iVar,+0))/dx
  sm    = (Q(iVar,+0) - Q(iVar,-1))/dx
  slope = MINMOD(sm,sp)

  WP(iVar) = Q(iVar,0) + 0.5*slope*dx
END DO

!-------------------------------------------------------------------------------!
END SUBROUTINE MUSCL
!===============================================================================!
!
!
!
!===============================================================================!
FUNCTION MINMOD(x,y)
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! FORMAL ARGUMENTS
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: x, y
REAL            :: MINMOD
!-------------------------------------------------------------------------------!
! LOCAL VARIABLES
!-------------------------------------------------------------------------------!

MINMOD = 0.5*(SIGN(1.0,x) + SIGN(1.0,y))*MIN(ABS(x),ABS(y))

!-------------------------------------------------------------------------------!
END FUNCTION MINMOD
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE WENO(V,Fn,Fp,WM,WP,dx,WhichReconstruction)
!-------------------------------------------------------------------------------!
USE MOD_FiniteDifference2D_vars,ONLY: nVar
USE MOD_FiniteDifference2D_vars,ONLY: nGhosts
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! FORMAL ARGUMENTS
!-------------------------------------------------------------------------------!
REAL,   INTENT(IN)  :: V(1:nVar,-nGhosts:nGhosts)
REAL,   INTENT(IN)  :: Fn(1:nVar,-nGhosts:nGhosts)
REAL,   INTENT(IN)  :: Fp(1:nVar,-nGhosts:nGhosts)
REAL,   INTENT(IN)  :: dx
INTEGER,INTENT(IN)  :: WhichReconstruction
REAL,   INTENT(OUT) :: WM(1:nVar)
REAL,   INTENT(OUT) :: WP(1:nVar)
!-------------------------------------------------------------------------------!
! LOCAL VARIABLES
!-------------------------------------------------------------------------------!
INTEGER             :: iVar, ii, jj, iGP
CHARACTER(LEN=255)  :: ErrorMessage
!-------------------------------------------------------------------------------!

SELECT CASE(WhichReconstruction)
  CASE(3)
    DO iVar = 1,nVar
      CALL WENO3(Fn(iVar,:),Fp(iVar,:),WM(iVar),WP(iVar))
    END DO
  CASE(4)
    DO iVar = 1,nVar
      CALL WENO5(Fn(iVar,:),Fp(iVar,:),WM(iVar),WP(iVar))
    END DO
  CASE DEFAULT
    ErrorMessage = "Reconstruction not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!-------------------------------------------------------------------------------!
END SUBROUTINE WENO
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE WENO3(Fn,Fp,WM,WP)
!-------------------------------------------------------------------------------!
USE MOD_FiniteDifference2D_vars,ONLY: nGhosts
USE MOD_FiniteDifference2D_vars,ONLY: WENOEPS, WENOEXP
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! FORMAL ARGUMENTS
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Fn(-nGhosts:nGhosts)
REAL,INTENT(IN)  :: Fp(-nGhosts:nGhosts)
REAL,INTENT(OUT) :: WM
REAL,INTENT(OUT) :: WP
!-------------------------------------------------------------------------------!
! LOCAL VARIABLES
!-------------------------------------------------------------------------------!
REAL             :: Q(-nGhosts:nGhosts)
REAL             :: alpha1, alpha2
REAL             :: beta1, beta2
REAL             :: gamma1, gamma2
REAL             :: omega1, omega2
REAL             :: W1, W2
!-------------------------------------------------------------------------------!

!------------------------------!
! WM: x_{i-1/2}                !
!------------------------------!
Q = Fn

! Smoothness Indicators
beta1 = (Q(-1) - Q(+0))**2.0
beta2 = (Q(+0) - Q(+1))**2.0

! Linear Weights
gamma1 = 2.0/3.0
gamma2 = 1.0/3.0

alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP

! Nonlinear Weights
omega1 = alpha1/(alpha1 + alpha2)
omega2 = alpha2/(alpha1 + alpha2)

W1 = 0.5*(    Q(-1) + Q(+0))
W2 = 0.5*(3.0*Q(+0) - Q(+1))
WM = omega1*W1 + omega2*W2


!------------------------------!
! WP: x_{i+1/2}                !
!------------------------------!
Q = Fp

! Smoothness Indicators
beta1 = (Q(-1) - Q(+0))**2.0
beta2 = (Q(+0) - Q(+1))**2.0

! Linear Weights
gamma1 = 1.0/3.0
gamma2 = 2.0/3.0

alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP

! Nonlinear Weights
omega1 = alpha1/(alpha1 + alpha2)
omega2 = alpha2/(alpha1 + alpha2)

! Reconstructed Polynomial
W1 = 0.5*(-Q(-1) + 3.0*Q(+0))
W2 = 0.5*( Q(+0) +     Q(+1))
WP = omega1*W1 + omega2*W2

!-------------------------------------------------------------------------------!
END SUBROUTINE WENO3
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE WENO5(Fn,Fp,WM,WP)
!-------------------------------------------------------------------------------!
USE MOD_FiniteDifference2D_vars,ONLY: nGhosts
USE MOD_FiniteDifference2D_vars,ONLY: WENOEPS, WENOEXP
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! FORMAL ARGUMENTS
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Fn(-nGhosts:nGhosts)
REAL,INTENT(IN)  :: Fp(-nGhosts:nGhosts)
REAL,INTENT(OUT) :: WM
REAL,INTENT(OUT) :: WP
!-------------------------------------------------------------------------------!
! LOCAL VARIABLES
!-------------------------------------------------------------------------------!
REAL             :: Q(-nGhosts:nGhosts)
REAL             :: alpha1, alpha2, alpha3
REAL             :: beta1,  beta2,  beta3
REAL             :: gamma1, gamma2, gamma3
REAL             :: omega1, omega2, omega3
REAL             :: W1, W2, W3
!-------------------------------------------------------------------------------!


!------------------------------!
! WM: x_{i-1/2}^{R}            !
!------------------------------!
Q = Fn

! Smoothness Indicators
beta1 = (1.0/3.0)*( 4.0*Q(-2)*Q(-2) - 19.0*Q(-2)*Q(-1) + 25.0*Q(-1)*Q(-1) &
      + 11.0*Q(-2)*Q(+0) - 31.0*Q(-1)*Q(+0) + 10.0*Q(+0)*Q(+0))
beta2 = (1.0/3.0)*( 4.0*Q(-1)*Q(-1) - 13.0*Q(-1)*Q(+0) + 13.0*Q(+0)*Q(+0) &
      +  5.0*Q(-1)*Q(+1) - 13.0*Q(+0)*Q(+1) +  4.0*Q(+1)*Q(+1))
beta3 = (1.0/3.0)*(10.0*Q(+0)*Q(+0) - 31.0*Q(+0)*Q(+1) + 25.0*Q(+1)*Q(+1) &
      + 11.0*Q(+0)*Q(+2) - 19.0*Q(+1)*Q(+2) +  4.0*Q(+2)*Q(+2))

! Linear Weights
gamma1 = 3.0/10.0
gamma2 = 3.0/5.0
gamma3 = 1.0/10.0

alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP
alpha3 = gamma3/(WENOEPS + beta3)**WENOEXP

! Nonlinear Weights
omega1 = alpha1/(alpha1 + alpha2 + alpha3)
omega2 = alpha2/(alpha1 + alpha2 + alpha3)
omega3 = alpha3/(alpha1 + alpha2 + alpha3)

W1 = (1.0/6.0)*(    -Q(-2) + 5.0*Q(-1) + 2.0*Q(+0))
W2 = (1.0/6.0)*( 2.0*Q(-1) + 5.0*Q(+0) -     Q(+1))
W3 = (1.0/6.0)*(11.0*Q(+0) - 7.0*Q(+1) + 2.0*Q(+2))
WM = omega1*W1 + omega2*W2 + omega3*W3


!------------------------------!
! WP: x_{i+1/2}^{L}            !
!------------------------------!
Q = Fp

! Smoothness Indicators
beta1 = (1.0/3.0)*( 4.0*Q(-2)*Q(-2) - 19.0*Q(-2)*Q(-1) + 25.0*Q(-1)*Q(-1) &
      + 11.0*Q(-2)*Q(+0) - 31.0*Q(-1)*Q(+0) + 10.0*Q(+0)*Q(+0))
beta2 = (1.0/3.0)*( 4.0*Q(-1)*Q(-1) - 13.0*Q(-1)*Q(+0) + 13.0*Q(+0)*Q(+0) &
      +  5.0*Q(-1)*Q(+1) - 13.0*Q(+0)*Q(+1) +  4.0*Q(+1)*Q(+1))
beta3 = (1.0/3.0)*(10.0*Q(+0)*Q(+0) - 31.0*Q(+0)*Q(+1) + 25.0*Q(+1)*Q(+1) &
      + 11.0*Q(+0)*Q(+2) - 19.0*Q(+1)*Q(+2) +  4.0*Q(+2)*Q(+2))

! Linear Weights
gamma1 = 1.0/10.0
gamma2 = 3.0/5.0
gamma3 = 3.0/10.0

alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP
alpha3 = gamma3/(WENOEPS + beta3)**WENOEXP

! Nonlinear Weights
omega1 = alpha1/(alpha1 + alpha2 + alpha3)
omega2 = alpha2/(alpha1 + alpha2 + alpha3)
omega3 = alpha3/(alpha1 + alpha2 + alpha3)

W1 = (1.0/6.0)*(2.0*Q(-2) - 7.0*Q(-1) + 11.0*Q(+0))
W2 = (1.0/6.0)*(   -Q(-1) + 5.0*Q(+0) +  2.0*Q(+1))
W3 = (1.0/6.0)*(2.0*Q(+0) + 5.0*Q(+1) -      Q(+2))
WP = omega1*W1 + omega2*W2 + omega3*W3

!-------------------------------------------------------------------------------!
END SUBROUTINE WENO5
!===============================================================================!
!
!
!
!===============================================================================!
END MODULE MOD_Reconstruction
!-------------------------------------------------------------------------------!
