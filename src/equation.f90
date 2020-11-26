!===============================================================================!
MODULE MOD_Equation
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
PRIVATE
!-------------------------------------------------------------------------------!
INTERFACE ExactFunction
  MODULE PROCEDURE ExactFunction
END INTERFACE

INTERFACE SourceTerms
  MODULE PROCEDURE SourceTerms
END INTERFACE

INTERFACE BoundaryConditions
  MODULE PROCEDURE BoundaryConditions
END INTERFACE

INTERFACE TimeStep
  MODULE PROCEDURE TimeStep
END INTERFACE

INTERFACE ConsToPrim
  MODULE PROCEDURE ConsToPrim
END INTERFACE

INTERFACE PrimToCons
  MODULE PROCEDURE PrimToCons
END INTERFACE

INTERFACE EvaluateFluxX
  MODULE PROCEDURE EvaluateFluxX
END INTERFACE

INTERFACE EvaluateFluxY
  MODULE PROCEDURE EvaluateFluxY
END INTERFACE

INTERFACE MixedGLM
  MODULE PROCEDURE MixedGLM
END INTERFACE
!-------------------------------------------------------------------------------!
PUBLIC :: ExactFunction
PUBLIC :: SourceTerms
PUBLIC :: BoundaryConditions
PUBLIC :: TimeStep
PUBLIC :: ConsToPrim
PUBLIC :: PrimToCons
PUBLIC :: EvaluateFluxX
PUBLIC :: EvaluateFluxY
PUBLIC :: MixedGLM
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
SUBROUTINE ExactFunction(WhichInitialCondition,t,x,Cons)
!-------------------------------------------------------------------------------!
USE MOD_FiniteDifference2D_vars,ONLY: nVar
USE MOD_FiniteDifference2D_vars,ONLY: nDims
USE MOD_FiniteDifference2D_vars,ONLY: PI
USE MOD_FiniteDifference2D_vars,ONLY: MESH_X0
USE MOD_FiniteDifference2D_vars,ONLY: MESH_X1
USE MOD_FiniteDifference2D_vars,ONLY: MESH_SX
USE MOD_FiniteDifference2D_vars,ONLY: Kappa
USE MOD_FiniteDifference2D_vars,ONLY: KappaM1
USE MOD_FiniteDifference2D_vars,ONLY: KappaP1
USE MOD_FiniteDifference2D_vars,ONLY: sKappaM1
USE MOD_FiniteDifference2D_vars,ONLY: PrimRefState1
USE MOD_FiniteDifference2D_vars,ONLY: PrimRefState2
USE MOD_FiniteDifference2D_vars,ONLY: PrimRefState3
USE MOD_FiniteDifference2D_vars,ONLY: PrimRefState4
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! FORMAL ARGUMENTS
!-------------------------------------------------------------------------------!
INTEGER,INTENT(IN) :: WhichInitialCondition
REAL,INTENT(IN)    :: t
REAL,INTENT(IN)    :: x(1:nDims)
REAL,INTENT(OUT)   :: Cons(1:nVar)
!-------------------------------------------------------------------------------!
! LOCAL VARIABLES
!-------------------------------------------------------------------------------!
REAL                :: Prim(1:nVar)
REAL                :: Cons_In(1:nVar),Cons_Out(1:nVar)
REAL                :: r, r0, rho0, p0, B0, vs, rho1, rho2, am, om 
REAL                :: xc(2), xm, ym
REAL                :: rho_in,rho_out,v_in,v_out,p_in,p_out,B_in,B_out
REAL                :: delta_rho, delta_vx, delta_vy, delta_p, delta_bx, delta_by
CHARACTER(LEN=255)  :: ErrorMessage
!-------------------------------------------------------------------------------!

Cons = 0.0
Prim = 0.0
SELECT CASE(WhichInitialCondition)
  !----------------------------------------------------------------------!
  ! [200] Constant State                                                 !
  !----------------------------------------------------------------------!
  CASE(200)
    Prim(1:nVar) = PrimRefState1(1:nVar)

    CALL PrimToCons(Prim,Cons)
  !----------------------------------------------------------------------!
  ! [211] Magnetic Field Loop Advection                                  !
  !----------------------------------------------------------------------!
  CASE(211)
    xm    = MESH_X0(1)+0.5*MESH_SX(1)
    ym    = MESH_X0(2)+0.5*MESH_SX(2)
    xc(1) = x(1)-xm
    xc(2) = x(2)-ym
    r     = SQRT(xc(1)**2 + xc(2)**2)

    rho0  = 1.0
    B0    = 1.0E-03
    p0    = 1.0

    Prim(1) = rho0
    Prim(2) = 2.00
    Prim(3) = 1.00
    Prim(4) = 1.00
    Prim(5) = p0
    Prim(6) = 0.0
    Prim(7) = 0.0
    Prim(8) = 0.0

    IF (r .LT. 0.3) THEN
      Prim(6) = -B0*x(2)/r
      Prim(7) =  B0*x(1)/r
      Prim(8) = 0.0
    END IF

    CALL PrimToCons(Prim,Cons)
  !----------------------------------------------------------------------!
  ! [213] Orszag-Tang Vortex                                             !
  !----------------------------------------------------------------------!
  CASE(213)
    Prim(1)= Kappa**2
    Prim(2)= -SIN(2.0*PI*x(2))
    Prim(3)= +SIN(2.0*PI*x(1))
    Prim(4)= 0.0
    Prim(5)= Kappa
    Prim(6)= -SIN(2.0*PI*x(2))
    Prim(7)= +SIN(4.0*PI*x(1))
    Prim(8)= 0.0

    CALL PrimToCons(Prim,Cons)
  !----------------------------------------------------------------------!
  ! [215] Rotor Problem                                                  !
  !----------------------------------------------------------------------!
  CASE(215)
    Cons_In  = 0.0
    Cons_Out = 0.0

    xm      = MESH_X0(1)+0.5*MESH_SX(1)
    ym      = MESH_X0(2)+0.5*MESH_SX(2)
    xc(1)   = x(1)-xm
    xc(2)   = x(2)-ym
    r       = SQRT(xc(1)**2 + xc(2)**2)

    r0      = 1.0E-1
    om      = 1.0E+0

    rho_in  = 1.0E+1
    rho_out = 1.0E+0
    v_in    = 0.0E+0
    v_out   = 0.0E+0
    p_in    = 1.0E+0
    p_out   = 1.0E+0
    B_in    = 0.5*SQRT(2.0)
    B_out   = 0.5*SQRT(2.0)

    Prim(1) = rho_in
    Prim(2) = -om*xc(2)/r0
    Prim(3) = +om*xc(1)/r0
    Prim(4) = v_in
    Prim(5) = p_in
    Prim(6) = B_in
    Prim(7) = 0.0
    Prim(8) = 0.0

    CALL PrimToCons(Prim,Cons_In)

    Prim(1) = rho_out
    Prim(2) = v_out
    Prim(3) = v_out
    Prim(4) = v_out
    Prim(5) = p_out
    Prim(6) = B_out
    Prim(7) = 0.0
    Prim(8) = 0.0
    CALL PrimToCons(Prim,Cons_Out)
    Cons    = -0.5*(Cons_In-Cons_Out)*TANH(80.0*(r-r0)) + Cons_Out + 0.5*(Cons_In-Cons_Out)
  !----------------------------------------------------------------------!
  ! [217] Kelvin-Helmholtz Instability                                   !
  !----------------------------------------------------------------------!
  CASE(217)
    vs = 0.50
    am = 0.01

    rho1 = 1.0
    rho2 = 2.0
    p0   = 2.5
    B0   = 0.2

    IF (ABS(x(2)) .GE. 0.25) THEN
      delta_rho = rho1
      delta_vx  = +(vs+am*SIN(2.0*PI*x(1)))
    ELSE
      delta_rho = rho2
      delta_vx  = -(vs+am*SIN(2.0*PI*x(1)))
    END IF
    IF (x(2) .GE. 0.00) THEN
      delta_vy  = +am*SIN(2.0*PI*x(1))
    ELSE
      delta_vy  = +am*SIN(2.0*PI*x(1))
    END IF

    Prim(1) = delta_rho
    Prim(2) = delta_vx
    Prim(3) = delta_vy
    Prim(4) = 0.0
    Prim(5) = p0
    Prim(6) = B0
    Prim(7) = 0.0
    Prim(8) = 0.0
    CALL PrimToCons(Prim,Cons)
  CASE DEFAULT
    ErrorMessage = "Exact function not specified"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!-------------------------------------------------------------------------------!
END SUBROUTINE ExactFunction
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE SourceTerms(t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteDifference2D_vars,ONLY: S
USE MOD_FiniteDifference2D_vars,ONLY: MeshNodes
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! FORMAL ARGUMENTS
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: t
!-------------------------------------------------------------------------------!
! LOCAL VARIABLES
!-------------------------------------------------------------------------------!
INTEGER          :: ii, jj
!-------------------------------------------------------------------------------!

S = 0.0

!-------------------------------------------------------------------------------!
END SUBROUTINE SourceTerms
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE BoundaryConditions(t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteDifference2D_vars,ONLY: nVar
USE MOD_FiniteDifference2D_vars,ONLY: nElemsX
USE MOD_FiniteDifference2D_vars,ONLY: nElemsY
USE MOD_FiniteDifference2D_vars,ONLY: nGhosts
USE MOD_FiniteDifference2D_vars,ONLY: MESH_X0
USE MOD_FiniteDifference2D_vars,ONLY: MESH_X1
USE MOD_FiniteDifference2D_vars,ONLY: MeshNodes
USE MOD_FiniteDifference2D_vars,ONLY: PrimRefState1
USE MOD_FiniteDifference2D_vars,ONLY: PrimRefState2
USE MOD_FiniteDifference2D_vars,ONLY: PrimRefState3
USE MOD_FiniteDifference2D_vars,ONLY: PrimRefState4
USE MOD_FiniteDifference2D_vars,ONLY: BoundaryConditionsType
USE MOD_FiniteDifference2D_vars,ONLY: U
USE MOD_FiniteDifference2D_vars,ONLY: V
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! FORMAL ARGUMENTS
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)    :: t
!-------------------------------------------------------------------------------!
! LOCAL VARIABLES
!-------------------------------------------------------------------------------!
INTEGER            :: ii, jj
INTEGER            :: idx_vx, idx_vy
REAL               :: x0, xc, xt
REAL               :: Prim_in(1:nVar), Prim_out(1:nVar)
REAL               :: Cons_in(1:nVar), Cons_out(1:nVar)
REAL               :: ConsRefState1(1:nVar), ConsRefState2(1:nVar)
CHARACTER(LEN=255) :: ErrorMessage
!-------------------------------------------------------------------------------!

idx_vx = 2
idx_vy = 3

!------------------------------!
! Left Boundary Conditions     !
!------------------------------!
SELECT CASE(BoundaryConditionsType(4))
  CASE(1) ! Periodic
    DO jj=1,nElemsY
      DO ii=1,nGhosts
        U(1:nVar,-nGhosts+ii,jj) = U(1:nVar,nElemsX-nGhosts+ii,jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO jj=1,nElemsY
      DO ii=1,nGhosts
        U(1:nVar,-nGhosts+ii,jj) = U(1:nVar,nGhosts-ii+1,jj)
      END DO
    END DO
  CASE(3) ! Inflow
    Prim_in(1:nVar) = PrimRefState4(1:nVar)
    CALL PrimToCons(Prim_in(1:nVar),Cons_in(1:nVar))
    DO jj=1,nElemsY
      DO ii=1,nGhosts
        U(1:nVar,-nGhosts+ii,jj) = Cons_in(1:nVar)
      END DO
    END DO
  CASE(4) ! Outflow
    Prim_out(1:nVar) = PrimRefState4(1:nVar)
    CALL PrimToCons(Prim_out(1:nVar),Cons_out(1:nVar))
    DO jj=1,nElemsY
      DO ii=1,nGhosts
        U(1:nVar,-nGhosts+ii,jj) = Cons_out(1:nVar)
      END DO
    END DO
  CASE(5) ! Reflecting
    DO jj=1,nElemsY
      DO ii=1,nGhosts
        U(1:nVar,-nGhosts+ii,jj) = U(1:nVar,nGhosts-ii+1,jj)
        U(idx_vx,-nGhosts+ii,jj) =-U(idx_vx,nGhosts-ii+1,jj)
      END DO
    END DO
  CASE(11) ! Double Mach Reflection
    ErrorMessage = "Boundary condition defined only for faces 1 and 3"
    WRITE(*,*) ErrorMessage
    STOP
  CASE DEFAULT
    ErrorMessage = "Boundary condition not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!------------------------------!
! Right Boundary Conditions    !
!------------------------------!
SELECT CASE(BoundaryConditionsType(2))
  CASE(1) ! Periodic
    DO jj=1,nElemsY
      DO ii=1,nGhosts
        U(1:nVar,nElemsX+ii,jj) = U(1:nVar,ii,jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO jj=1,nElemsY    
      DO ii=1,nGhosts
        U(1:nVar,nElemsX+ii,jj) = U(1:nVar,nElemsX-ii+1,jj)
      END DO
    END DO
  CASE(3) ! Inflow
    Prim_in(1:nVar) = PrimRefState2(1:nVar)
    CALL PrimToCons(Prim_in(1:nVar),Cons_in(1:nVar))
    DO jj=1,nElemsY
      DO ii=1,nGhosts
        U(1:nVar,nElemsX+ii,jj) = Cons_in(1:nVar)
      END DO
    END DO
  CASE(4) ! Outflow
    Prim_out(1:nVar) = PrimRefState2(1:nVar)
    CALL PrimToCons(Prim_out(1:nVar),Cons_out(1:nVar))
    DO jj=1,nElemsY
      DO ii=1,nGhosts
        U(1:nVar,nElemsX+ii,jj) = Cons_out(1:nVar)
      END DO
    END DO
  CASE(5) ! Reflecting
    DO jj=1,nElemsY
      DO ii=1,nGhosts
        U(1:nVar,nElemsX+ii,jj) = U(1:nVar,nElemsX-ii+1,jj)
        U(idx_vx,nElemsX+ii,jj) =-U(idx_vx,nElemsX-ii+1,jj)
      END DO
    END DO
  CASE(11) ! Double Mach Reflection
    ErrorMessage = "Boundary condition defined only for faces 1 and 3"
    WRITE(*,*) ErrorMessage
    STOP
  CASE DEFAULT
    ErrorMessage = "Boundary condition not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!------------------------------!
! Top Boundary Conditions      !
!------------------------------!
SELECT CASE(BoundaryConditionsType(3))
  CASE(1) ! Periodic
    DO ii=1,nElemsX
      DO jj=1,nGhosts
        U(1:nVar,ii,nElemsY+jj) = U(1:nVar,ii,jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO ii=1,nElemsX
      DO jj=1,nGhosts
        U(1:nVar,ii,nElemsY+jj) = U(1:nVar,ii,nElemsY-jj+1)
      END DO
    END DO
  CASE(3) ! Inflow
    Prim_in(1:nVar) = PrimRefState3(1:nVar)
    CALL PrimToCons(Prim_in(1:nVar),Cons_in(1:nVar))
    DO ii=1,nElemsX
      DO jj=1,nGhosts
        U(1:nVar,ii,nElemsY+jj) = Cons_in(1:nVar)
      END DO
    END DO
  CASE(4) ! Outflow
    Prim_out(1:nVar) = PrimRefState3(1:nVar)
    CALL PrimToCons(Prim_out(1:nVar),Cons_out(1:nVar))
    DO ii=1,nElemsX
      DO jj=1,nGhosts
        U(1:nVar,ii,nElemsY+jj) = Cons_out(1:nVar)
      END DO
    END DO
  CASE(5) ! Reflecting
    DO ii=1,nElemsX
      DO jj=1,nGhosts
        U(1:nVar,ii,nElemsY+jj) = U(1:nVar,ii,nElemsY-jj+1)
        U(idx_vy,ii,nElemsY+jj) =-U(idx_vy,ii,nElemsY-jj+1)
      END DO
    END DO
  CASE(11) ! Double Mach Reflection
    Prim_in(1:nVar)  = PrimRefState1(1:nVar)
    Prim_out(1:nVar) = PrimRefState3(1:nVar)
    CALL PrimToCons(Prim_in(1:nVar),Cons_in(1:nVar))
    CALL PrimToCons(Prim_out(1:nVar),Cons_out(1:nVar))
    DO ii=1,nElemsX
      DO jj=1,nGhosts
        x0 = 1.0/6.0
        xc = MeshNodes(1,ii,nElemsY)
        xt = x0 + (1.0+2.0*10.0*t)/SQRT(3.0)
        IF (xc .LT. xt) THEN
          U(1:nVar,ii,nElemsY+jj) = Cons_in(1:nVar)
        ELSE
          U(1:nVar,ii,nElemsY+jj) = Cons_out(1:nVar)
        END IF
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!------------------------------!
! Bottom Boundary Conditions   !
!------------------------------!
SELECT CASE(BoundaryConditionsType(1))
  CASE(1) ! Periodic
    DO ii=1,nElemsX
      DO jj=1,nGhosts
        U(1:nVar,ii,-nGhosts+jj) = U(1:nVar,ii,nElemsY-nGhosts+jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO ii=1,nElemsX
      DO jj=1,nGhosts
        U(1:nVar,ii,-nGhosts+jj) = U(1:nVar,ii,nGhosts-jj+1)
      END DO
    END DO
  CASE(3) ! Inflow
    Prim_in(1:nVar) = PrimRefState1(1:nVar)
    CALL PrimToCons(Prim_in(1:nVar),Cons_in(1:nVar))
    DO ii=1,nElemsX
      DO jj=1,nGhosts
        U(1:nVar,ii,-nGhosts+jj) = Cons_in(1:nVar)
      END DO
    END DO
  CASE(4) ! Outflow
    Prim_out(1:nVar) = PrimRefState1(1:nVar)
    CALL PrimToCons(Prim_out(1:nVar),Cons_out(1:nVar))
    DO ii=1,nElemsX
      DO jj=1,nGhosts
        U(1:nVar,ii,-nGhosts+jj) = Cons_out(1:nVar)
      END DO
    END DO
  CASE(5) ! Reflecting
    DO ii=1,nElemsX
      DO jj=1,nGhosts
        U(1:nVar,ii,-nGhosts+jj) = U(1:nVar,ii,nGhosts-jj+1)
        U(idx_vy,ii,-nGhosts+jj) =-U(idx_vy,ii,nGhosts-jj+1)
      END DO
    END DO
  CASE(11) ! Double Mach Reflection
    Prim_in(1:nVar) = PrimRefState1(1:nVar)
    CALL PrimToCons(Prim_in(1:nVar),Cons_in(1:nVar))
    DO ii=1,nElemsX
      DO jj=1,nGhosts
        x0 = 1.0/6.0
        xc = MeshNodes(1,ii,1)
        IF (xc .LT. x0) THEN
          U(1:nVar,ii,-nGhosts+jj) = Cons_in(1:nVar)
        ELSE
          U(1:nVar,ii,-nGhosts+jj) = U(1:nVar,ii,nGhosts-jj+1)
          U(idx_vy,ii,-nGhosts+jj) =-U(idx_vy,ii,nGhosts-jj+1)
        END IF
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

DO jj=-nGhosts,nElemsY+nGhosts+1
  DO ii=-nGhosts,nElemsX+nGhosts+1
    CALL ConsToPrim(U(1:nVar,ii,jj),V(1:nVar,ii,jj))
  END DO
END DO

!-------------------------------------------------------------------------------!
END SUBROUTINE BoundaryConditions
!===============================================================================!
!
!
!
!===============================================================================!
FUNCTION TimeStep()
!-------------------------------------------------------------------------------!
USE MOD_FiniteDifference2D_vars,ONLY: U
USE MOD_FiniteDifference2D_vars,ONLY: CFL
USE MOD_FiniteDifference2D_vars,ONLY: MESH_DX
USE MOD_FiniteDifference2D_vars,ONLY: nVar
USE MOD_FiniteDifference2D_vars,ONLY: nElemsX
USE MOD_FiniteDifference2D_vars,ONLY: nElemsY
USE MOD_FiniteDifference2D_vars,ONLY: LambdaMaxX
USE MOD_FiniteDifference2D_vars,ONLY: LambdaMaxY
USE MOD_FiniteDifference2D_vars,ONLY: MIN_TIMESTEP
USE MOD_FiniteDifference2D_vars,ONLY: GLM_ch
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! FORMAL ARGUMENTS
!-------------------------------------------------------------------------------!
REAL    :: TimeStep
!-------------------------------------------------------------------------------!
! LOCAL VARIABLES
!-------------------------------------------------------------------------------!
REAL    :: FastestWaveX, FastestWaveY
REAL    :: Prim(1:nVar)
INTEGER :: ii, jj
!-------------------------------------------------------------------------------!

GLM_ch = 0.0
LambdaMaxX = 0.0
LambdaMaxY = 0.0
TimeStep = HUGE(1.0)

DO jj=1,nElemsY
  DO ii=1,nElemsX
    CALL ConsToPrim(U(1:nVar,ii,jj),Prim(1:nVar))
    CALL WaveSpeeds2D(Prim(1:nVar),FastestWaveX,FastestWaveY)
    LambdaMaxX = MAX(LambdaMaxX,ABS(FastestWaveX))
    LambdaMaxY = MAX(LambdaMaxY,ABS(FastestWaveY))
    TimeStep  = MIN(TimeStep,MESH_DX(1)/LambdaMaxX,MESH_DX(2)/LambdaMaxY)
    GLM_ch = MAX(GLM_ch,LambdaMaxX,LambdaMaxY)
  END DO
END DO

TimeStep = CFL*TimeStep

IF (TimeStep .LT. MIN_TIMESTEP) THEN
  TimeStep = MIN_TIMESTEP
END IF

!-------------------------------------------------------------------------------!
END FUNCTION TimeStep
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE WaveSpeeds2D(Prim,fastestx,fastesty)
!-------------------------------------------------------------------------------!
USE MOD_FiniteDifference2D_vars,ONLY: nVar
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! FORMAL ARGUMENTS
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Prim(1:nVar)
REAL,INTENT(OUT) :: fastestx
REAL,INTENT(OUT) :: fastesty
!-------------------------------------------------------------------------------!
! LOCAL VARIABLES
!-------------------------------------------------------------------------------!
REAL             :: rho, vx, vy, vz, p, Bx, By, Bz
REAL             :: c, v2, B2, c2
REAL             :: cfx, cfy, cax2, cay2, cfx2, cfy2
!-------------------------------------------------------------------------------!

rho = Prim(1)
vx  = Prim(2)
vy  = Prim(3)
vz  = Prim(4)
p   = Prim(5)
Bx  = Prim(6)
By  = Prim(7)
Bz  = Prim(8)

c  = EOS_SoundSpeed(rho,p)
c2 = c**2
v2 = vx*vx + vy*vy + vz*vz
B2 = Bx*Bx + By*By + Bz*Bz

cax2 = Bx**2/rho
cfx2 = 0.5*((c2+B2/rho) + SQRT((c2+B2/rho)**2.0 - 4.0*c2*cax2))
cfx  = SQRT(cfx2)

cay2 = By**2/rho
cfy2 = 0.5*((c2+B2/rho) + SQRT((c2+B2/rho)**2.0 - 4.0*c2*cay2))
cfy  = SQRT(cfy2)

fastestx = vx + cfx
fastesty = vy + cfy

!-------------------------------------------------------------------------------!
END SUBROUTINE WaveSpeeds2D
!===============================================================================!
!
!
!
!===============================================================================!
FUNCTION EOS_SoundSpeed(rho,p)
!-------------------------------------------------------------------------------!
USE MOD_FiniteDifference2D_vars,ONLY: Kappa
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! FORMAL ARGUMENTS
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: rho, p
REAL            :: EOS_SoundSpeed
!-------------------------------------------------------------------------------!
! LOCAL VARIABLES
!-------------------------------------------------------------------------------!
CHARACTER(LEN=255) :: ErrorMessage
!-------------------------------------------------------------------------------!

EOS_SoundSpeed = ABS(Kappa*p/rho)
IF (EOS_SoundSpeed .LT. 0.0) THEN
  ErrorMessage = "Negative speed of sound"
  WRITE(*,*) ErrorMessage
  STOP
END IF
EOS_SoundSpeed = SQRT(EOS_SoundSpeed)

!-------------------------------------------------------------------------------!
END FUNCTION EOS_SoundSpeed
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE ConsToPrim(Cons,Prim)
!-------------------------------------------------------------------------------!
USE MOD_FiniteDifference2D_vars,ONLY: nVar
USE MOD_FiniteDifference2D_vars,ONLY: Kappa
USE MOD_FiniteDifference2D_vars,ONLY: KappaM1
USE MOD_FiniteDifference2D_vars,ONLY: sKappaM1
USE MOD_FiniteDifference2D_vars,ONLY: MIN_DENSITY, MIN_ENERGY
USE MOD_FiniteDifference2D_vars,ONLY: MIN_MOMENTUM, MIN_MAGNETIC
USE MOD_FiniteDifference2D_vars,ONLY: MIN_PSI
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! FORMAL ARGUMENTS
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Cons(1:nVar)
REAL,INTENT(OUT) :: Prim(1:nVar)
!-------------------------------------------------------------------------------!
! LOCAL VARIABLES
!-------------------------------------------------------------------------------!
REAL             :: rho, Sx, Sy, Sz, E, Bx, By, Bz
REAL             :: S2, B2
REAL             :: Psi
!-------------------------------------------------------------------------------!

rho = Cons(1)
Sx  = Cons(2)
Sy  = Cons(3)
Sz  = Cons(4)
E   = Cons(5)
Bx  = Cons(6)
By  = Cons(7)
Bz  = Cons(8)
Psi = Cons(9)

S2 = Sx**2 + Sy**2 + Sz**2
B2 = Bx**2 + By**2 + Bz**2

IF (rho .LT. MIN_DENSITY) THEN
  rho = MIN_DENSITY
END IF
IF (E .LT. MIN_ENERGY) THEN
  E = MIN_ENERGY
END IF
IF (ABS(Sx) .LT. MIN_MOMENTUM) THEN
  Sx = 0.0
END IF
IF (ABS(Sy) .LT. MIN_MOMENTUM) THEN
  Sy = 0.0
END IF
IF (ABS(Sz) .LT. MIN_MOMENTUM) THEN
  Sz = 0.0
END IF
IF (ABS(Bx) .LT. MIN_MAGNETIC) THEN
  Bx = 0.0
END IF
IF (ABS(By) .LT. MIN_MAGNETIC) THEN
  By = 0.0
END IF
IF (ABS(Bz) .LT. MIN_MAGNETIC) THEN
  Bz = 0.0
END IF
IF (ABS(Psi) .LT. MIN_PSI) THEN
  Psi = 0.0
END IF

Prim(1) = rho
Prim(2) = Sx/rho
Prim(3) = Sy/rho
Prim(4) = Sz/rho
Prim(5) = KappaM1*(E - 0.5*S2/rho - 0.5*B2)
Prim(6) = Bx
Prim(7) = By
Prim(8) = Bz
Prim(9) = Psi


!-------------------------------------------------------------------------------!
END SUBROUTINE ConsToPrim
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE PrimToCons(Prim,Cons)
!-------------------------------------------------------------------------------!
USE MOD_FiniteDifference2D_vars,ONLY: nVar
USE MOD_FiniteDifference2D_vars,ONLY: Kappa
USE MOD_FiniteDifference2D_vars,ONLY: KappaM1
USE MOD_FiniteDifference2D_vars,ONLY: sKappaM1
USE MOD_FiniteDifference2D_vars,ONLY: MIN_DENSITY, MIN_PRESSURE
USE MOD_FiniteDifference2D_vars,ONLY: MIN_SPEED, MIN_MAGNETIC
USE MOD_FiniteDifference2D_vars,ONLY: MIN_PSI
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! FORMAL ARGUMENTS
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Prim(1:nVar)
REAL,INTENT(OUT) :: Cons(1:nVar)
!-------------------------------------------------------------------------------!
! LOCAL VARIABLES
!-------------------------------------------------------------------------------!
REAL             :: rho, vx, vy, vz, p, Bx, By, Bz
REAL             :: v2, B2
REAL             :: Psi
!-------------------------------------------------------------------------------!

rho = Prim(1)
vx  = Prim(2)
vy  = Prim(3)
vz  = Prim(4)
p   = Prim(5)
Bx  = Prim(6)
By  = Prim(7)
Bz  = Prim(8)
Psi = Prim(9)

v2 = vx**2 + vy**2 + vz**2
B2 = Bx**2 + By**2 + Bz**2

IF (rho .LT. MIN_DENSITY) THEN
  rho = MIN_DENSITY
END IF
IF (p .LT. MIN_PRESSURE) THEN
  p = MIN_PRESSURE
END IF
IF (ABS(vx) .LT. MIN_SPEED) THEN
  vx = 0.0
END IF
IF (ABS(vy) .LT. MIN_SPEED) THEN
  vy = 0.0
END IF
IF (ABS(vz) .LT. MIN_SPEED) THEN
  vz = 0.0
END IF
IF (ABS(Bx) .LT. MIN_MAGNETIC) THEN
  Bx = 0.0
END IF
IF (ABS(By) .LT. MIN_MAGNETIC) THEN
  By = 0.0
END IF
IF (ABS(Bz) .LT. MIN_MAGNETIC) THEN
  Bz = 0.0
END IF
IF (ABS(Psi) .LT. MIN_PSI) THEN
  Psi = 0.0
END IF

Cons(1) = rho
Cons(2) = rho*vx
Cons(3) = rho*vy
Cons(4) = rho*vz
Cons(5) = sKappaM1*p + 0.5*rho*v2 + 0.5*B2
Cons(6) = Bx
Cons(7) = By
Cons(8) = Bz
Cons(9) = Psi

!-------------------------------------------------------------------------------!
END SUBROUTINE PrimToCons
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE EvaluateFluxX(Prim,Flux)
!-------------------------------------------------------------------------------!
USE MOD_FiniteDifference2D_vars,ONLY: nVar
USE MOD_FiniteDifference2D_vars,ONLY: Kappa
USE MOD_FiniteDifference2D_vars,ONLY: GLM_ch
USE MOD_FiniteDifference2D_vars,ONLY: sKappaM1
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! FORMAL ARGUMENTS
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Prim(1:nVar)
REAL,INTENT(OUT) :: Flux(1:nVar)
!-------------------------------------------------------------------------------!
! LOCAL VARIABLES
!-------------------------------------------------------------------------------!
REAL             :: rho, vx, vy, vz, p, Bx, By, Bz
REAL             :: v2, B2, pt
REAL             :: Psi
!-------------------------------------------------------------------------------!

rho = Prim(1)
vx  = Prim(2)
vy  = Prim(3)
vz  = Prim(4)
p   = Prim(5)
Bx  = Prim(6)
By  = Prim(7)
Bz  = Prim(8)
Psi = Prim(9)

v2 = vx**2 + vy**2 + vz**2
B2 = Bx**2 + By**2 + Bz**2
pt = (p + 0.5*B2) - Bx*Bx

Flux(1) = rho*vx
Flux(2) = rho*vx*vx + pt
Flux(3) = rho*vy*vx - By*Bx
Flux(4) = rho*vz*vx - Bz*Bx
Flux(5) = (sKappaM1*p + 0.5*rho*v2 + 0.5*B2)*vx + pt*vx - Bx*By*vy - Bx*Bz*vz
Flux(6) = Psi
Flux(7) = By*vx - Bx*vy
Flux(8) = Bz*vx - Bx*vz
Flux(9) = GLM_ch*GLM_ch*Bx

!-------------------------------------------------------------------------------!
END SUBROUTINE EvaluateFluxX
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE EvaluateFluxY(Prim,Flux)
!-------------------------------------------------------------------------------!
USE MOD_FiniteDifference2D_vars,ONLY: nVar
USE MOD_FiniteDifference2D_vars,ONLY: Kappa
USE MOD_FiniteDifference2D_vars,ONLY: GLM_ch
USE MOD_FiniteDifference2D_vars,ONLY: sKappaM1
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! FORMAL ARGUMENTS
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Prim(1:nVar)
REAL,INTENT(OUT) :: Flux(1:nVar)
!-------------------------------------------------------------------------------!
! LOCAL VARIABLES
!-------------------------------------------------------------------------------!
REAL             :: rho, vx, vy, vz, p, Bx, By, Bz
REAL             :: v2, B2, pt
REAL             :: Psi
!-------------------------------------------------------------------------------!

rho = Prim(1)
vx  = Prim(2)
vy  = Prim(3)
vz  = Prim(4)
p   = Prim(5)
Bx  = Prim(6)
By  = Prim(7)
Bz  = Prim(8)
Psi = Prim(9)

v2 = vx**2 + vy**2 + vz**2
B2 = Bx**2 + By**2 + Bz**2
pt = (p + 0.5*B2) - By*By

Flux(1) = rho*vy
Flux(2) = rho*vx*vy - Bx*By
Flux(3) = rho*vy*vy + pt
Flux(4) = rho*vz*vy - Bz*By
Flux(5) = (sKappaM1*p + 0.5*rho*v2 + 0.5*B2)*vy + pt*vy - By*Bx*vx - By*Bz*vz
Flux(6) = Bx*vy - By*vx
Flux(7) = Psi
Flux(8) = Bz*vy - By*vz
Flux(9) = GLM_ch*GLM_ch*By

!-------------------------------------------------------------------------------!
END SUBROUTINE EvaluateFluxY
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE MixedGLM
!-------------------------------------------------------------------------------!
! Divergence Cleaning with mixed GLM
!-------------------------------------------------------------------------------!
USE MOD_FiniteDifference2D_vars,ONLY: U
USE MOD_FiniteDifference2D_vars,ONLY: CFL
USE MOD_FiniteDifference2D_vars,ONLY: nVar
USE MOD_FiniteDifference2D_vars,ONLY: GLM_alpha
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! FORMAL ARGUMENTS
!-------------------------------------------------------------------------------!
! LOCAL VARIABLES
!-------------------------------------------------------------------------------!

U(nVar,:,:) = U(nVar,:,:)*EXP(-CFL*GLM_alpha)

!===============================================================================!
END SUBROUTINE MixedGLM
!-------------------------------------------------------------------------------!
!
!
!
!===============================================================================!
END MODULE MOD_Equation
!-------------------------------------------------------------------------------!
