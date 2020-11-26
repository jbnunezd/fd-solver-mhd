!===============================================================================!
MODULE MOD_FiniteDifference2D
!===============================================================================!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
PRIVATE
!-------------------------------------------------------------------------------!
INTERFACE InitializeFiniteDifference
  MODULE PROCEDURE InitializeFiniteDifference
END INTERFACE

INTERFACE FillInitialConditions
  MODULE PROCEDURE FillInitialConditions
END INTERFACE

INTERFACE FDTimeDerivative
  MODULE PROCEDURE FDTimeDerivative
END INTERFACE

INTERFACE FinalizeFiniteDifference
  MODULE PROCEDURE FinalizeFiniteDifference
END INTERFACE
!-------------------------------------------------------------------------------!
PUBLIC :: InitializeFiniteDifference
PUBLIC :: FillInitialConditions
PUBLIC :: FDTimeDerivative
PUBLIC :: FinalizeFiniteDifference
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
SUBROUTINE InitializeFiniteDifference()
!-------------------------------------------------------------------------------!
USE MOD_FiniteDifference2D_vars,ONLY: nVar
USE MOD_FiniteDifference2D_vars,ONLY: nDims
USE MOD_FiniteDifference2D_vars,ONLY: nGhosts
USE MOD_FiniteDifference2D_vars,ONLY: nElemsX
USE MOD_FiniteDifference2D_vars,ONLY: nElemsY
USE MOD_FiniteDifference2D_vars,ONLY: MeshNodes
USE MOD_FiniteDifference2D_vars,ONLY: U
USE MOD_FiniteDifference2D_vars,ONLY: V
USE MOD_FiniteDifference2D_vars,ONLY: Ut
USE MOD_FiniteDifference2D_vars,ONLY: S
USE MOD_FiniteDifference2D_vars,ONLY: WM
USE MOD_FiniteDifference2D_vars,ONLY: WP
USE MOD_FiniteDifference2D_vars,ONLY: FX
USE MOD_FiniteDifference2D_vars,ONLY: FY
USE MOD_FiniteDifference2D_vars,ONLY: Fn
USE MOD_FiniteDifference2D_vars,ONLY: Fp
USE MOD_FiniteDifference2D_vars,ONLY: Ind
USE MOD_FiniteDifference2D_vars,ONLY: Reconstruction
USE MOD_FiniteDifference2D_vars,ONLY: K0
USE MOD_FiniteDifference2D_vars,ONLY: K1
USE MOD_FiniteDifference2D_vars,ONLY: K2
USE MOD_FiniteDifference2D_vars,ONLY: K3
USE MOD_FiniteDifference2D_vars,ONLY: K4
USE MOD_FiniteDifference2D_vars,ONLY: K5
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
INTEGER :: ii
CHARACTER(LEN=255) :: ErrorMessage
!-------------------------------------------------------------------------------!

SELECT CASE (Reconstruction)
  CASE(1) ! NONE
    nGhosts = 1
  CASE(2) ! MUSCL
    nGhosts = 1
  CASE(3) ! WENO3
    nGhosts = 1
  CASE(4) ! WENO5
    nGhosts = 2
  CASE DEFAULT
    ErrorMessage = "Reconstruction not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

ALLOCATE(MeshNodes(1:nDims,1:nElemsX,1:nElemsY))

ALLOCATE( U(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1))
ALLOCATE( V(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1))
ALLOCATE(Fn(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1))
ALLOCATE(Fp(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1))
ALLOCATE(WM(1:nVar,0:nElemsX+1,0:nElemsY+1))
ALLOCATE(WP(1:nVar,0:nElemsX+1,0:nElemsY+1))
ALLOCATE( S(1:nVar,1:nElemsX,1:nElemsY))
ALLOCATE(Ut(1:nVar,1:nElemsX,1:nElemsY))
ALLOCATE(FX(1:nVar,0:nElemsX,1:nElemsY))
ALLOCATE(FY(1:nVar,1:nElemsX,0:nElemsY))
ALLOCATE(Ind(1:2,0:nElemsX+1,0:nElemsY+1))

ALLOCATE(K0(1:nVar,1:nElemsX,1:nElemsY))
ALLOCATE(K1(1:nVar,1:nElemsX,1:nElemsY))
ALLOCATE(K2(1:nVar,1:nElemsX,1:nElemsY))
ALLOCATE(K3(1:nVar,1:nElemsX,1:nElemsY))
ALLOCATE(K4(1:nVar,1:nElemsX,1:nElemsY))
ALLOCATE(K5(1:nVar,1:nElemsX,1:nElemsY))

U   = 0.0
V   = 0.0
S   = 0.0
Ut  = 0.0
Fn  = 0.0
Fp  = 0.0
WM  = 0.0
WP  = 0.0
FX  = 0.0
FY  = 0.0
Ind = .FALSE.

K0 = 0.0
K1 = 0.0
K2 = 0.0
K3 = 0.0
K4 = 0.0
K5 = 0.0

!-------------------------------------------------------------------------------!
END SUBROUTINE InitializeFiniteDifference
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE FillInitialConditions()
!-------------------------------------------------------------------------------!
USE MOD_Equation,               ONLY: ConsToPrim
USE MOD_Equation,               ONLY: ExactFunction
USE MOD_FiniteDifference2D_vars,ONLY: U
USE MOD_FiniteDifference2D_vars,ONLY: V
USE MOD_FiniteDifference2D_vars,ONLY: nVar
USE MOD_FiniteDifference2D_vars,ONLY: nDims
USE MOD_FiniteDifference2D_vars,ONLY: nElemsX
USE MOD_FiniteDifference2D_vars,ONLY: nElemsY
USE MOD_FiniteDifference2D_vars,ONLY: MeshNodes
USE MOD_FiniteDifference2D_vars,ONLY: InitialCondition
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! FORMAL ARGUMENTS
!-------------------------------------------------------------------------------!
! LOCAL VARIABLES
!-------------------------------------------------------------------------------!
INTEGER :: ii, jj
!-------------------------------------------------------------------------------!

DO jj=1,nElemsY
  DO ii=1,nElemsX
    CALL ExactFunction(&
      InitialCondition,0.0,MeshNodes(1:nDims,ii,jj),U(1:nVar,ii,jj))
    CALL ConsToPrim(U(1:nVar,ii,jj),V(1:nVar,ii,jj))
  END DO
END DO

!-------------------------------------------------------------------------------!
END SUBROUTINE FillInitialConditions
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE FDTimeDerivative(t)
!-------------------------------------------------------------------------------!
USE MOD_Equation,       ONLY: SourceTerms
USE MOD_Equation,       ONLY: BoundaryConditions
USE MOD_Reconstruction, ONLY: ReconstructionX
USE MOD_Reconstruction, ONLY: ReconstructionY
USE MOD_Reconstruction, ONLY: ReconstructionFixX
USE MOD_Reconstruction, ONLY: ReconstructionFixY
USE MOD_ShocksIndicator,ONLY: ShocksIndicatorX
USE MOD_ShocksIndicator,ONLY: ShocksIndicatorY
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! FORMAL ARGUMENTS
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: t
!-------------------------------------------------------------------------------!
! LOCAL VARIABLES
!-------------------------------------------------------------------------------!

CALL BoundaryConditions(t)

CALL ShocksIndicatorX()
CALL FluxSplittingX()
CALL ReconstructionX()
CALL ReconstructionFixX()
CALL TotalFluxX()

CALL ShocksIndicatorY()
CALL FluxSplittingY()
CALL ReconstructionY()
CALL ReconstructionFixY()
CALL TotalFluxY()

CALL SourceTerms(t)
CALL UpdateTimeDerivative()

!-------------------------------------------------------------------------------!
END SUBROUTINE FDTimeDerivative
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE UpdateTimeDerivative()
!-------------------------------------------------------------------------------!
USE MOD_FiniteDifference2D_vars,ONLY: Ut, S
USE MOD_FiniteDifference2D_vars,ONLY: MESH_DX
USE MOD_FiniteDifference2D_vars,ONLY: nVar
USE MOD_FiniteDifference2D_vars,ONLY: nElemsX
USE MOD_FiniteDifference2D_vars,ONLY: nElemsY
USE MOD_FiniteDifference2D_vars,ONLY: FX
USE MOD_FiniteDifference2D_vars,ONLY: FY
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! FORMAL ARGUMENTS
!-------------------------------------------------------------------------------!
! LOCAL VARIABLES
!-------------------------------------------------------------------------------!
INTEGER :: ii, jj
!-------------------------------------------------------------------------------!

DO jj=1,nElemsY
  DO ii=1,nElemsX
    Ut(1:nVar,ii,jj) = S(1:nVar,ii,jj) &
                     - (FX(1:nVar,ii+0,jj+0) - FX(1:nVar,ii-1,jj+0))/Mesh_DX(1) &
                     - (FY(1:nVar,ii+0,jj+0) - FY(1:nVar,ii+0,jj-1))/Mesh_DX(2)
  END DO
END DO

!-------------------------------------------------------------------------------!
END SUBROUTINE UpdateTimeDerivative
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE FluxSplittingX()
!-------------------------------------------------------------------------------!
USE MOD_Equation               ,ONLY: EvaluateFluxX
USE MOD_FiniteDifference2D_vars,ONLY: Fn, Fp, U, V
USE MOD_FiniteDifference2D_vars,ONLY: nVar
USE MOD_FiniteDifference2D_vars,ONLY: nElemsX
USE MOD_FiniteDifference2D_vars,ONLY: nElemsY
USE MOD_FiniteDifference2D_vars,ONLY: nGhosts
USE MOD_FiniteDifference2D_vars,ONLY: LambdaMaxX
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! FORMAL ARGUMENTS
!-------------------------------------------------------------------------------!
! LOCAL VARIABLES
!-------------------------------------------------------------------------------!
REAL    :: Cons(1:nVar)
REAL    :: Prim(1:nVar)
REAL    :: Flux(1:nVar)
INTEGER :: ii, jj
!-------------------------------------------------------------------------------!

DO jj=-nGhosts,nElemsY+nGhosts+1
  DO ii=-nGhosts,nElemsX+nGhosts+1
    Cons(1:nVar) = U(1:nVar,ii,jj)
    Prim(1:nVar) = V(1:nVar,ii,jj)
    CALL EvaluateFluxX(Prim(1:nVar),Flux(1:nVar))
    Fn(1:nVar,ii,jj) = 0.5*(Flux(1:nVar)-LambdaMaxX*Cons(1:nVar))
    Fp(1:nVar,ii,jj) = 0.5*(Flux(1:nVar)+LambdaMaxX*Cons(1:nVar))
  END DO
END DO

!-------------------------------------------------------------------------------!
END SUBROUTINE FluxSplittingX
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE FluxSplittingY()
!-------------------------------------------------------------------------------!
USE MOD_Equation               ,ONLY: EvaluateFluxY
USE MOD_FiniteDifference2D_vars,ONLY: Fn, Fp, U, V
USE MOD_FiniteDifference2D_vars,ONLY: nVar
USE MOD_FiniteDifference2D_vars,ONLY: nElemsX
USE MOD_FiniteDifference2D_vars,ONLY: nElemsY
USE MOD_FiniteDifference2D_vars,ONLY: nGhosts
USE MOD_FiniteDifference2D_vars,ONLY: LambdaMaxY
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! FORMAL ARGUMENTS
!-------------------------------------------------------------------------------!
! LOCAL VARIABLES
!-------------------------------------------------------------------------------!
REAL    :: Cons(1:nVar)
REAL    :: Prim(1:nVar)
REAL    :: Flux(1:nVar)
INTEGER :: ii, jj
!-------------------------------------------------------------------------------!

DO jj=-nGhosts,nElemsY+nGhosts+1
  DO ii=-nGhosts,nElemsX+nGhosts+1
    Cons(1:nVar) = U(1:nVar,ii,jj)
    Prim(1:nVar) = V(1:nVar,ii,jj)
    CALL EvaluateFluxY(Prim(1:nVar),Flux(1:nVar))
    Fn(1:nVar,ii,jj) = 0.5*(Flux(1:nVar)-LambdaMaxY*Cons(1:nVar))
    Fp(1:nVar,ii,jj) = 0.5*(Flux(1:nVar)+LambdaMaxY*Cons(1:nVar))
  END DO
END DO

!-------------------------------------------------------------------------------!
END SUBROUTINE FluxSplittingY
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE TotalFluxX()
!-------------------------------------------------------------------------------!
USE MOD_FiniteDifference2D_vars,ONLY: WM, WP, FX
USE MOD_FiniteDifference2D_vars,ONLY: nVar
USE MOD_FiniteDifference2D_vars,ONLY: nElemsX
USE MOD_FiniteDifference2D_vars,ONLY: nElemsY
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! FORMAL ARGUMENTS
!-------------------------------------------------------------------------------!
! LOCAL VARIABLES
!-------------------------------------------------------------------------------!
INTEGER :: ii, jj
!-------------------------------------------------------------------------------!

DO jj=1,nElemsY
  DO ii=0,nElemsX
    FX(1:nVar,ii,jj) = WM(1:nVar,ii+1,jj) + WP(1:nVar,ii,jj)
  END DO
END DO

!-------------------------------------------------------------------------------!
END SUBROUTINE TotalFluxX
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE TotalFluxY()
!-------------------------------------------------------------------------------!
USE MOD_FiniteDifference2D_vars,ONLY: WM, WP, FY
USE MOD_FiniteDifference2D_vars,ONLY: nVar
USE MOD_FiniteDifference2D_vars,ONLY: nElemsX
USE MOD_FiniteDifference2D_vars,ONLY: nElemsY
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! FORMAL ARGUMENTS
!-------------------------------------------------------------------------------!
! LOCAL VARIABLES
!-------------------------------------------------------------------------------!
INTEGER :: ii, jj
!-------------------------------------------------------------------------------!

DO jj=0,nElemsY
  DO ii=1,nElemsX
    FY(1:nVar,ii,jj) = WM(1:nVar,ii,jj+1) + WP(1:nVar,ii,jj)
  END DO
END DO

!-------------------------------------------------------------------------------!
END SUBROUTINE TotalFluxY
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE FinalizeFiniteDifference()
!-------------------------------------------------------------------------------!
USE MOD_FiniteDifference2D_vars,ONLY: MeshNodes
USE MOD_FiniteDifference2D_vars,ONLY: U
USE MOD_FiniteDifference2D_vars,ONLY: V
USE MOD_FiniteDifference2D_vars,ONLY: Ut
USE MOD_FiniteDifference2D_vars,ONLY: S
USE MOD_FiniteDifference2D_vars,ONLY: WM
USE MOD_FiniteDifference2D_vars,ONLY: WP
USE MOD_FiniteDifference2D_vars,ONLY: FX
USE MOD_FiniteDifference2D_vars,ONLY: FY
USE MOD_FiniteDifference2D_vars,ONLY: Fn
USE MOD_FiniteDifference2D_vars,ONLY: Fp
USE MOD_FiniteDifference2D_vars,ONLY: Ind
USE MOD_FiniteDifference2D_vars,ONLY: K0
USE MOD_FiniteDifference2D_vars,ONLY: K1
USE MOD_FiniteDifference2D_vars,ONLY: K2
USE MOD_FiniteDifference2D_vars,ONLY: K3
USE MOD_FiniteDifference2D_vars,ONLY: K4
USE MOD_FiniteDifference2D_vars,ONLY: K5
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! FORMAL ARGUMENTS
!-------------------------------------------------------------------------------!
! LOCAL VARIABLES
!-------------------------------------------------------------------------------!

DEALLOCATE(MeshNodes)
DEALLOCATE(U)
DEALLOCATE(V)
DEALLOCATE(Fn)
DEALLOCATE(Fp)
DEALLOCATE(WM)
DEALLOCATE(WP)
DEALLOCATE(S)
DEALLOCATE(Ut)
DEALLOCATE(FX)
DEALLOCATE(FY)
DEALLOCATE(Ind)
DEALLOCATE(K0)
DEALLOCATE(K1)
DEALLOCATE(K2)
DEALLOCATE(K3)
DEALLOCATE(K4)
DEALLOCATE(K5)

!-------------------------------------------------------------------------------!
END SUBROUTINE FinalizeFiniteDifference
!===============================================================================!
!
!
!
!===============================================================================!
END MODULE MOD_FiniteDifference2D
!-------------------------------------------------------------------------------!
