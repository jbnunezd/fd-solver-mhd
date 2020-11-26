!===============================================================================!
PROGRAM FiniteDifference2D
!-------------------------------------------------------------------------------!
USE MOD_Mesh,              ONLY: BuildMesh
USE MOD_Parameters,        ONLY: InitializeParameters
USE MOD_FiniteDifference2D,ONLY: FillInitialConditions
USE MOD_FiniteDifference2D,ONLY: FinalizeFiniteDifference
USE MOD_FiniteDifference2D,ONLY: InitializeFiniteDifference
USE MOD_TimeDiscretization,ONLY: TimeDiscretization
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!

CALL InitializeParameters()
CALL InitializeFiniteDifference()
CALL BuildMesh()
CALL FillInitialConditions()
CALL TimeDiscretization()
CALL FinalizeFiniteDifference()

!-------------------------------------------------------------------------------!
END PROGRAM FiniteDifference2D
!===============================================================================!
