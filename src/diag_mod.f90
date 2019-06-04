MODULE diag_mod
  use constants_mod
  use mesh_mod
  use parameters_mod
  use stat_mod
  
  implicit none
  
    contains
    subroutine calc_total_mass(total_mass,stat)
      type(stat_field), intent(in ) :: stat
      real           , intent(out) :: total_mass
      
      real    massOnCell(Nx,Ny,Nf)
      integer iCell,jCell,iPatch,iPV,jPV,iDOF,jDOF
      
      massOnCell = 0
      total_mass = 0.
      do iPatch = ifs, ife
        do jCell = 1, Ny
          do iCell = 1, Nx
            do jDOF = 1, DOF
              do iDOF = 1, DOF
                iPV = pvIdx(iDOF,iCell)
                jPV = pvIdx(jDOF,jCell)
                
                massOnCell(iCell,jCell,iPatch) = massOnCell(iCell,jCell,iPatch) + mesh%weightsOnPV(iDOF,jDOF) * stat%phi(iPV,jPV,iPatch)
              enddo
            enddo
          enddo
        enddo
        massOnCell(:,:,iPatch) = massOnCell(:,:,iPatch) * mesh%areaCell
      enddo
      
      total_mass = sum(massOnCell)
    
    end subroutine calc_total_mass
    

END MODULE diag_mod

