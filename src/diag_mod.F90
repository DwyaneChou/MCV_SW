MODULE diag_mod
  use constants_mod
  use mesh_mod
  use parameters_mod
  use stat_mod
  
  implicit none
  
    contains
    subroutine calc_total_mass(total_mass,stat)
      type(stat_field), intent(in ) :: stat
      real            , intent(out) :: total_mass
      
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
                jPV = pvIdy(jDOF,jCell)
                
                massOnCell(iCell,jCell,iPatch) = massOnCell(iCell,jCell,iPatch) + mesh%weightsOnPV(iDOF,jDOF) * stat%phiG(iPV,jPV,iPatch)
              enddo
            enddo
          enddo
        enddo
      enddo
      
      total_mass = sum(massOnCell)
    
    end subroutine calc_total_mass
    
    subroutine calc_total_energy(total_energy,stat)
      type(stat_field), intent(in ) :: stat
      real            , intent(out) :: total_energy
      
      real energyOnCell(Nx,Ny,Nf)
      real KE,PE
      integer iCell,jCell,iPatch,iPV,jPV,iDOF,jDOF
      
      energyOnCell = 0.
      total_energy = 0.
      do iPatch = ifs, ife
        do jCell = 1, Ny
          do iCell = 1, Nx
            do jDOF = 1, DOF
              do iDOF = 1, DOF
                iPV = pvIdx(iDOF,iCell)
                jPV = pvIdy(jDOF,jCell)
                
                KE = 0.5 * stat%phi(iPV,jPV,iPatch) * (stat%u(iPV,jPV,iPatch) * stat%contraU(iPV,jPV,iPatch) + stat%v(iPV,jPV,iPatch) * stat%contraV(iPV,jPV,iPatch))
                PE = 0.5 * (stat%phi(iPV,jPV,iPatch) + mesh%phi_s(iPV,jPV,iPatch))**2
                
                energyOnCell(iCell,jCell,iPatch) = energyOnCell(iCell,jCell,iPatch) + mesh%weightsOnPV(iDOF,jDOF) * mesh%sqrtG(iPV,jPV,iPatch) * (KE + PE)
              enddo
            enddo
          enddo
        enddo
      enddo
      
      total_energy = sum(energyOnCell)
    end subroutine calc_total_energy
    
    subroutine calc_VIA(fieldOnCell,fieldOnPoint)
      real, intent(out) :: fieldOnCell (ics:ice,jcs:jce,ifs:ife)
      real, intent(in ) :: fieldOnPoint(ips:ipe,jps:jpe,ifs:ife)
      
      integer iCell,jCell,iPatch,iPV,jPV,iDOF,jDOF
      
      fieldOnCell = 0.
      
      do iPatch = ifs, ife
        do jCell = 1, Ny
          do iCell = 1, Nx
            do jDOF = 1, DOF
              do iDOF = 1, DOF
                iPV = pvIdx(iDOF,iCell)
                jPV = pvIdy(jDOF,jCell)
                
                fieldOnCell(iCell,jCell,iPatch) = fieldOnCell(iCell,jCell,iPatch) + mesh%weightsOnPV(iDOF,jDOF) * fieldOnPoint(iPV,jPV,iPatch)
              enddo
            enddo
          enddo
        enddo
      enddo
      
    
    end subroutine calc_VIA
END MODULE diag_mod

