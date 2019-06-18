MODULE spatial_operators_mod
  use mesh_mod
  use parameters_mod
  use stat_mod
  use tend_mod
  use projection_mod
  
  implicit none
  
    contains
    
    subroutine spatial_operator(stat,tend)
      type(stat_field), intent(in ) :: stat
      type(tend_field), intent(out) :: tend
    
      real div_x (ids:ide,jds:jde,ifs:ife)
      real div_y (ids:ide,jds:jde,ifs:ife)
      
      real upstream_x(ips:ipe) ! upstream coefficient along x direction
      real upstream_y(jps:jpe) ! upstream coefficient along y direction
      
      real ux   (ips:ipe) ! u array along x direction
      real vy   (jps:jpe) ! u array along x direction
      
      real phiGu (ips:ipe,jps:jpe,ifs:ife) ! phiGu = phi * contraU
      real phiGux(ips:ipe                ) ! phiGu array along x direction
      
      real phiGv (ips:ipe,jps:jpe,ifs:ife) ! phiGv = phi * contraV
      real phiGvy(jps:jpe                ) ! phiGv array along y direction
      
      integer i,j,iPatch
      integer P1
      
      phiGu = stat%phiG * stat%contraU
      phiGv = stat%phiG * stat%contraV
      
      !!$OMP PARALLEL SECTIONS PRIVATE(iPatch)
      !!$OMP SECTION
        ! calculate tend in x direction
        !!$OMP PARALLEL DO PRIVATE(i,j,phiGux,upstream_x)
        do iPatch = ifs, ife
          do j = jds, jde
            phiGux= phiGu    (:,j,iPatch)
            
            where(ux==0)ux=1
            upstream_x = ux/abs(ux)
            
            call calc_tendP(div_x (:,j,iPatch),phiGux,phiGux,upstream_x)
          enddo
        enddo
        !!$OMP END PARALLEL DO
      
      !!$OMP SECTION
        ! calculate tend in y direction
        !!$OMP PARALLEL DO PRIVATE(i,j,phiGvy,upstream_y)
        do iPatch = ifs, ife
          do i = ids, ide
            phiGvy = phiGv    (i,:,iPatch)
            
            where(vy==0)vy=1
            upstream_y = vy/abs(vy)
            
            call calc_tendP(div_y (i,:,iPatch),phiGvy,phiGvy,upstream_y)
          enddo
        enddo
        !!$OMP END PARALLEL DO
      !!$OMP END PARALLEL SECTIONS
      
      tend%phiG = -div_x - div_y
      
      !print*,'min/max value of u              : ', minval(stat%u      ), maxval(stat%u      )
      !print*,'min/max value of v              : ', minval(stat%v      ), maxval(stat%v      )
      !print*,'min/max value of phi            : ', minval(stat%phi    ), maxval(stat%phi    )
      !print*,'min/max value of stat%contraU   : ', minval(stat%contraU), maxval(stat%contraU)
      !print*,'min/max value of stat%contraV   : ', minval(stat%contraV), maxval(stat%contraV)
      !print*,'min/max value of flux_x         : ', minval(flux_x      ), maxval(flux_x      )
      !print*,'min/max value of flux_y         : ', minval(flux_y      ), maxval(flux_y      )
      !print*,'min/max value of div_x          : ', minval(div_x       ), maxval(div_x       )
      !print*,'min/max value of div_y          : ', minval(div_y       ), maxval(div_y       )
      !print*,'min/max value of dudt           : ', minval(tend%u      ), maxval(tend%u      )
      !print*,'min/max value of dvdt           : ', minval(tend%v      ), maxval(tend%v      )
      !print*,'min/max value of dphidt         : ', minval(tend%phiG   ), maxval(tend%phiG   )
      !print*,'min/max value of dudy           : ', minval(dudy        ), maxval(dudy        )
      !print*,'min/max value of dvdx           : ', minval(dvdx        ), maxval(dvdx        )
      !print*,'min/max value of vorticity      : ', minval(vorticity   ), maxval(vorticity   )
      
    end subroutine spatial_operator
    
    subroutine calc_tendP(derivP,f,q,eigenvalue)
      real   , intent(out) :: derivP    (ids:ide)
      real   , intent(in ) :: f         (ips:ipe)
      real   , intent(in ) :: q         (ips:ipe)
      real   , intent(in ) :: eigenvalue(ips:ipe)
      
      real fxL    (ics:ice)
      real fxR    (ics:ice)
      real qxL    (ics:ice)
      real qxR    (ics:ice)
      real fxL_fit(ics:ice)
      real fxR_fit(ics:ice)
      real qxL_fit(ics:ice)
      real qxR_fit(ics:ice)
      real fCell  (DOF)
      real qCell  (DOF)
      
      real fxVIA       ! The 1st order derivative of f on cell
      real fxxc        ! The 2nd order derivative of f
      
      real lambda_max  ! max eigenvalue
      
      integer P1,P2,P3,P4 ! points in cell
      
      integer iCell
      
      do iCell = ics, ice
        fCell = f(pvIdx(:,iCell))
        qCell = q(pvIdx(:,iCell))
        
        call polyfit_tend(fxL_fit(iCell),fxR_fit(iCell),fCell,dx)
        call polyfit_tend(qxL_fit(iCell),qxR_fit(iCell),qCell,dx)
      enddo
      
      do iCell = ils, ile
        fxL(iCell) = fxR_fit(iCell-1)
        fxR(iCell) = fxL_fit(iCell)
        qxL(iCell) = qxR_fit(iCell-1)
        qxR(iCell) = qxL_fit(iCell)
        
        P1         = pvIdx(1,iCell)
        lambda_max = eigenvalue(P1)
        !lambda_max = maxval(eigenvalue(P1-DOF+1:P1+DOF-1))
        !lambda_max = maxval(eigenvalue)
        call riemann_solver(derivP(P1),fxL(iCell),fxR(iCell),qxL(iCell),qxR(iCell),lambda_max)
        
        ! derivP(P4) = derivP(P1(iCell+1))
      enddo
      
      ! Compute tend of inner point(s) on cells
#    ifdef MCV3
      ! For MCV3 only
      do iCell = 1, Nx
        P1    = pvIdx(1,iCell)
        P2    = pvIdx(2,iCell)
        P3    = pvIdx(3,iCell)
        
        fxVIA     = (f(P3) - f(P1)) / dx
        derivP(P2) = 1.5 * fxVIA - 0.25 * (derivP(P1) + derivP(P3))
      enddo
#    endif
#    ifdef MCV4
      ! For MCV4 only
      do iCell = 1, Nx
        P1    = pvIdx(1,iCell)
        P2    = pvIdx(2,iCell)
        P3    = pvIdx(3,iCell)
        P4    = pvIdx(4,iCell)
        
        fxVIA = (f(P4) - f(P1)) / dx
        fxxc  = 4.5 * ( f(P1) - f(P2) - f(P3) + f(P4) ) / (dx**2)
        
        derivP(P2) = 4./3. * fxVIA - (4. * derivP(P1) + 5. * derivP(P4)) / 27. - 4. / 27. * dx * fxxc
        derivP(P3) = 4./3. * fxVIA - (5. * derivP(P1) + 4. * derivP(P4)) / 27. + 4. / 27. * dx * fxxc
      enddo
#    endif
      
      !derivP = -derivP
      
    end subroutine calc_tendP
    
    subroutine polyfit_tend(fitDerivCL,fitDerivCR,pv,dh)
      real, intent(in ) :: pv(DOF)
      real, intent(in ) :: dh
      real, intent(out) :: fitDerivCL ! left side tend of the cell
      real, intent(out) :: fitDerivCR ! right side tend of the cell
      
#    ifdef MCV3
      ! For MCV3 only
      fitDerivCL = -( pv(3) - 4. * pv(2) + 3. * pv(1)) / dh
      fitDerivCR =  ( pv(1) - 4. * pv(2) + 3. * pv(3)) / dh
#    endif

#    ifdef MCV4
      ! For MCV4 only
      fitDerivCL = (-11. * pv(1) + 18. * pv(2) - 9. * pv(3) + 2. * pv(4))/(2. * dh);
      fitDerivCR = ( 11. * pv(4) - 18. * pv(3) + 9. * pv(2) - 2. * pv(1))/(2. * dh);
#    endif
    end subroutine polyfit_tend
    
    subroutine riemann_solver(derivP,fxL,fxR,qxL,qxR,eigenvalue)
      real, intent(in ) :: fxL
      real, intent(in ) :: fxR
      real, intent(in ) :: qxL
      real, intent(in ) :: qxR
      real, intent(in ) :: eigenvalue
      real, intent(out) :: derivP
      
      derivP = 0.5 * ((fxL + fxR) - eigenvalue * (qxR - qxL))
    
    end subroutine riemann_solver
    
    subroutine unify_bdy_stat(stat)
      type(stat_field), intent(inout) :: stat
      
      call unify_bdy_field(stat%phi            (ids:ide,jds:jde,ifs:ife))
      
    end subroutine unify_bdy_stat
    
    subroutine unify_bdy_field(field)
      real,intent(inout) :: field(ids:ide,jds:jde,ifs:ife)
      
      real field_raw(ids:ide,jds:jde,ifs:ife)
      
      field_raw = field
      
      ! Bdy lines
      ! edge(1,4)
      field(ids,jds:jde,1) = 0.5 * (field_raw(ide,jds:jde,4) + field_raw(ids,jds:jde,1)) ! left bdy of patch 1
      field(ide,jds:jde,4) = field(ids,jds:jde,1)
      
      ! edge(1,2)
      field(ids,jds:jde,2) = 0.5 * (field_raw(ide,jds:jde,1) + field_raw(ids,jds:jde,2)) ! left bdy of patch 2
      field(ide,jds:jde,1) = field(ids,jds:jde,2)
      
      ! edge(2,3)
      field(ids,jds:jde,3) = 0.5 * (field_raw(ide,jds:jde,2) + field_raw(ids,jds:jde,3)) ! left bdy of patch 3
      field(ide,jds:jde,2) = field(ids,jds:jde,3)
      
      ! edge(3,4)
      field(ids,jds:jde,4) = 0.5 * (field_raw(ide,jds:jde,3) + field_raw(ids,jds:jde,4)) ! left bdy of patch 4
      field(ide,jds:jde,3) = field(ids,jds:jde,4)
      
      ! edge(1,5)
      field(ids:ide,jde,1) = 0.5 * (field_raw(ids:ide,jds,5) + field_raw(ids:ide,jde,1)) ! top bdy of patch 1
      field(ids:ide,jds,5) = field(ids:ide,jde,1)
      
      ! edge(1,6)
      field(ids:ide,jds,1) = 0.5 * (field_raw(ids:ide,jds,1) + field_raw(ids:ide,jde,6)) ! bottom bdy of patch 1
      field(ids:ide,jde,6) = field(ids:ide,jds,1)
      
      ! edge(2,5)
      field(ids:ide,jde,2) = 0.5 * (field_raw(ide,jds:jde,5) + field_raw(ids:ide,jde,2)) ! top bdy of patch 2
      field(ide,jds:jde,5) = field(ids:ide,jde,2)
      
      ! edge(2,6)
      field(ids:ide,jds   ,2) = 0.5 * (field_raw(ids:ide,jds,2) + field_raw(ide,jde:jds:-1,6)) ! bottom bdy of patch 2
      field(ide,jde:jds:-1,6) = field(ids:ide,jds,2)
      
      ! edge(3,5)
      field(ids:ide,jde   ,3) = 0.5 * (field_raw(ids:ide,jde,3) + field_raw(ide:ids:-1,jde,5)) ! top bdy of patch 3
      field(ide:ids:-1,jde,5) = field(ids:ide,jde,3)
      
      ! edge(3,6)
      field(ids:ide,jds   ,3) = 0.5 * (field_raw(ids:ide,jds,3) + field_raw(ide:ids:-1,jds,6)) ! bottom bdy of patch 3
      field(ide:ids:-1,jds,6) = field(ids:ide,jds,3)
      
      ! edge(4,5)
      field(ids:ide,jde   ,4) = 0.5 * (field_raw(ids:ide,jde,4) + field_raw(ids,jde:jds:-1,5)) ! top bdy of patch 4
      field(ids,jde:jds:-1,5) = field(ids:ide,jde,4)
      
      ! edge(4,6)
      field(ids:ide,jds,4) = 0.5 * (field_raw(ids:ide,jds,4) + field_raw(ids,jds:jde,6)) ! bottom bdy of patch 4
      field(ids,jds:jde,6) = field(ids:ide,jds,4)
      
      ! Vertex points
      ! Point (1,4,5)
      field(ids,jde,1) = (field_raw(ids,jde,1) + field_raw(ide,jde,4) + field_raw(ids,jds,5)) / 3. ! up-left point of patch 1
      field(ids,jds,5) = field(ids,jde,1)
      field(ide,jde,4) = field(ids,jde,1)
      
      ! Point(1,2,5)
      field(ide,jde,1) = (field_raw(ide,jde,1) + field_raw(ids,jde,2) + field_raw(ide,jds,5)) / 3. ! up-right point of patch 1
      field(ids,jde,2) = field(ide,jde,1)
      field(ide,jds,5) = field(ide,jde,1)
      
      ! Point(1,4,6)
      field(ids,jds,1) = (field_raw(ids,jds,1) + field_raw(ide,jds,4) + field_raw(ids,jde,6)) / 3. ! low-left point of patch 1
      field(ide,jds,4) = field(ids,jds,1)
      field(ids,jde,6) = field(ids,jds,1)
      
      ! Point(1,2,6)
      field(ide,jds,1) = (field_raw(ide,jds,1) + field_raw(ids,jds,2) + field_raw(ide,jde,6)) / 3. ! low-right point of patch 1
      field(ids,jds,2) = field(ide,jds,1)
      field(ide,jde,6) = field(ide,jds,1)
      
      ! Point(2,3,5)
      field(ids,jde,3) = (field_raw(ids,jde,3) + field_raw(ide,jde,2) + field_raw(ide,jde,5)) / 3. ! up-left point of patch 3
      field(ide,jde,2) = field(ids,jde,3)
      field(ide,jde,5) = field(ids,jde,3)
      
      ! Point(3,4,5)
      field(ide,jde,3) = (field_raw(ide,jde,3) + field_raw(ids,jde,4) + field_raw(ids,jde,5)) / 3. ! up-right point of patch 3
      field(ids,jde,4) = field(ide,jde,3)
      field(ids,jde,5) = field(ide,jde,3)
      
      ! Point(2,3,6)
      field(ids,jds,3) = (field_raw(ide,jds,2) + field_raw(ids,jds,3) + field_raw(ide,jds,6)) / 3. ! low-left point of patch 3
      field(ide,jds,2) = field(ids,jds,3)
      field(ide,jds,6) = field(ids,jds,3)
      
      ! Point(3,4,6)
      field(ide,jds,3) = (field_raw(ide,jds,3) + field_raw(ids,jds,4) + field_raw(ids,jds,6)) / 3. ! low-right point of patch 3
      field(ids,jds,4) = field(ide,jds,3)
      field(ids,jds,6) = field(ide,jds,3)
      
    end subroutine unify_bdy_field
    
    ! convert vector from patch to sphere
    subroutine convert_wind_P2SP(stat)
      type(stat_field), intent(inout) :: stat
      
      integer i,j,iPatch
      
      !!$OMP PARALLEL DO PRIVATE(i,j)
      do iPatch = ifs, ife
        do j = jps, jpe
          do i = ips, ipe
            call covProjPlane2Sphere(stat%zonal_wind(i,j,iPatch), stat%meridional_wind(i,j,iPatch), stat%u(i,j,iPatch), stat%v(i,j,iPatch), mesh%matrixA(:,:,i,j,iPatch), mesh%matrixIG(:,:,i,j,iPatch))
          enddo
        enddo
      enddo
      !!$OMP END PARALLEL DO
    end subroutine convert_wind_P2SP
    
    subroutine convert_wind_SP2P(stat)
      type(stat_field), intent(inout) :: stat
      
      integer i,j,iPatch
      
      !!$OMP PARALLEL DO PRIVATE(i,j)
      do iPatch = ifs, ife
        do j = jps, jpe
          do i = ips, ipe
            call covProjSphere2Plane(stat%u(i,j,iPatch), stat%v(i,j,iPatch), stat%zonal_wind(i,j,iPatch), stat%meridional_wind(i,j,iPatch), mesh%matrixIA(:,:,i,j,iPatch), mesh%matrixG(:,:,i,j,iPatch))
          enddo
        enddo
      enddo
      !!$OMP END PARALLEL DO
    end subroutine convert_wind_SP2P
    
    subroutine convert_wind_cov2contrav(stat)
      type(stat_field), intent(inout) :: stat
      
      integer i,j,iPatch
      
      !!$OMP PARALLEL DO PRIVATE(i,j)
      do iPatch = ifs, ife
        do j = jps, jpe
          do i = ips, ipe
            call cov2contrav(stat%contraU(i,j,iPatch),stat%contraV(i,j,iPatch),stat%u(i,j,iPatch),stat%v(i,j,iPatch),mesh%matrixIG(:,:,i,j,iPatch))
          enddo
        enddo
      enddo
      !!$OMP END PARALLEL DO
    end subroutine convert_wind_cov2contrav
    
END MODULE spatial_operators_mod

