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
    
      real flux_x(ids:ide,jds:jde,ifs:ife)
      real flux_y(ids:ide,jds:jde,ifs:ife)
      real div_x (ids:ide,jds:jde,ifs:ife)
      real div_y (ids:ide,jds:jde,ifs:ife)
      
      real vorticity(ids:ide,jds:jde,ifs:ife)
      
      real lambda_x(ils:ile) ! eigenvalue along x direction
      real lambda_y(jls:jle) ! eigenvalue along y direction
      
      real ux  (ips:ipe) ! u array along x direction
      real vy  (jps:jpe) ! v array along y direction
      real phix(ips:ipe) ! phi array along x direction
      real phiy(jps:jpe) ! phi array along y direction
      
      real E (ips:ipe,jps:jpe,ifs:ife) ! E = phi + K
      real Ex(ips:ipe                ) ! E array along x direction
      real Ey(jps:jpe                ) ! E array along y direction
      
      real phiu (ips:ipe,jps:jpe,ifs:ife) ! phiu = phi * contraU
      real phiux(ips:ipe                ) ! phiu array along x direction
      
      real phiv (ips:ipe,jps:jpe,ifs:ife) ! phiv = phi * contraV
      real phivy(jps:jpe                ) ! phiv array along y direction
      
      integer i,j,iPatch
      integer P1
      
      E    = stat%phi + 0.5 * (stat%contraU * stat%u + stat%contraV * stat%v)
      phiu = stat%phi * stat%contraU
      phiv = stat%phi * stat%contraV
      
      ! calculate tend in x direction
      do iPatch = ifs, ife
        do j = jds, jde
          Ex    = E       (:,j,iPatch)
          ux    = stat%u  (:,j,iPatch)
          phiux = phiu    (:,j,iPatch)
          phix  = stat%phi(:,j,iPatch)
          
          do i = ils, ile
            P1          = pvIdx(1,i)
            lambda_x(i) = eigenvalue_x(stat%contraU(P1,j,iPatch),stat%phi(P1,j,iPatch),mesh%matrixIG(:,:,P1,j,iPatch))
          enddo
          
          call calc_tendP(flux_x(:,j,iPatch),Ex   ,ux  ,lambda_x)
          call calc_tendP(div_x (:,j,iPatch),phiux,phix,lambda_x)
        enddo
      enddo
      
      ! calculate tend in y direction
      do iPatch = ifs, ife
        do i = ids, ide
          Ey    = E       (i,:,iPatch)
          vy    = stat%v  (i,:,iPatch)
          phivy = phiv    (i,:,iPatch)
          phiy  = stat%phi(i,:,iPatch)
          
          do j = jls, jle
            P1          = pvIdx(1,j)
            lambda_y(j) = eigenvalue_y(stat%contraV(i,P1,iPatch),stat%phi(i,P1,iPatch),mesh%matrixIG(:,:,i,P1,iPatch))
          enddo
            
          call calc_tendP(flux_y(i,:,iPatch),Ey   ,vy  ,lambda_y)
          call calc_tendP(div_y (i,:,iPatch),phivy,phiy,lambda_y)
        enddo
      enddo
      
      call calc_vorticity(vorticity,stat%u,stat%v)
      
      tend%u   = -flux_x + mesh%sqrtG(ids:ide,jds:jde,:) * (vorticity + mesh%f(ids:ide,jds:jde,:)) * stat%contraV(ids:ide,jds:jde,:) - mesh%dphisdx(ids:ide,jds:jde,:)
      tend%v   = -flux_y - mesh%sqrtG(ids:ide,jds:jde,:) * (vorticity + mesh%f(ids:ide,jds:jde,:)) * stat%contraU(ids:ide,jds:jde,:) - mesh%dphisdy(ids:ide,jds:jde,:)
      tend%phi = -div_x - div_y - ( stat%phi  (ids:ide,jds:jde,:) * stat%contraU(ids:ide,jds:jde,:) * mesh%dsqrtGdx(ids:ide,jds:jde,:)   &
                                  + stat%phi  (ids:ide,jds:jde,:) * stat%contraV(ids:ide,jds:jde,:) * mesh%dsqrtGdy(ids:ide,jds:jde,:) ) &
                                  / mesh%sqrtG(ids:ide,jds:jde,:)
      
      !print*,'min/max value of u         : ', minval(stat%u   ), maxval(stat%u   )
      !print*,'min/max value of v         : ', minval(stat%v   ), maxval(stat%v   )
      !print*,'min/max value of phi       : ', minval(stat%phi ), maxval(stat%phi )
      !print*,'min/max value of stat%contraU   : ', minval(stat%contraU  ), maxval(stat%contraU  )
      !print*,'min/max value of stat%contraV   : ', minval(stat%contraV  ), maxval(stat%contraV  )
      !print*,'min/max value of flux_x    : ', minval(flux_x   ), maxval(flux_x   )
      !print*,'min/max value of flux_y    : ', minval(flux_y   ), maxval(flux_y   )
      !print*,'min/max value of div_x     : ', minval(div_x    ), maxval(div_x    )
      !print*,'min/max value of div_y     : ', minval(div_y    ), maxval(div_y    )
      !print*,'min/max value of dudt      : ', minval(tend%u   ), maxval(tend%u   )
      !print*,'min/max value of dvdt      : ', minval(tend%v   ), maxval(tend%v   )
      !print*,'min/max value of dphidt    : ', minval(tend%phi ), maxval(tend%phi )
      !print*,'min/max value of vorticity : ', minval(vorticity), maxval(vorticity)
      
    end subroutine spatial_operator
    
    subroutine calc_tendP(tendP,f,q,eigenvalue)
      real   , intent(out) :: tendP     (ids:ide)
      real   , intent(in ) :: f         (ips:ipe)
      real   , intent(in ) :: q         (ips:ipe)
      real   , intent(in ) :: eigenvalue(ils:ile)
      
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
        
        call riemann_solver(tendP(P1),fxL(iCell),fxR(iCell),qxL(iCell),qxR(iCell),eigenvalue(iCell))
        
        ! tendP(P4) = tendP(P1(iCell+1))
      enddo
      
      ! For 4th order MCV
      do iCell = 1, Nx
        P1    = pvIdx(1,iCell)
        P2    = pvIdx(2,iCell)
        P3    = pvIdx(3,iCell)
        P4    = pvIdx(4,iCell)
        
        fxVIA = (f(P4) - f(P1)) / dx
        fxxc  = 4.5 * ( f(P1) - f(P2) - f(P3) + f(P4) ) / (dx**2)
        
        tendP(P2) = 0.75 * fxVIA - (4. * tendP(P1) + 5. * tendP(P4)) / 27. - 4. / 27. * dx * fxxc
        tendP(P3) = 0.75 * fxVIA - (5. * tendP(P1) + 4. * tendP(P4)) / 27. + 4. / 27. * dx * fxxc
      enddo
      
      !tendP = -tendP
      
    end subroutine calc_tendP
    
    subroutine polyfit_tend(fitTendCL,fitTendCR,pv,dx)
      real, intent(in ) :: pv(DOF)
      real, intent(in ) :: dx
      real, intent(out) :: fitTendCL ! left side tend of the cell
      real, intent(out) :: fitTendCR ! right side tend of the cell
      
      fitTendCL = (-11. * pv(1) + 18. * pv(2) - 9. * pv(3) + 2. * pv(4))/(2. * dx);
      fitTendCR = ( 11. * pv(4) - 18. * pv(3) + 9. * pv(2) - 2. * pv(1))/(2. * dx);
    
    end subroutine polyfit_tend
    
    subroutine riemann_solver(tendP,fxL,fxR,qxL,qxR,eigenvalue)
      real, intent(in ) :: fxL
      real, intent(in ) :: fxR
      real, intent(in ) :: qxL
      real, intent(in ) :: qxR
      real, intent(in ) :: eigenvalue
      real, intent(out) :: tendP
      
      tendP = 0.5 * (fxL + fxR) - eigenvalue * (qxR - qxL)
    
    end subroutine riemann_solver
    
    function eigenvalue_x(contraU,phi,matrixIG) result(lambda)
      real contraU      ! contravariant wind on x direction
      real phi          ! geopotential height
      real lambda
      real matrixIG(1,1)
      real G11
      
      real lambda1,lambda2
      
      G11 = matrixIG(1,1)
      
      lambda1 = contraU - sqrt(G11*phi)
      lambda2 = contraU + sqrt(G11*phi)
      lambda  = max(abs(lambda1),abs(lambda2))
      
    end function eigenvalue_x
    
    function eigenvalue_y(contraV,phi,matrixIG) result(lambda)
      real contraV      ! contravariant wind on y direction
      real phi          ! geopotential height
      real lambda
      real matrixIG(2,2)
      real G22
      
      real lambda1,lambda2
      
      G22 = matrixIG(2,2)
      
      lambda1 = contraV - sqrt(G22*phi)
      lambda2 = contraV + sqrt(G22*phi)
      lambda  = max(abs(lambda1),abs(lambda2))
      
    end function eigenvalue_y
    
    subroutine calc_vorticity(vorticity,u,v)
      real, intent(out) :: vorticity(ids:ide,jds:jde,ifs:ife)
      real, intent(in ) :: u        (ips:ipe,jps:jpe,ifs:ife)
      real, intent(in ) :: v        (ips:ipe,jps:jpe,ifs:ife)
      
      real dvdx(ids:ide,jds:jde,ifs:ife)
      real dudy(ids:ide,jds:jde,ifs:ife)
      
      integer i,j,iPatch
      
      do iPatch = ifs, ife
        do j = jds, jde
          call CD4(dvdx(:,j,iPatch),v(:,j,iPatch),dx,ips,ipe,ids,ide)
        enddo
      enddo
      
      do iPatch = ifs, ife
        do i = ids, ide
          call CD4(dudy(i,:,iPatch),u(i,:,iPatch),dy,jps,jpe,jds,jde)
        enddo
      enddo
      
      vorticity = (dvdx - dudy) / mesh%sqrtG(ids:ide,jds:jde,:)
      
    end subroutine calc_vorticity
    
    subroutine CD4(dqdh,q,dh,ims,ime,its,ite)
      real   , intent(in ) :: q   (ims:ime)
      real   , intent(in ) :: dh
      integer, intent(in ) :: ims,ime,its,ite
      real   , intent(out) :: dqdh(its:ite)
      
      integer i
      
      do i = ids, ide
        dqdh(i) = (q(i-2) - 8. * q(i-1) + 8. * q(i+1) - q(i+2)) / 12. / dh
      enddo
    
    end subroutine CD4
  
    subroutine unify_bdy_stat(stat)
      type(stat_field), intent(inout) :: stat
      
      call unify_bdy_field(stat%phi            (ids:ide,jds:jde,:))
      call unify_bdy_field(stat%zonal_wind     (ids:ide,jds:jde,:))
      call unify_bdy_field(stat%meridional_wind(ids:ide,jds:jde,:))
      
    end subroutine unify_bdy_stat
    
    subroutine unify_bdy_field(field)
      real ,intent(inout) :: field(ids:ide,jds:jde,ifs:ife)
      
      ! Bdy lines
      field(ids,jds:jde,1) = 0.5 * (field(ids,jds:jde,4) + field(ide,jds:jde,1)) ! left bdy of patch 1
      field(ids,jds:jde,2) = 0.5 * (field(ids,jds:jde,1) + field(ide,jds:jde,2)) ! left bdy of patch 2
      field(ids,jds:jde,3) = 0.5 * (field(ids,jds:jde,2) + field(ide,jds:jde,3)) ! left bdy of patch 3
      field(ids,jds:jde,4) = 0.5 * (field(ids,jds:jde,3) + field(ide,jds:jde,4)) ! left bdy of patch 4
      
      field(ids:ide,jde,1) = 0.5 * (field(ids:ide,jds,5) + field(ids:ide,jde,1)) ! top bdy of patch 1
      field(ids:ide,jds,1) = 0.5 * (field(ids:ide,jds,1) + field(ids:ide,jde,6)) ! bottom bdy of patch 1
      
      field(ids:ide,jde,2) = 0.5 * (field(ide,jds:jde,5) + field(ids:ide,jde   ,2)) ! top bdy of patch 2
      field(ids:ide,jds,2) = 0.5 * (field(ids:ide,jds,2) + field(ide,jde:jds:-1,6)) ! bottom bdy of patch 2
      
      field(ids:ide,jde,3) = 0.5 * (field(ids:ide,jde,3) + field(ide:ids:-1,jde,5)) ! top bdy of patch 3
      field(ids:ide,jds,3) = 0.5 * (field(ids:ide,jds,3) + field(ide:ids:-1,jds,6)) ! bottom bdy of patch 3
      
      field(ids:ide,jde,4) = 0.5 * (field(ids:ide,jde,4) + field(ids,jde:jds:-1,5)) ! top bdy of patch 4
      field(ids:ide,jds,4) = 0.5 * (field(ids:ide,jds,4) + field(ids,jds:jde   ,6)) ! bottom bdy of patch 4
      
      field(ide,jds:jde,1) = field(ids,jds:jde,2) ! right bdy of patch 1
      field(ide,jds:jde,2) = field(ids,jds:jde,3) ! right bdy of patch 2
      field(ide,jds:jde,3) = field(ids,jds:jde,4) ! right bdy of patch 3
      field(ide,jds:jde,4) = field(ids,jds:jde,1) ! right bdy of patch 4
      
      field(ids,jds:jde,5) = field(ide:ids:-1,jde,4) ! left bdy of patch 5
      field(ids,jds:jde,6) = field(ids:ide   ,jds,4) ! left bdy of patch 6
      
      field(ide,jds:jde,5) = field(ids:ide,jde   ,2) ! right bdy of 5
      field(ide,jds:jde,6) = field(ide:ids:-1,jds,2) ! right bdy of 6
      
      field(ids:ide,jde,5) = field(ide:ids:-1,jde,3) ! top bdy of patch 5
      field(ids:ide,jds,5) = field(ids:ide   ,jde,1) ! bottom bdy of patch 5
      
      field(ids:ide,jde,6) = field(ids:ide   ,jds,1) ! top bdy of patch 6
      field(ids:ide,jds,6) = field(ide:ids:-1,jds,3) ! bottom bdy of patch 6
      
      ! Bdy points
      field(ids,jds,1) = (field(ids,jds,1) + field(ide,jds,4) + field(ids,jde,6)) / 3. ! low-left point of patch 1
      field(ide,jds,1) = (field(ide,jds,1) + field(ids,jds,2) + field(ide,jde,6)) / 3. ! low-right point of patch 1
      field(ide,jde,1) = (field(ide,jde,1) + field(ids,jde,2) + field(ide,jds,5)) / 3. ! up-right point of patch 1
      field(ids,jde,1) = (field(ids,jde,1) + field(ide,jde,4) + field(ids,jds,5)) / 3. ! up-left point of patch 1
      
      field(ids,jds,3) = (field(ids,jds,3) + field(ide,jds,2) + field(ide,jds,6)) / 3. ! low-left point of patch 3
      field(ide,jds,3) = (field(ide,jds,3) + field(ids,jds,4) + field(ide,jds,6)) / 3. ! low-right point of patch 3
      field(ide,jde,3) = (field(ide,jde,3) + field(ids,jde,5) + field(ids,jde,4)) / 3. ! up-right point of patch 3
      field(ids,jde,3) = (field(ids,jde,3) + field(ide,jde,2) + field(ide,jds,5)) / 3. ! up-left point of patch 3
      
      field(ids,jds,2) = field(ide,jds,1) ! low-left point of patch 2
      field(ide,jds,2) = field(ids,jds,3) ! low-right point of patch 2
      field(ide,jde,2) = field(ids,jde,3) ! up-right point of patch 2
      field(ids,jde,2) = field(ide,jde,1) ! up-left point of patch 2

      field(ids,jds,4) = field(ide,jds,3) ! low-left point of patch 4
      field(ide,jds,4) = field(ids,jds,1) ! low-right point of patch 4
      field(ide,jde,4) = field(ids,jde,1) ! up-right point of patch 4
      field(ids,jde,4) = field(ide,jde,3) ! up-left point of patch 4
      
      field(ids,jds,5) = field(ids,jde,1) ! low-left point of patch 5
      field(ide,jds,5) = field(ide,jde,1) ! low-right point of patch 5
      field(ide,jde,5) = field(ids,jde,3) ! up-right point of patch 5
      field(ids,jde,5) = field(ide,jde,3) ! up-left point of patch 5

      field(ids,jds,6) = field(ide,jds,3) ! low-left point of patch 6
      field(ide,jds,6) = field(ids,jds,3) ! low-right point of patch 6
      field(ide,jde,6) = field(ide,jds,1) ! up-right point of patch 6
      field(ids,jde,6) = field(ids,jds,1) ! up-left point of patch 6
      
    end subroutine unify_bdy_field
    
    ! convert vector from patch to sphere
    subroutine convert_wind_P2SP(stat)
      type(stat_field), intent(inout) :: stat
      
      integer i,j,iPatch
      
      do iPatch = ifs, ife
        do j = jps, jpe
          do i = ips, ipe
            call covProjPlane2Sphere(stat%zonal_wind(i,j,iPatch), stat%meridional_wind(i,j,iPatch), stat%u(i,j,iPatch), stat%v(i,j,iPatch), mesh%matrixA(:,:,i,j,iPatch), mesh%matrixIG(:,:,i,j,iPatch))
          enddo
        enddo
      enddo
    
    end subroutine convert_wind_P2SP
    
    subroutine convert_wind_SP2P(stat)
      type(stat_field), intent(inout) :: stat
      
      integer i,j,iPatch
      
      do iPatch = ifs, ife
        do j = jps, jpe
          do i = ips, ipe
            call covProjSphere2Plane(stat%u(i,j,iPatch), stat%v(i,j,iPatch), stat%zonal_wind(i,j,iPatch), stat%meridional_wind(i,j,iPatch), mesh%matrixIA(:,:,i,j,iPatch), mesh%matrixG(:,:,i,j,iPatch))
          enddo
        enddo
      enddo
    
    end subroutine convert_wind_SP2P
    
    subroutine convert_wind_cov2contrav(stat)
      type(stat_field), intent(inout) :: stat
      
      integer i,j,iPatch
      
      do iPatch = ifs, ife
        do j = jps, jpe
          do i = ips, ipe
            call cov2contrav(stat%contraU(i,j,iPatch),stat%contraV(i,j,iPatch),stat%u(i,j,iPatch),stat%v(i,j,iPatch),mesh%matrixIG(:,:,i,j,iPatch))
          enddo
        enddo
      enddo
    end subroutine convert_wind_cov2contrav
END MODULE spatial_operators_mod

