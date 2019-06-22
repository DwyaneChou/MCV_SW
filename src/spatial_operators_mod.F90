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
      real dvdx  (ids:ide,jds:jde,ifs:ife)
      real dudy  (ids:ide,jds:jde,ifs:ife)
      
      real vorticity(ids:ide,jds:jde,ifs:ife)
      
      real lambdaL   ! eigenvalue on the left side of cell
      real lambdaR   ! eigenvalue on the right side of cell
      real upstreamL ! upstream coefficient on the left side of cell
      real upstreamR ! upstream coefficient on the right side of cell
      
      real ux      (3*DOF-2) ! u array along x direction
      real uy      (3*DOF-2) ! u array along y direction
      real vx      (3*DOF-2) ! v array along x direction
      real vy      (3*DOF-2) ! v array along y direction
      real contraUx(3*DOF-2) ! contravariant u array along x direction
      real contraVy(3*DOF-2) ! contravariant v array along x direction
      real phiGx   (3*DOF-2) ! phiG array along x direction
      real phiGy   (3*DOF-2) ! phiG array along y direction
      
      real K (ips:ipe,jps:jpe,ifs:ife) ! K = 0.5 * (contraU*coU + contraV*coV)
      real E (ips:ipe,jps:jpe,ifs:ife) ! E = phi + phi_s + K
      real Ex(3*DOF-2                ) ! E array along x direction
      real Ey(3*DOF-2                ) ! E array along y direction
      
      real phiGu (ips:ipe,jps:jpe,ifs:ife) ! phiGu = phi * contraU
      real phiGux(3*DOF-2                ) ! phiGu array along x direction
      
      real phiGv (ips:ipe,jps:jpe,ifs:ife) ! phiGv = phi * contraV
      real phiGvy(3*DOF-2                ) ! phiGv array along y direction
      
      integer Ps,Pe ! Index of first point on left cell, index of last point on right cell
      integer Pl,Pr ! Index of left point on the center cell, index of right point on the center cell
      integer i,j,iPatch
      
      ! Ps            Pl           Pr             Pe
      ! |------o------|------o------|------o------|
      !    left cell    center cell   right cell
      
      phiGu = stat%phiG * stat%contraU
      phiGv = stat%phiG * stat%contraV
      K     = 0.5 * (stat%contraU * stat%u + stat%contraV * stat%v)
      E     = stat%phi + mesh%phi_s + K
      
      ! calculate tend in x direction
      do iPatch = ifs, ife
        !$OMP PARALLEL DO PRIVATE(i,Ex,ux,vx,contraUx,phiGux,PhiGx,Pl,Pr,Ps,Pe, lambdaL, lambdaR, upstreamL, upstreamR)
        do j = jds, jde
          do i = its, ite
            Pl       = pvIdx(1  ,i)
            Pr       = pvIdx(DOF,i)
            Ps       = Pl - DOF + 1
            Pe       = Pr + DOF - 1
            
            Ex       = E           (Ps:Pe,j,iPatch)
            ux       = stat%u      (Ps:Pe,j,iPatch)
            vx       = stat%v      (Ps:Pe,j,iPatch)
            contraUx = stat%contraU(Ps:Pe,j,iPatch)
            phiGux   = phiGu       (Ps:Pe,j,iPatch)
            phiGx    = stat%phiG   (Ps:Pe,j,iPatch)
            
            lambdaL = eigenvalue_x(stat%contraU(Pl,j,iPatch),stat%phi(Pl,j,iPatch),mesh%matrixIG(:,:,Pl,j,iPatch))
            lambdaR = eigenvalue_x(stat%contraU(Pr,j,iPatch),stat%phi(Pr,j,iPatch),mesh%matrixIG(:,:,Pr,j,iPatch))
            
            where(contraUx==0)contraUx=1.
            upstreamL = contraUx(DOF    ) / abs(contraUx(DOF    ))
            upstreamR = contraUx(2*DOF-1) / abs(contraUx(2*DOF-1))
            
            call calc_tendP(flux_x(Pl:Pr,j,iPatch),Ex    ,ux   ,lambdaL  , lambdaR  )
            call calc_tendP(div_x (Pl:Pr,j,iPatch),phiGux,phiGx,lambdaL  , lambdaR  )
            call calc_tendP(dvdx  (Pl:Pr,j,iPatch),vx    ,vx   ,upstreamL, upstreamR)
          enddo
        enddo
        !$OMP END PARALLEL DO
      enddo
      
      ! calculate tend in y direction
      do iPatch = ifs, ife
        !$OMP PARALLEL DO PRIVATE(j,Ey,uy,vy,contraVy,phiGvy,PhiGy,Pl,Pr,Ps,Pe, lambdaL, lambdaR, upstreamL, upstreamR)
        do i = ids, ide
          do j = jts, jte
            Pl       = pvIdy(1  ,j)
            Pr       = pvIdy(DOF,j)
            Ps       = Pl - DOF + 1
            Pe       = Pr + DOF - 1
            
            Ey       = E           (i,Ps:Pe,iPatch)
            uy       = stat%u      (i,Ps:Pe,iPatch)
            vy       = stat%v      (i,Ps:Pe,iPatch)
            contraVy = stat%contraV(i,Ps:Pe,iPatch)
            phiGvy   = phiGv       (i,Ps:Pe,iPatch)
            phiGy    = stat%phiG   (i,Ps:Pe,iPatch)
          
            lambdaL = eigenvalue_y(stat%contraV(i,Pl,iPatch),stat%phi(i,Pl,iPatch),mesh%matrixIG(:,:,i,Pl,iPatch))
            lambdaR = eigenvalue_y(stat%contraV(i,Pr,iPatch),stat%phi(i,Pr,iPatch),mesh%matrixIG(:,:,i,Pr,iPatch))
          
            where(contraVy==0)contraVy=1.
            upstreamL = contraVy(DOF    ) / abs(contraVy(DOF    ))
            upstreamR = contraVy(2*DOF-1) / abs(contraVy(2*DOF-1))
          
            call calc_tendP(flux_y(i,Pl:Pr,iPatch),Ey    ,vy   ,lambdaL  , lambdaR  )
            call calc_tendP(div_y (i,Pl:Pr,iPatch),phiGvy,phiGy,lambdaL  , lambdaR  )
            call calc_tendP(dudy  (i,Pl:Pr,iPatch),uy    ,uy   ,upstreamL, upstreamR)
          enddo
        enddo
        !$OMP END PARALLEL DO
      enddo
      
      vorticity = (dvdx - dudy) / mesh%sqrtG(ids:ide,jds:jde,ifs:ife)
      
      tend%phiG = -div_x - div_y
      tend%u    = -flux_x + mesh%sqrtG(ids:ide,jds:jde,:) * (vorticity + mesh%f(ids:ide,jds:jde,:)) * stat%contraV(ids:ide,jds:jde,:)
      tend%v    = -flux_y - mesh%sqrtG(ids:ide,jds:jde,:) * (vorticity + mesh%f(ids:ide,jds:jde,:)) * stat%contraU(ids:ide,jds:jde,:)
      
    end subroutine spatial_operator
    
    ! Input 3 cells' infomation to compute the derivative in the center cell
    subroutine calc_tendP(derivP,f,q,eigenvalueL,eigenvalueR)
      real   , intent(out) :: derivP(DOF)
      real   , intent(in ) :: f     (3*DOF-2)
      real   , intent(in ) :: q     (3*DOF-2)
      real   , intent(in ) :: eigenvalueL
      real   , intent(in ) :: eigenvalueR
      
      real fxL  (2)
      real fxR  (2)
      real qxL  (2)
      real qxR  (2)
      real fCell(DOF,3)
      real qCell(DOF,3)
      
      real fxVIA       ! The 1st order derivative of f on cell
      real fxxc        ! The 2nd order derivative of f
      
      real lambdaL  ! left eigenvalue
      real lambdaR  ! right eigenvalue
      
      integer P1,P2,P3,P4 ! points in cell
      
      integer iCell
      
      do iCell = 1, 3
        fCell(:,iCell) = f(pvIdx(:,iCell))
        qCell(:,iCell) = q(pvIdx(:,iCell))
      enddo
      
      call polyfitR(fxL(1),fCell(:,1),dx)
      call polyfitR(fxL(2),fCell(:,2),dx)
      call polyfitL(fxR(1),fCell(:,2),dx)
      call polyfitL(fxR(2),fCell(:,3),dx)
      
      call polyfitR(qxL(1),qCell(:,1),dx)
      call polyfitR(qxL(2),qCell(:,2),dx)
      call polyfitL(qxR(1),qCell(:,2),dx)
      call polyfitL(qxR(2),qCell(:,3),dx)
      
      call riemann_solver(derivP(1  ),fxL(1),fxR(1),qxL(1),qxR(1),eigenvalueL)
      call riemann_solver(derivP(DOF),fxL(2),fxR(2),qxL(2),qxR(2),eigenvalueR)
      
      ! Compute tend of inner point(s) on cells
#    ifdef MCV3
      ! For MCV3 only
      P1        = pvIdx(1,2)
      P2        = pvIdx(2,2)
      P3        = pvIdx(3,2)
      
      fxVIA     = (f(P3) - f(P1)) / dx
      derivP(2) = 1.5 * fxVIA - 0.25 * (derivP(1) + derivP(DOF))
#    endif
#    ifdef MCV4
      ! For MCV4 only
      P1    = pvIdx(1,2)
      P2    = pvIdx(2,2)
      P3    = pvIdx(3,2)
      P4    = pvIdx(4,2)
      
      fxVIA = (f(P4) - f(P1)) / dx
      fxxc  = 4.5 * ( f(P1) - f(P2) - f(P3) + f(P4) ) / (dx**2)
      
      derivP(2) = 4./3. * fxVIA - (4. * derivP(1) + 5. * derivP(DOF)) / 27. - 4. / 27. * dx * fxxc
      derivP(3) = 4./3. * fxVIA - (5. * derivP(1) + 4. * derivP(DOF)) / 27. + 4. / 27. * dx * fxxc
#    endif
      
    end subroutine calc_tendP
    
    subroutine polyfitL(fitDerivCL,pv,dh)
      real, intent(in ) :: pv(DOF)
      real, intent(in ) :: dh
      real, intent(out) :: fitDerivCL ! left side tend of the cell
      
#    ifdef MCV3
      ! For MCV3 only
      fitDerivCL = -( pv(3) - 4. * pv(2) + 3. * pv(1)) / dh
#    endif

#    ifdef MCV4
      ! For MCV4 only
      fitDerivCL = (-11. * pv(1) + 18. * pv(2) - 9. * pv(3) + 2. * pv(4))/(2. * dh);
#    endif
    end subroutine polyfitL
    
    subroutine polyfitR(fitDerivCR,pv,dh)
      real, intent(in ) :: pv(DOF)
      real, intent(in ) :: dh
      real, intent(out) :: fitDerivCR ! right side tend of the cell
      
#    ifdef MCV3
      ! For MCV3 only
      fitDerivCR =  ( pv(1) - 4. * pv(2) + 3. * pv(3)) / dh
#    endif

#    ifdef MCV4
      ! For MCV4 only
      fitDerivCR = ( 11. * pv(4) - 18. * pv(3) + 9. * pv(2) - 2. * pv(1))/(2. * dh);
#    endif
    end subroutine polyfitR
    
    subroutine riemann_solver(derivP,fxL,fxR,qxL,qxR,eigenvalue)
      real, intent(in ) :: fxL
      real, intent(in ) :: fxR
      real, intent(in ) :: qxL
      real, intent(in ) :: qxR
      real, intent(in ) :: eigenvalue
      real, intent(out) :: derivP
      
      derivP = 0.5 * ((fxL + fxR) - eigenvalue * (qxR - qxL))
    
    end subroutine riemann_solver
    
    function eigenvalue_x(contraU,phi,matrixIG) result(lambda)
      real contraU      ! contravariant wind on x direction
      real phi          ! geopotential height
      real lambda
      real matrixIG(2,2)
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
    
    subroutine unify_bdy_stat(stat)
      type(stat_field), intent(inout) :: stat
      
      call unify_bdy_field(stat%phi            (ids:ide,jds:jde,ifs:ife))
      call unify_bdy_field(stat%zonal_wind     (ids:ide,jds:jde,ifs:ife))
      call unify_bdy_field(stat%meridional_wind(ids:ide,jds:jde,ifs:ife))
      
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
      
      !$OMP PARALLEL DO PRIVATE(i,j)
      do iPatch = ifs, ife
        do j = jps, jpe
          do i = ips, ipe
            call covProjPlane2Sphere(stat%zonal_wind(i,j,iPatch), stat%meridional_wind(i,j,iPatch), stat%u(i,j,iPatch), stat%v(i,j,iPatch), mesh%matrixA(:,:,i,j,iPatch), mesh%matrixIG(:,:,i,j,iPatch))
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO
    end subroutine convert_wind_P2SP
    
    ! convert vector from patch to sphere, for ghost zone
    subroutine convert_ghost_wind_P2SP(stat)
      type(stat_field), intent(inout) :: stat
      
      integer i,j,iPatch
      
      !$OMP PARALLEL DO PRIVATE(i,j)
      do iPatch = ifs, ife
        ! left
        do j = jds, jde
          do i = ids, ids+nPVHalo
            call covProjPlane2Sphere(stat%zonal_wind(i,j,iPatch), stat%meridional_wind(i,j,iPatch), stat%u(i,j,iPatch), stat%v(i,j,iPatch), mesh%matrixA(:,:,i,j,iPatch), mesh%matrixIG(:,:,i,j,iPatch))
          enddo
        enddo
        ! right
        do j = jds, jde
          do i = ide-nPVHalo, ide
            call covProjPlane2Sphere(stat%zonal_wind(i,j,iPatch), stat%meridional_wind(i,j,iPatch), stat%u(i,j,iPatch), stat%v(i,j,iPatch), mesh%matrixA(:,:,i,j,iPatch), mesh%matrixIG(:,:,i,j,iPatch))
          enddo
        enddo
        ! top
        do j = jde-nPVHalo, jde
          do i = ids+nPVHalo, ide-nPVHalo
            call covProjPlane2Sphere(stat%zonal_wind(i,j,iPatch), stat%meridional_wind(i,j,iPatch), stat%u(i,j,iPatch), stat%v(i,j,iPatch), mesh%matrixA(:,:,i,j,iPatch), mesh%matrixIG(:,:,i,j,iPatch))
          enddo
        enddo
        ! bottom
        do j = jds, jds+nPVHalo
          do i = ids+nPVHalo, ide-nPVHalo
            call covProjPlane2Sphere(stat%zonal_wind(i,j,iPatch), stat%meridional_wind(i,j,iPatch), stat%u(i,j,iPatch), stat%v(i,j,iPatch), mesh%matrixA(:,:,i,j,iPatch), mesh%matrixIG(:,:,i,j,iPatch))
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO
    end subroutine convert_ghost_wind_P2SP
    
    subroutine convert_wind_SP2P(stat)
      type(stat_field), intent(inout) :: stat
      
      integer i,j,iPatch
      
      !$OMP PARALLEL DO PRIVATE(i,j)
      do iPatch = ifs, ife
        do j = jps, jpe
          do i = ips, ipe
            call covProjSphere2Plane(stat%u(i,j,iPatch), stat%v(i,j,iPatch), stat%zonal_wind(i,j,iPatch), stat%meridional_wind(i,j,iPatch), mesh%matrixIA(:,:,i,j,iPatch), mesh%matrixG(:,:,i,j,iPatch))
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO
    end subroutine convert_wind_SP2P
    
    subroutine convert_ghost_wind_SP2P(stat)
      type(stat_field), intent(inout) :: stat
      
      integer i,j,iPatch
      
      !$OMP PARALLEL DO PRIVATE(i,j)
      do iPatch = ifs, ife
        ! left
        do j = jds, jde
          do i = ips, ids
            call covProjSphere2Plane(stat%u(i,j,iPatch), stat%v(i,j,iPatch), stat%zonal_wind(i,j,iPatch), stat%meridional_wind(i,j,iPatch), mesh%matrixIA(:,:,i,j,iPatch), mesh%matrixG(:,:,i,j,iPatch))
          enddo
        enddo
        ! right
        do j = jds, jde
          do i = ide, ipe
            call covProjSphere2Plane(stat%u(i,j,iPatch), stat%v(i,j,iPatch), stat%zonal_wind(i,j,iPatch), stat%meridional_wind(i,j,iPatch), mesh%matrixIA(:,:,i,j,iPatch), mesh%matrixG(:,:,i,j,iPatch))
          enddo
        enddo
        ! top
        do j = jde, jpe
          do i = ids, ide
            call covProjSphere2Plane(stat%u(i,j,iPatch), stat%v(i,j,iPatch), stat%zonal_wind(i,j,iPatch), stat%meridional_wind(i,j,iPatch), mesh%matrixIA(:,:,i,j,iPatch), mesh%matrixG(:,:,i,j,iPatch))
          enddo
        enddo
        ! bottom
        do j = jps, jds
          do i = ids, ide
            call covProjSphere2Plane(stat%u(i,j,iPatch), stat%v(i,j,iPatch), stat%zonal_wind(i,j,iPatch), stat%meridional_wind(i,j,iPatch), mesh%matrixIA(:,:,i,j,iPatch), mesh%matrixG(:,:,i,j,iPatch))
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO
    end subroutine convert_ghost_wind_SP2P
    
    subroutine convert_wind_cov2contrav(stat)
      type(stat_field), intent(inout) :: stat
      
      integer i,j,iPatch
      
      !$OMP PARALLEL DO PRIVATE(i,j)
      do iPatch = ifs, ife
        do j = jps, jpe
          do i = ips, ipe
            call cov2contrav(stat%contraU(i,j,iPatch),stat%contraV(i,j,iPatch),stat%u(i,j,iPatch),stat%v(i,j,iPatch),mesh%matrixIG(:,:,i,j,iPatch))
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO
    end subroutine convert_wind_cov2contrav
    
END MODULE spatial_operators_mod

