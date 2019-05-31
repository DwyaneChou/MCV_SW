module output_mod
  use parameters_mod
  use mesh_mod
  use stat_mod
  use projection_mod
  use io_mod
  implicit none
  
    contains
    subroutine write_netCDF
      character(13) ncFile
    
      real, dimension(ips:ipe,jps:jpe,ifs:ife) :: u ! zonal wind
      real, dimension(ips:ipe,jps:jpe,ifs:ife) :: v ! merdional wind
      
      real contraU
      real contraV
      
      integer i,j,iPatch
      
      do iPatch = ifs, ife
        do j = jps, jpe
          do i = ips, ipe
            call cov2contrav            (contraU      , contraV      , stat(0)%uP(i,j,iPatch), stat(0)%vP(i,j,iPatch), mesh%matrixIG(:,:,i,j,iPatch))
            call contravProjPlane2Sphere(u(i,j,iPatch), v(i,j,iPatch), contraU               , contraV               , mesh%matrixA (:,:,i,j,iPatch))
          enddo
        enddo
      enddo
    
      ncFile = 'mcv_output.nc'
      
      call io_init
      
      ! output data into uniform netCDF file
      call io_create_dataset( name='mcv_output', file_path=ncFile, mode='output' )
      
      call io_add_dim(name='lonC'    , dataset_name='mcv_output', size=Nx_halo  )
      call io_add_dim(name='latC'    , dataset_name='mcv_output', size=Ny_halo  )
      call io_add_dim(name='lonP'    , dataset_name='mcv_output', size=nPVx_halo)
      call io_add_dim(name='latP'    , dataset_name='mcv_output', size=nPVy_halo)
      call io_add_dim(name='nPatch'  , dataset_name='mcv_output', size=Nf       )
      
      ! attributes
      call io_add_att(name='DOF'  , value=DOF         , dataset_name='mcv_output')
      call io_add_att(name='x_min', value=x_min * R2D , dataset_name='mcv_output')
      call io_add_att(name='x_max', value=x_max * R2D , dataset_name='mcv_output')
      call io_add_att(name='y_min', value=y_min * R2D , dataset_name='mcv_output')
      call io_add_att(name='y_max', value=y_max * R2D , dataset_name='mcv_output')
      call io_add_att(name='dx'   , value=dx    * R2D , dataset_name='mcv_output')
      call io_add_att(name='dy'   , value=dy    * R2D , dataset_name='mcv_output')
      call io_add_att(name='Nx'   , value=Nx          , dataset_name='mcv_output')
      call io_add_att(name='Ny'   , value=Ny          , dataset_name='mcv_output')
      call io_add_att(name='Nf'   , value=Nf          , dataset_name='mcv_output')
      call io_add_att(name='nPVx' , value=nPVx        , dataset_name='mcv_output')
      call io_add_att(name='nPVy' , value=nPVy        , dataset_name='mcv_output')
      call io_add_att(name='xhalo', value=xhalo       , dataset_name='mcv_output')
      call io_add_att(name='yhalo', value=yhalo       , dataset_name='mcv_output')
      call io_add_att(name='ids'  , value=ids         , dataset_name='mcv_output')
      call io_add_att(name='ide'  , value=ide         , dataset_name='mcv_output')
      call io_add_att(name='jds'  , value=jds         , dataset_name='mcv_output')
      call io_add_att(name='jde'  , value=jde         , dataset_name='mcv_output')
      call io_add_att(name='ics'  , value=ics         , dataset_name='mcv_output')
      call io_add_att(name='ice'  , value=ice         , dataset_name='mcv_output')
      call io_add_att(name='jcs'  , value=jcs         , dataset_name='mcv_output')
      call io_add_att(name='jce'  , value=jce         , dataset_name='mcv_output')
      call io_add_att(name='ips'  , value=ips         , dataset_name='mcv_output')
      call io_add_att(name='ipe'  , value=ipe         , dataset_name='mcv_output')
      call io_add_att(name='jps'  , value=jps         , dataset_name='mcv_output')
      call io_add_att(name='jpe'  , value=jpe         , dataset_name='mcv_output')
      call io_add_att(name='ifs'  , value=ifs         , dataset_name='mcv_output')
      call io_add_att(name='ife'  , value=ife         , dataset_name='mcv_output')
      
      ! variables
      call io_add_var(name='xC'  , dataset_name='mcv_output', long_name='central angle on x direction of cell center', units='radian'      , dim_names=['lonC  ', 'latC  ', 'nPatch'], data_type='real(8)')
      call io_add_var(name='yC'  , dataset_name='mcv_output', long_name='central angle on y direction of cell center', units='radian'      , dim_names=['lonC  ', 'latC  ', 'nPatch'], data_type='real(8)')
      call io_add_var(name='lonC', dataset_name='mcv_output', long_name='longitude on sphere coordinate for Cells'   , units='degree_east' , dim_names=['lonC  ', 'latC  ', 'nPatch'], data_type='real(8)')
      call io_add_var(name='latC', dataset_name='mcv_output', long_name='latitude on sphere coordinate for Cells'    , units='degree_north', dim_names=['lonC  ', 'latC  ', 'nPatch'], data_type='real(8)')
      call io_add_var(name='xP'  , dataset_name='mcv_output', long_name='central angle on x direction of Points'     , units='radian'      , dim_names=['lonP  ', 'latP  ', 'nPatch'], data_type='real(8)')
      call io_add_var(name='yP'  , dataset_name='mcv_output', long_name='central angle on y direction of Points'     , units='radian'      , dim_names=['lonP  ', 'latP  ', 'nPatch'], data_type='real(8)')
      call io_add_var(name='lonP', dataset_name='mcv_output', long_name='longitude on sphere coordinate for Points'  , units='degree_east' , dim_names=['lonP  ', 'latP  ', 'nPatch'], data_type='real(8)')
      call io_add_var(name='latP', dataset_name='mcv_output', long_name='latitude on sphere coordinate for Points'   , units='degree_north', dim_names=['lonP  ', 'latP  ', 'nPatch'], data_type='real(8)')
      
      call io_add_var(name='uP'  , dataset_name='mcv_output', long_name='covariant u wind on cube points'            , units='m/s'         , dim_names=['lonP  ', 'latP  ', 'nPatch'], data_type='real(8)')
      call io_add_var(name='vP'  , dataset_name='mcv_output', long_name='covariant v wind on cube points'            , units='m/s'         , dim_names=['lonP  ', 'latP  ', 'nPatch'], data_type='real(8)')
      call io_add_var(name='phiP', dataset_name='mcv_output', long_name='geopotential height on cube points'         , units='gpm'         , dim_names=['lonP  ', 'latP  ', 'nPatch'], data_type='real(8)')
                                                                                                                                           
      call io_add_var(name='uC'  , dataset_name='mcv_output', long_name='covariant u wind on cube cells'             , units='m/s'         , dim_names=['lonC  ', 'latC  ', 'nPatch'], data_type='real(8)')
      call io_add_var(name='vC'  , dataset_name='mcv_output', long_name='covariant v wind on cube cells'             , units='m/s'         , dim_names=['lonC  ', 'latC  ', 'nPatch'], data_type='real(8)')
      call io_add_var(name='phiC', dataset_name='mcv_output', long_name='geopotential height on cube cells'          , units='gpm'         , dim_names=['lonC  ', 'latC  ', 'nPatch'], data_type='real(8)')
                                                                                                                                           
      call io_add_var(name='u'   , dataset_name='mcv_output', long_name='zonal wind'                                 , units='m/s'         , dim_names=['lonP  ', 'latP  ', 'nPatch'], data_type='real(8)')
      call io_add_var(name='v'   , dataset_name='mcv_output', long_name='merdional wind'                             , units='m/s'         , dim_names=['lonP  ', 'latP  ', 'nPatch'], data_type='real(8)')
      
      call io_start_output(dataset_name='mcv_output')
      
      call io_output(name='xC'          , array=mesh%xC         ,dataset_name='mcv_output')
      call io_output(name='yC'          , array=mesh%yC         ,dataset_name='mcv_output')
      call io_output(name='xP'          , array=mesh%xP         ,dataset_name='mcv_output')
      call io_output(name='yP'          , array=mesh%yP         ,dataset_name='mcv_output')
      call io_output(name='lonC'        , array=mesh%lonC * R2D ,dataset_name='mcv_output')
      call io_output(name='latC'        , array=mesh%latC * R2D ,dataset_name='mcv_output')
      call io_output(name='lonP'        , array=mesh%lonP * R2D ,dataset_name='mcv_output')
      call io_output(name='latP'        , array=mesh%latP * R2D ,dataset_name='mcv_output')
      
      call io_output(name='uP'          , array=stat(0)%uP   ,dataset_name='mcv_output')
      call io_output(name='vP'          , array=stat(0)%vP   ,dataset_name='mcv_output')
      call io_output(name='phiP'        , array=stat(0)%phiP ,dataset_name='mcv_output')
      !call io_output(name='uC'          , array=stat(0)%uC   ,dataset_name='mcv_output')
      !call io_output(name='vC'          , array=stat(0)%vC   ,dataset_name='mcv_output')
      !call io_output(name='phiC'        , array=stat(0)%phiC ,dataset_name='mcv_output')
      call io_output(name='u'           , array=u            ,dataset_name='mcv_output')
      call io_output(name='v'           , array=v            ,dataset_name='mcv_output')
      
      call io_end_output('mcv_output')
    end subroutine write_netCDF
    
end module output_mod
    