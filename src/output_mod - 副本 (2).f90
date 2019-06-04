module output_mod
  use netcdf
  use constants_mod
  use parameters_mod
  use mesh_mod
  use stat_mod
  use projection_mod
  implicit none
    
    character(13) :: ncFile = 'mcv_output.nc'
    
    contains
    subroutine history_init(stat)
      type(stat_field), intent(in) :: stat
      
      integer status
      integer ncid
      integer lonC_dim_id,latC_dim_id
      integer lonP_dim_id,latP_dim_id
      integer patch_dim_id
      integer time_dim_id
      integer lonC_id,latC_id
      integer lonP_id,latP_id
      integer xC_id,yC_id
      integer xP_id,yP_id
      integer time_id
      integer u_id,v_id
      integer zonal_wind_id,meridional_wind_id
      integer phi_id
      
      status = nf90_create(ncFile, NF90_CLOBBER + NF90_64BIT_OFFSET , ncid)
      if(status/=nf90_noerr) call handle_err(status)
      
      status = nf90_def_dim(ncid,'lonC'  ,Nx            ,lonC_dim_id )
      status = nf90_def_dim(ncid,'latC'  ,Ny            ,latC_dim_id )
      status = nf90_def_dim(ncid,'lonP'  ,nPVx          ,lonP_dim_id )
      status = nf90_def_dim(ncid,'latP'  ,nPVy          ,latP_dim_id )
      status = nf90_def_dim(ncid,'nPatch',Nf            ,patch_dim_id)
      status = nf90_def_dim(ncid,'time'  ,NF90_UNLIMITED,time_dim_id )
      if(status/=nf90_noerr) call handle_err(status)
      
      status = nf90_def_var(ncid,'lonC'           ,NF90_DOUBLE,(/lonC_dim_id,latC_dim_id,patch_dim_id            /),lonC_id           )
      status = nf90_def_var(ncid,'latC'           ,NF90_DOUBLE,(/lonC_dim_id,latC_dim_id,patch_dim_id            /),latC_id           )
      status = nf90_def_var(ncid,'lonP'           ,NF90_DOUBLE,(/lonP_dim_id,latP_dim_id,patch_dim_id            /),lonP_id           )
      status = nf90_def_var(ncid,'latP'           ,NF90_DOUBLE,(/lonP_dim_id,latP_dim_id,patch_dim_id            /),latP_id           )
      status = nf90_def_var(ncid,'uP'             ,NF90_DOUBLE,(/lonP_dim_id,latP_dim_id,patch_dim_id,time_dim_id/),u_id              )
      status = nf90_def_var(ncid,'vP'             ,NF90_DOUBLE,(/lonP_dim_id,latP_dim_id,patch_dim_id,time_dim_id/),v_id              )
      status = nf90_def_var(ncid,'zonal_wind'     ,NF90_DOUBLE,(/lonP_dim_id,latP_dim_id,patch_dim_id,time_dim_id/),zonal_wind_id     )
      status = nf90_def_var(ncid,'meridional_wind',NF90_DOUBLE,(/lonP_dim_id,latP_dim_id,patch_dim_id,time_dim_id/),meridional_wind_id)
      status = nf90_def_var(ncid,'phiP'           ,NF90_DOUBLE,(/lonP_dim_id,latP_dim_id,patch_dim_id,time_dim_id/),phi_id            )
      status = nf90_def_var(ncid,'time'           ,NF90_INT   ,(/                                     time_dim_id/),time_id           )
      if(status/=nf90_noerr) call handle_err(status)
      
      !print*,'nf90_put_att'
      status = nf90_put_att(ncid,lonC_id          ,'units'     ,'degree_east' )
      status = nf90_put_att(ncid,latC_id          ,'units'     ,'degree_north')
      status = nf90_put_att(ncid,lonP_id          ,'units'     ,'degree_east' )
      status = nf90_put_att(ncid,latP_id          ,'units'     ,'degree_north')
      status = nf90_put_att(ncid,time_id          ,'units'     ,'seconds'     )
      status = nf90_put_att(ncid,zonal_wind_id    ,'units'     ,'m/s'         )
      status = nf90_put_att(ncid,meridional_wind_id,'units'     ,'m/s'         )
      status = nf90_put_att(ncid,lonC_id          ,'long_name' ,'longitude on sphere coordinate for Cells' )
      status = nf90_put_att(ncid,latC_id          ,'long_name' ,'latitude on sphere coordinate for Cells'  )
      status = nf90_put_att(ncid,lonP_id          ,'long_name' ,'longitude on sphere coordinate for Points')
      status = nf90_put_att(ncid,latP_id          ,'long_name' ,'latitude on sphere coordinate for Points' )
      status = nf90_put_att(ncid,time_id          ,'long_name' ,'time'                                     )
      status = nf90_put_att(ncid,zonal_wind_id    ,'long_name' ,'zonal wind'                               )
      status = nf90_put_att(ncid,meridional_wind_id,'long_name','meridional_wind'                          )
      !if(status/=nf90_noerr) call handle_err(status)
      
      status = nf90_enddef(ncid)
      if(status/=nf90_noerr) call handle_err(status)
      
      status = nf90_put_var(ncid,lonC_id, mesh%lonC(1:Nx  ,1:Ny  ,:) * R2D)
      status = nf90_put_var(ncid,latC_id, mesh%latC(1:Nx  ,1:Ny  ,:) * R2D)
      status = nf90_put_var(ncid,lonP_id, mesh%lonP(1:nPVx,1:nPVy,:) * R2D)
      status = nf90_put_var(ncid,latP_id, mesh%latP(1:nPVx,1:nPVy,:) * R2D)
      if(status/=nf90_noerr) call handle_err(status)
      
      status = nf90_close(ncid)
      if(status/=nf90_noerr) call handle_err(status)
      
    end subroutine history_init
    
    subroutine history_write_stat(stat,time_slot_num)
      type(stat_field), intent(in) :: stat
      integer         , intent(in) :: time_slot_num
      
      integer status
      integer ncid
      integer time_id
      integer u_id,v_id
      integer zonal_wind_id,meridional_wind_id
      integer phi_id
      
      integer :: time(1)
      
      real, dimension(ids:ide,jds:jde,ifs:ife) :: u ! zonal wind
      real, dimension(ids:ide,jds:jde,ifs:ife) :: v ! merdional wind
      
      real contraU
      real contraV
      
      integer i,j,iPatch
      
      time(1) = time_slot_num
      
      do iPatch = ifs, ife
        do j = jds, jde
          do i = ids, ide
            call cov2contrav            (contraU      , contraV      , stat%u(i,j,iPatch), stat%v(i,j,iPatch), mesh%matrixIG(:,:,i,j,iPatch))
            call contravProjPlane2Sphere(u(i,j,iPatch), v(i,j,iPatch), contraU           , contraV           , mesh%matrixA (:,:,i,j,iPatch))
          enddo
        enddo
      enddo
      
      !print*,'nf90_open'
      status = nf90_open(ncFile,NF90_WRITE,ncid)
      if(status/=nf90_noerr) call handle_err(status)
      
      !print*,'nf90_inq_varid'
      status = nf90_inq_varid(ncid,'time'           , time_id           )
      status = nf90_inq_varid(ncid,'uP'             , u_id              )
      status = nf90_inq_varid(ncid,'vP'             , v_id              )
      status = nf90_inq_varid(ncid,'zonal_wind'     , zonal_wind_id     )
      status = nf90_inq_varid(ncid,'meridional_wind', meridional_wind_id)
      status = nf90_inq_varid(ncid,'phiP'           , phi_id            )
      if(status/=nf90_noerr) call handle_err(status)
      
      !print*,'nf90_put_var'
      status = nf90_put_var(ncid,time_id           ,time                             ,start=(/      time_slot_num/),count=(/             1/))
      status = nf90_put_var(ncid,u_id              ,stat%u  (ids:ide,jds:jde,ifs:ife),start=(/1,1,1,time_slot_num/),count=(/nPVx,nPVy,Nf,1/))
      status = nf90_put_var(ncid,v_id              ,stat%v  (ids:ide,jds:jde,ifs:ife),start=(/1,1,1,time_slot_num/),count=(/nPVx,nPVy,Nf,1/))
      status = nf90_put_var(ncid,zonal_wind_id     ,u       (ids:ide,jds:jde,ifs:ife),start=(/1,1,1,time_slot_num/),count=(/nPVx,nPVy,Nf,1/))
      status = nf90_put_var(ncid,meridional_wind_id,v       (ids:ide,jds:jde,ifs:ife),start=(/1,1,1,time_slot_num/),count=(/nPVx,nPVy,Nf,1/))
      status = nf90_put_var(ncid,phi_id            ,stat%phi(ids:ide,jds:jde,ifs:ife),start=(/1,1,1,time_slot_num/),count=(/nPVx,nPVy,Nf,1/))
      if(status/=nf90_noerr) call handle_err(status)
      
      !print*,'nf90_close'
      status = nf90_close(ncid)
      if(status/=nf90_noerr) call handle_err(status)
      
    end subroutine history_write_stat
    
    subroutine handle_err(status)
      implicit none
      integer,intent(in)::status
            
      if(status/=nf90_noerr)then
          print*, trim(nf90_strerror(status))
          stop "Stopped by netCDF"
      endif  
    endsubroutine handle_err
end module output_mod
    