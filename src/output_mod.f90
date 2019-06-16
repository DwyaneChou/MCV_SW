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
      integer areaCell_id
      integer time_id
      integer u_id,v_id
      integer phi_id
      integer zonal_wind_id,meridional_wind_id
      integer uC_id,vC_id
      integer phiC_id
      
      status = nf90_create(ncFile, NF90_CLOBBER + NF90_64BIT_OFFSET , ncid)
      if(status/=nf90_noerr) call handle_err(status)
      
      status = nf90_def_dim(ncid,'lonC'  ,Nx_halo       ,lonC_dim_id )
      status = nf90_def_dim(ncid,'latC'  ,Ny_halo       ,latC_dim_id )
      status = nf90_def_dim(ncid,'lonP'  ,nPVx_halo     ,lonP_dim_id )
      status = nf90_def_dim(ncid,'latP'  ,nPVy_halo     ,latP_dim_id )
      status = nf90_def_dim(ncid,'nPatch',Nf            ,patch_dim_id)
      status = nf90_def_dim(ncid,'time'  ,NF90_UNLIMITED,time_dim_id )
      if(status/=nf90_noerr) call handle_err(status)
      
      status = nf90_def_var(ncid,'lonC'           ,NF90_DOUBLE,(/lonC_dim_id,latC_dim_id,patch_dim_id            /),lonC_id           )
      status = nf90_def_var(ncid,'latC'           ,NF90_DOUBLE,(/lonC_dim_id,latC_dim_id,patch_dim_id            /),latC_id           )
      status = nf90_def_var(ncid,'lonP'           ,NF90_DOUBLE,(/lonP_dim_id,latP_dim_id,patch_dim_id            /),lonP_id           )
      status = nf90_def_var(ncid,'latP'           ,NF90_DOUBLE,(/lonP_dim_id,latP_dim_id,patch_dim_id            /),latP_id           )
      status = nf90_def_var(ncid,'areaCell'       ,NF90_DOUBLE,(/lonC_dim_id,latC_dim_id,patch_dim_id            /),areaCell_id       )
      status = nf90_def_var(ncid,'time'           ,NF90_INT   ,(/                                     time_dim_id/),time_id           )
      status = nf90_def_var(ncid,'uP'             ,NF90_DOUBLE,(/lonP_dim_id,latP_dim_id,patch_dim_id,time_dim_id/),u_id              )
      status = nf90_def_var(ncid,'vP'             ,NF90_DOUBLE,(/lonP_dim_id,latP_dim_id,patch_dim_id,time_dim_id/),v_id              )
      status = nf90_def_var(ncid,'phiP'           ,NF90_DOUBLE,(/lonP_dim_id,latP_dim_id,patch_dim_id,time_dim_id/),phi_id            )
      status = nf90_def_var(ncid,'zonal_wind'     ,NF90_DOUBLE,(/lonP_dim_id,latP_dim_id,patch_dim_id,time_dim_id/),zonal_wind_id     )
      status = nf90_def_var(ncid,'meridional_wind',NF90_DOUBLE,(/lonP_dim_id,latP_dim_id,patch_dim_id,time_dim_id/),meridional_wind_id)
      status = nf90_def_var(ncid,'uC'             ,NF90_DOUBLE,(/lonC_dim_id,latC_dim_id,patch_dim_id,time_dim_id/),uC_id             )
      status = nf90_def_var(ncid,'vC'             ,NF90_DOUBLE,(/lonC_dim_id,latC_dim_id,patch_dim_id,time_dim_id/),vC_id             )
      status = nf90_def_var(ncid,'phiC'           ,NF90_DOUBLE,(/lonC_dim_id,latC_dim_id,patch_dim_id,time_dim_id/),phiC_id           )
      if(status/=nf90_noerr) call handle_err(status)
      
      !print*,'nf90_put_att'
      status = nf90_put_att(ncid,nf90_global       ,'MCV_ORDER',DOF)
      status = nf90_put_att(ncid,nf90_global       ,'dx'       ,dx*R2D)
      status = nf90_put_att(ncid,nf90_global       ,'dy'       ,dy*R2D)
      status = nf90_put_att(ncid,nf90_global       ,'dt'       ,dt)
      status = nf90_put_att(ncid,nf90_global       ,'xhalo'    ,xhalo)
      status = nf90_put_att(ncid,nf90_global       ,'yhalo'    ,yhalo)
      status = nf90_put_att(ncid,nf90_global       ,'case_num' ,case_num)
      status = nf90_put_att(ncid,nf90_global       ,'ics'      ,ics)
      status = nf90_put_att(ncid,nf90_global       ,'ice'      ,ice)
      status = nf90_put_att(ncid,nf90_global       ,'jcs'      ,jcs)
      status = nf90_put_att(ncid,nf90_global       ,'jce'      ,jce)
      status = nf90_put_att(ncid,nf90_global       ,'ips'      ,ips)
      status = nf90_put_att(ncid,nf90_global       ,'ipe'      ,ipe)
      status = nf90_put_att(ncid,nf90_global       ,'jps'      ,jps)
      status = nf90_put_att(ncid,nf90_global       ,'jpe'      ,jpe)
      status = nf90_put_att(ncid,nf90_global       ,'ids'      ,ids)
      status = nf90_put_att(ncid,nf90_global       ,'ide'      ,ide)
      status = nf90_put_att(ncid,nf90_global       ,'jds'      ,jds)
      status = nf90_put_att(ncid,nf90_global       ,'jde'      ,jde)
      status = nf90_put_att(ncid,nf90_global       ,'ifs'      ,ifs)
      status = nf90_put_att(ncid,nf90_global       ,'ife'      ,ife)
      
      status = nf90_put_att(ncid,lonC_id           ,'units'    ,'degree_east' )
      status = nf90_put_att(ncid,latC_id           ,'units'    ,'degree_north')
      status = nf90_put_att(ncid,lonP_id           ,'units'    ,'degree_east' )
      status = nf90_put_att(ncid,latP_id           ,'units'    ,'degree_north')
      status = nf90_put_att(ncid,areaCell_id       ,'units'    ,'m^2'         )
      status = nf90_put_att(ncid,time_id           ,'units'    ,'seconds'     )
      status = nf90_put_att(ncid,u_id              ,'units'    ,'m/s'         )
      status = nf90_put_att(ncid,v_id              ,'units'    ,'m/s'         )
      status = nf90_put_att(ncid,phi_id            ,'units'    ,'m^2/s^2'     )
      status = nf90_put_att(ncid,zonal_wind_id     ,'units'    ,'m/s'         )
      status = nf90_put_att(ncid,meridional_wind_id,'units'    ,'m/s'         )
      status = nf90_put_att(ncid,uC_id             ,'units'    ,'m/s'         )
      status = nf90_put_att(ncid,vC_id             ,'units'    ,'m/s'         )
      status = nf90_put_att(ncid,phiC_id           ,'units'    ,'m^2/s^2'     )
      
      status = nf90_put_att(ncid,lonC_id           ,'long_name','longitude on sphere coordinate for Cells' )
      status = nf90_put_att(ncid,latC_id           ,'long_name','latitude on sphere coordinate for Cells'  )
      status = nf90_put_att(ncid,lonP_id           ,'long_name','longitude on sphere coordinate for Points')
      status = nf90_put_att(ncid,latP_id           ,'long_name','latitude on sphere coordinate for Points' )
      status = nf90_put_att(ncid,areaCell_id       ,'long_name','area of cells'                            )
      status = nf90_put_att(ncid,time_id           ,'long_name','time'                                     )
      status = nf90_put_att(ncid,u_id              ,'long_name','covariant wind vector on x direcvtion'    )
      status = nf90_put_att(ncid,v_id              ,'long_name','covariant wind vector on y direcvtion'    )
      status = nf90_put_att(ncid,phi_id            ,'long_name','geopotential height on points'            )
      status = nf90_put_att(ncid,zonal_wind_id     ,'long_name','zonal wind'                               )
      status = nf90_put_att(ncid,meridional_wind_id,'long_name','meridional wind'                          )
      status = nf90_put_att(ncid,uC_id             ,'long_name','covariant wind VIA on x direcvtion'       )
      status = nf90_put_att(ncid,vC_id             ,'long_name','covariant wind VIA on y direcvtion'       )
      status = nf90_put_att(ncid,phiC_id           ,'long_name','geopotential height VIA'                  )
      
      status = nf90_put_att(ncid,u_id              ,'_FillValue',FillValue)
      status = nf90_put_att(ncid,v_id              ,'_FillValue',FillValue)
      status = nf90_put_att(ncid,phi_id            ,'_FillValue',FillValue)
      status = nf90_put_att(ncid,zonal_wind_id     ,'_FillValue',FillValue)
      status = nf90_put_att(ncid,meridional_wind_id,'_FillValue',FillValue)
      status = nf90_put_att(ncid,uC_id             ,'_FillValue',FillValue)
      status = nf90_put_att(ncid,vC_id             ,'_FillValue',FillValue)
      status = nf90_put_att(ncid,phiC_id           ,'_FillValue',FillValue)
      if(status/=nf90_noerr) call handle_err(status)
      
      status = nf90_enddef(ncid)
      if(status/=nf90_noerr) call handle_err(status)
      
      status = nf90_put_var(ncid,lonC_id    , mesh%lonC * R2D)
      status = nf90_put_var(ncid,latC_id    , mesh%latC * R2D)
      status = nf90_put_var(ncid,lonP_id    , mesh%lonP * R2D)
      status = nf90_put_var(ncid,latP_id    , mesh%latP * R2D)
      status = nf90_put_var(ncid,areaCell_id, mesh%areaCell  )
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
      integer phi_id
      integer zonal_wind_id,meridional_wind_id
      integer uC_id,vC_id
      integer phiC_id
      
      integer :: time(1)
      
      integer i,j,iPatch
      
      time(1) = time_slot_num
      !print*,'nf90_open'
      status = nf90_open(ncFile,NF90_WRITE,ncid)
      if(status/=nf90_noerr) call handle_err(status)
      
      !print*,'nf90_inq_varid'
      status = nf90_inq_varid(ncid,'time'           , time_id           )
      status = nf90_inq_varid(ncid,'uP'             , u_id              )
      status = nf90_inq_varid(ncid,'vP'             , v_id              )
      status = nf90_inq_varid(ncid,'phiP'           , phi_id            )
      status = nf90_inq_varid(ncid,'zonal_wind'     , zonal_wind_id     )
      status = nf90_inq_varid(ncid,'meridional_wind', meridional_wind_id)
      status = nf90_inq_varid(ncid,'uC'             , uC_id             )
      status = nf90_inq_varid(ncid,'vC'             , vC_id             )
      status = nf90_inq_varid(ncid,'phiC'           , phiC_id           )
      if(status/=nf90_noerr) call handle_err(status)
      
      !print*,'nf90_put_var'
      status = nf90_put_var(ncid,time_id           ,time                ,start=(/      time_slot_num/),count=(/                       1/))
      status = nf90_put_var(ncid,u_id              ,stat%u              ,start=(/1,1,1,time_slot_num/),count=(/nPVx_halo,nPVy_halo,Nf,1/))
      status = nf90_put_var(ncid,v_id              ,stat%v              ,start=(/1,1,1,time_slot_num/),count=(/nPVx_halo,nPVy_halo,Nf,1/))
      status = nf90_put_var(ncid,phi_id            ,stat%phi            ,start=(/1,1,1,time_slot_num/),count=(/nPVx_halo,nPVy_halo,Nf,1/))
      status = nf90_put_var(ncid,zonal_wind_id     ,stat%zonal_wind     ,start=(/1,1,1,time_slot_num/),count=(/nPVx_halo,nPVy_halo,Nf,1/))
      status = nf90_put_var(ncid,meridional_wind_id,stat%meridional_wind,start=(/1,1,1,time_slot_num/),count=(/nPVx_halo,nPVy_halo,Nf,1/))
      status = nf90_put_var(ncid,uC_id             ,stat%uC             ,start=(/1,1,1,time_slot_num/),count=(/Nx_halo  ,Ny_halo  ,Nf,1/))
      status = nf90_put_var(ncid,vC_id             ,stat%vC             ,start=(/1,1,1,time_slot_num/),count=(/Nx_halo  ,Ny_halo  ,Nf,1/))
      status = nf90_put_var(ncid,phiC_id           ,stat%phiC           ,start=(/1,1,1,time_slot_num/),count=(/Nx_halo  ,Ny_halo  ,Nf,1/))
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
    
    ! add FillValue for output
    subroutine addFillValue(stat)
      type(stat_field), intent(inout) :: stat
      
      ! low left corner
      stat%u(ips:ids-1,jps:jds-1,:) = FillValue
      stat%v(ips:ids-1,jps:jds-1,:) = FillValue
      
      ! low right corner
      stat%u(ide+1:ipe,jps:jds-1,:) = FillValue
      stat%v(ide+1:ipe,jps:jds-1,:) = FillValue
      
      ! up left corner
      stat%u(ips:ids-1,jde+1:jpe,:) = FillValue
      stat%v(ips:ids-1,jde+1:jpe,:) = FillValue
      
      ! up right corner
      stat%u(ide+1:ipe,jde+1:jpe,:) = FillValue
      stat%v(ide+1:ipe,jde+1:jpe,:) = FillValue
    end subroutine addFillValue
end module output_mod
    