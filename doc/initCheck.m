

nc_file = 'E:\Study\Models\MCV\MCV_SW\run\mcv_output.nc';

u  =squeeze(ncread(nc_file,'zonal_wind'));
v  =squeeze(ncread(nc_file,'meridional_wind'));
phi=squeeze(ncread(nc_file,'phiP'));
lat=squeeze(ncread(nc_file,'latP'));

ids=4;
ide=364;
jds=4;
jde=184;
it = 2;

phi1l  =         u(:,jds:94,it);
phi1u  = fliplr(u(:,94:jde,it));
diff   = phi1l-phi1u;
pcolor(diff')
shading interp