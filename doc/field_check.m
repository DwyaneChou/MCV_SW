% This program is the test of Saturated Vapor Pressure
% clc
% clear

var_name = 'phiP_t';
it       = 0*24+1;

nc_file = '..\run\mcv_output.nc';

ids        = ncreadatt(nc_file,'/','ids');
ide        = ncreadatt(nc_file,'/','ide');
jds        = ncreadatt(nc_file,'/','jds');
jde        = ncreadatt(nc_file,'/','jde');
ips        = ncreadatt(nc_file,'/','ips');
ipe        = ncreadatt(nc_file,'/','ipe');
jps        = ncreadatt(nc_file,'/','jps');
jpe        = ncreadatt(nc_file,'/','jpe');
Fill_Value = ncreadatt(nc_file,var_name,'_FillValue');

nHalo = ids-ips;
its   = ids-ips+1;
ite   = ipe;
jts   = jds-jps+1;
jte   = jpe;
nPVx  = ite-its+1;
nPVy  = jte-jts+1;

lon = ncread(nc_file,'lonP'  ,[its,jts,1   ],[nPVx,nPVy,6  ]);
lat = ncread(nc_file,'latP'  ,[its,jts,1   ],[nPVx,nPVy,6  ]);
var = ncread(nc_file,var_name,[its,jts,1,it],[nPVx,nPVy,6,1]);

lon(lon<0) = 360 + lon((lon<0));

lon1d = reshape(lon,[],1);
lat1d = reshape(lat,[],1);
var1d = reshape(var,[],1);

res = 1;
x   = 0:res:360;
y   = -90:res:90;

[lon2d,lat2d] = meshgrid(x,y);

var_plot = griddata(lon1d,lat1d,var1d,lon2d,lat2d,'linear');

figure
pcolor(lon2d,lat2d,var_plot)
shading interp
% set(gca,'CLim',[-16,38])
colormap(jet)