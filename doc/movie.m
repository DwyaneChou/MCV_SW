clc
clear

picture_output = 'E:\Study\Models\MCV\MCV_SW\doc\picture';
nc_file = 'E:\Study\Models\MCV\MCV_SW\run\mcv_output.nc';

nt      = 600;
varname = 'phiP_t';
var     = ncread(nc_file,varname,[1 1 1 1],[Inf Inf Inf nt]);
lonP    = ncread(nc_file,'lonP');
latP    = ncread(nc_file,'latP');

ids = ncreadatt(nc_file,'/','ids');
ide = ncreadatt(nc_file,'/','ide');
jds = ncreadatt(nc_file,'/','jds');
jde = ncreadatt(nc_file,'/','jde');
ips = ncreadatt(nc_file,'/','ips');
ipe = ncreadatt(nc_file,'/','ipe');
jps = ncreadatt(nc_file,'/','jps');
jpe = ncreadatt(nc_file,'/','jpe');
ifs = ncreadatt(nc_file,'/','ifs');
ife = ncreadatt(nc_file,'/','ife');

ims=ids-ips+1;
ime=ipe;
jms=jds-jps+1;
jme=jpe;

resPlot = 1;
lon1d   = 0:resPlot:360;
lat1d   = -90:resPlot:90;

[lon2d,lat2d] = meshgrid(lon1d,lat1d);

lonP1d = reshape(lonP(ims:ime,jms:jme,ifs:ife),[],1);
latP1d = reshape(latP(ims:ime,jms:jme,ifs:ife),[],1);

lonP1d(lonP1d<0) = 360+lonP1d(lonP1d<0);
    
for it = 1:nt
    disp(['Plotting fig# ',num2str(it),'/',num2str(nt)])
    
    var1d  = reshape(var(ims:ime,jms:jme,ifs:ife,it),[],1);

    method = 'linear';
    [lon_ext,lat_ext,var_ext] = extend_field(lonP1d,latP1d,var1d);
    F = scatteredInterpolant(double(lon_ext),double(lat_ext),double(var_ext),method,method);
    var_plot = F(lon2d,lat2d);
    
    figure('visible','off')
    pcolor(lon2d,lat2d,var_plot)
    shading interp
    set(gca,'Clim',[49000,59000])
    colormap('jet')
    colorbar
    title(['Geopotential height at hour',num2str(it)])
    
    print('-dpng','-opengl','-r300',[picture_output,'\','MCV4_',num2str(it),'.png']);
end

function [lon_ext,lat_ext,var_ext] = extend_field(lon_src,lat_src,var_src)

n = size(var_src,1);

lon_ext            = lon_src;
lon_ext(n+1  :2*n) = lon_src - 360;
lon_ext(2*n+1:3*n) = lon_src + 360;

lat_ext            = lat_src;
lat_ext(n+1  :2*n) = lat_src;
lat_ext(2*n+1:3*n) = lat_src;

var_ext            = var_src;
var_ext(n+1  :2*n) = var_src;
var_ext(2*n+1:3*n) = var_src;

end