clear;
%% load ascii data

% t0 = importdata('advectGaussHills_triUnif_3__0000_Interp_tracer1.txt');
% t1 = importdata('advectGaussHills_triUnif_3__0001_Interp_tracer1.txt');
% t2 = importdata('advectGaussHills_triUnif_3__0002_Interp_tracer1.txt');

t0 = importdata('advectGaussHills_triUnif_2__0000_Interp_tracer1.txt');
t1 = importdata('advectGaussHills_triUnif_2__0001_Interp_tracer1.txt');
t2 = importdata('advectGaussHills_triUnif_2__0002_Interp_tracer1.txt');

outputFile = 'gaussHills_tri2.nc';

%% set up coordinates

nLon = 720;
nLat = 361;
nT = 3;
lon = 0:0.5:359.5;
lat = -90:0.5:90;

ncid = netcdf.create(outputFile,'NC_WRITE');

%% define coordinates

dimid_lon = netcdf.defDim(ncid,'longitude',nLon);
dimid_lat = netcdf.defDim(ncid,'latitude',nLat);

varid_lon = netcdf.defVar(ncid,'lon','double',dimid_lon);
netcdf.putAtt(ncid,varid_lon,'long_name','Longitude');
netcdf.putAtt(ncid,varid_lon,'units','degrees_east');

varid_lat = netcdf.defVar(ncid,'lat','double',dimid_lat);
netcdf.putAtt(ncid,varid_lat,'long_name','Latitude');
netcdf.putAtt(ncid,varid_lat,'units','degrees_north');

%% define variables

varid_tracer0 = netcdf.defVar(ncid,'tracer0','double',[dimid_lon,dimid_lat]);
varid_tracer1 = netcdf.defVar(ncid,'tracer1','double',[dimid_lon,dimid_lat]);
varid_tracer2 = netcdf.defVar(ncid,'tracer2','double',[dimid_lon,dimid_lat]);

netcdf.endDef(ncid);

%% write variables and coordinates

netcdf.putVar(ncid, varid_lon, lon);
netcdf.putVar(ncid, varid_lat, lat);
netcdf.putVar(ncid, varid_tracer0, t0');
netcdf.putVar(ncid, varid_tracer1, t1');
netcdf.putVar(ncid, varid_tracer2, t2');


netcdf.close(ncid);