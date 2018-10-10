% This function does the searching for locations
function [i,j]=findloc(lon,lat)
load('C:\research\paper_conference\data_matlab\nc_lon.mat');
load('C:\research\paper_conference\data_matlab\nc_lat.mat');
nc_lon=nc_lon-360;

for ii=1:length(nc_lon)
    if nc_lon(ii)>lon
        
        break
    end
end
i=ii-1;


for jj=1:length(nc_lat)
    if nc_lat(jj)>lat
        break
    end
end
j=jj-1;


